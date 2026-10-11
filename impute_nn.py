import math
import typing as T
import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.utils.data import Dataset, DataLoader

# --------------------------
# Utilities
# --------------------------
def _select_numeric(df: pd.DataFrame) -> T.Tuple[pd.DataFrame, pd.Index]:
    num_cols = df.select_dtypes(include=[np.number]).columns
    return df[num_cols].copy(), num_cols

def _standardize(df_num: pd.DataFrame) -> T.Tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    means = df_num.mean(skipna=True).to_numpy(dtype=np.float32)
    stds = df_num.std(skipna=True).replace(0, 1.0).to_numpy(dtype=np.float32)
    arr = df_num.to_numpy(dtype=np.float32)

    # Make a copy and fill only observed entries
    arr_std = np.zeros_like(arr)
    for j in range(arr.shape[1]):
        col = arr[:, j]
        obs = ~np.isnan(col)
        arr_std[obs, j] = (col[obs] - means[j]) / stds[j]
        arr_std[~obs, j] = np.nan

    return pd.DataFrame(arr_std, columns=df_num.columns, index=df_num.index), means, stds

def _destandardize(arr: np.ndarray, means: np.ndarray, stds: np.ndarray) -> np.ndarray:
    return arr * stds[np.newaxis, :] + means[np.newaxis, :]

# --------------------------
# Dataset
# --------------------------
class TabularMaskDataset(Dataset):
    def __init__(self, X: np.ndarray, mask_prob: float = 0.15, for_remasker: bool = True):
        self.X = X.astype(np.float32)
        self.obs_mask = (~np.isnan(self.X)).astype(np.float32)
        self.mask_prob = mask_prob
        self.for_remasker = for_remasker
        self.X_in = np.nan_to_num(self.X, nan=0.0)

    def __len__(self):
        return self.X.shape[0]

    def __getitem__(self, idx):
        x = self.X_in[idx]
        obs = self.obs_mask[idx]
        if self.for_remasker:
            addl_mask = (np.random.rand(*obs.shape) < self.mask_prob).astype(np.float32) * obs
        else:
            addl_mask = np.zeros_like(obs, dtype=np.float32)
        return (
            torch.tensor(x, dtype=torch.float32),
            torch.tensor(obs, dtype=torch.float32),
            torch.tensor(addl_mask, dtype=torch.float32),
        )

# --------------------------
# Models
# --------------------------
class Autoencoder(nn.Module):
    def __init__(self, d_in: int, hidden: T.Tuple[int, ...] = (128, 64, 32)):
        super().__init__()
        layers = []
        prev = d_in
        for h in hidden:
            layers += [nn.Linear(prev, h), nn.ReLU()]
            prev = h
        dec_layers = []
        prev = hidden[-1]
        for h in reversed(hidden[:-1]):
            dec_layers += [nn.Linear(prev, h), nn.ReLU()]
            prev = h
        dec_layers += [nn.Linear(prev, d_in)]
        self.encoder = nn.Sequential(*layers)
        self.decoder = nn.Sequential(*dec_layers)

    def forward(self, x):
        z = self.encoder(x)
        out = self.decoder(z)
        return out

class RemaskerTransformer(nn.Module):
    def __init__(self, d_in: int, d_model: int = 128, nhead: int = 8, num_layers: int = 4, dim_feedforward: int = 256):
        super().__init__()
        self.d_in = d_in
        self.value_embed = nn.Linear(1, d_model)
        self.pos_embed = nn.Parameter(torch.randn(1, d_in, d_model) * 0.02)
        encoder_layer = nn.TransformerEncoderLayer(d_model=d_model, nhead=nhead,
                                                   dim_feedforward=dim_feedforward, batch_first=True)
        self.encoder = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)
        self.head = nn.Linear(d_model, 1)

    def forward(self, x):
        B, D = x.shape
        x_tok = x.view(B, D, 1)
        h = self.value_embed(x_tok) + self.pos_embed[:, :D, :]
        h = self.encoder(h)
        out = self.head(h).squeeze(-1)
        return out

# --------------------------
# Training & Imputation
# --------------------------
def _masked_mse(pred, target, mask):
    denom = torch.clamp(mask.sum(), min=1.0)
    return ((pred - target) ** 2 * mask).sum() / denom

def _train_autoencoder(X: np.ndarray, epochs=200, batch_size=64, lr=1e-3, device="cpu"):
    model = Autoencoder(d_in=X.shape[1]).to(device)
    ds = TabularMaskDataset(X, mask_prob=0.0, for_remasker=False)
    dl = DataLoader(ds, batch_size=batch_size, shuffle=True)
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    model.train()
    for _ in range(epochs):
        for xb, obs, _ in dl:
            xb, obs = xb.to(device), obs.to(device)
            pred = model(xb)
            loss = _masked_mse(pred, xb, obs)
            opt.zero_grad(); loss.backward(); opt.step()
    return model

def _train_remasker(X: np.ndarray, epochs=300, batch_size=64, lr=1e-3, mask_prob=0.15, device="cpu",
                    d_model=128, nhead=8, num_layers=4, dim_feedforward=256):
    model = RemaskerTransformer(d_in=X.shape[1], d_model=d_model, nhead=nhead,
                                num_layers=num_layers, dim_feedforward=dim_feedforward).to(device)
    ds = TabularMaskDataset(X, mask_prob=mask_prob, for_remasker=True)
    dl = DataLoader(ds, batch_size=batch_size, shuffle=True)
    opt = torch.optim.AdamW(model.parameters(), lr=lr)
    model.train()
    for _ in range(epochs):
        for xb, obs, addl in dl:
            xb, obs, addl = xb.to(device), obs.to(device), addl.to(device)
            with torch.no_grad():
                target = xb.clone()
            train_mask = torch.clamp((1.0 - obs) + addl, 0, 1)
            pred = model(xb)
            loss = _masked_mse(pred, target, train_mask)
            opt.zero_grad(); loss.backward(); opt.step()
    return model

def _impute_with_model(model: nn.Module, X_std: np.ndarray, means: np.ndarray, stds: np.ndarray, device="cpu") -> np.ndarray:
    model.eval()
    with torch.no_grad():
        xb = torch.tensor(np.nan_to_num(X_std, nan=0.0), dtype=torch.float32, device=device)
        pred_std = model(xb).cpu().numpy()
    obs = ~np.isnan(X_std)
    X_completed_std = X_std.copy()
    X_completed_std[~obs] = pred_std[~obs]
    return _destandardize(X_completed_std, means, stds)

# --------------------------
# Public API
# --------------------------
def impute_with_autoencoder(df: pd.DataFrame, epochs=200, batch_size=64, lr=1e-3, device=None):
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    df = df.copy()
    df_num, num_cols = _select_numeric(df)
    if len(num_cols) == 0:
        return df
    df_std, means, stds = _standardize(df_num)
    X_std = df_std.to_numpy(dtype=np.float32)
    model = _train_autoencoder(X_std, epochs, batch_size, lr, device)
    X_imp = _impute_with_model(model, X_std, means, stds, device)
    df[num_cols] = X_imp
    return df

def impute_with_remasker(df: pd.DataFrame, epochs=300, batch_size=64, lr=1e-3, mask_prob=0.15,
                         d_model=128, nhead=8, num_layers=4, dim_feedforward=256, device=None):
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    df = df.copy()
    df_num, num_cols = _select_numeric(df)
    if len(num_cols) == 0:
        return df
    df_std, means, stds = _standardize(df_num)
    X_std = df_std.to_numpy(dtype=np.float32)
    model = _train_remasker(X_std, epochs, batch_size, lr, mask_prob, device,
                            d_model, nhead, num_layers, dim_feedforward)
    X_imp = _impute_with_model(model, X_std, means, stds, device)
    df[num_cols] = X_imp
    return df

