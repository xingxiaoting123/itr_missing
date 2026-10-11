# itr_missing

Code for the paper **"Linear biosignatures for individualized treatment rules with missing baseline covariates"** (Xiaoting Xing and Thaddeus Tarpey).


This package contains the **actual functions used in the paper** together
with a short driver that runs one representative example end to end:

```
simulate (AR(1)) -> 50% MAR missingness -> impute -> fit each ITR learner
        -> assign treatment -> compute PCD and Value
```

It lets a reader run the core methodology on a simulated example. It does
**not** rerun the full 500-replicate Monte-Carlo study or the real-data
(EMBARC) analysis: the EMBARC data are obtained from the NIMH Data Archive
under a data-use agreement and are not redistributed here.

## Software

- R version 4.x (developed under R 4.3+)

Required R packages:
`MASS`, `dplyr`, `tidyr`, `mice`, `missForest`, `glmnet`, `randomForest`, `e1071`.

```r
install.packages(c("MASS","dplyr","tidyr","mice","missForest",
                   "glmnet","randomForest","e1071"))
```

ReMasker is optional. It additionally needs Python + PyTorch and is called
from R through `reticulate` using `impute_nn.py`. If PyTorch is not
available the example detects this automatically and simply skips ReMasker;
everything else still runs.

## Files

Driver / example:

- `reproduce_example.R` — runs one replicate (AR(1), 50% MAR) and prints
  PCD and Value for every ITR learner under no-missing, each single
  imputation, and multiple imputation.
- `missingness.R` — MCAR/MAR/MNAR missingness generator (from the driver).

Paper's own functions (used unchanged):

- `dat_gen_new_add_star_28Dec.R` — `generate_simulation_data()`: baseline
  covariates (AR(1) / equal / real covariance), potential outcomes,
  treatment assignment, and the optimal rule.
- `compute_c_star_value_star.R` — population-level c\* / value\* helper.
- `imputation_function_new_dat_updated_22Oct_addMAE.R` — imputation:
  complete case, mean, single, multiple (`mice`), MissForest, and ReMasker
  (`impute_nn.py`).
- `calculate_pcd_value*.R` — the ITR learners and value function for
  linear model, lasso, ridge, Naive Bayes, and random forest, each with a
  single-imputation version and a `_rubins2` version that pools multiple
  imputations by Rubin's rule.
- `impute_nn.py` — ReMasker (PyTorch), called via `reticulate` (optional).

## To run

From within this folder:

```r
source("reproduce_example.R")
```

Output is a table of `PCD` (proportion of correct decisions vs the optimal
rule) and `Value` (estimated mean outcome under the estimated rule) for
each learner and imputation strategy.

## Notes

- The example setting (AR(1), p = 10, rho = 0.3, 50% MAR) is at the top of
  `reproduce_example.R`; change `cov_type`, `p`, `rho`, sample sizes, or
  `overall_missing` there to explore other settings from the paper. To run
  the correctly-specified, misspecified, and low-signal designs across many
  replicates, the same functions are looped over seeds (as in the full
  study); only one replicate is shown here for brevity.

## Note on the real-data covariance setting

`generate_simulation_data(cov_type = "real")` reads `cortical_thickness_data.csv`,
derived from the EMBARC study. Those data are not redistributed here; they are
available from the NIMH Data Archive (Collection #2199) under a data use
agreement. The `"ar1"` and `"equal"` settings need no external data.
