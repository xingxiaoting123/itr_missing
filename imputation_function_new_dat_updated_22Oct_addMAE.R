library(dplyr)
library(mice)

# 1. Function for Complete Case Analysis
# Corrected function for Complete Case Imputation
complete_case_imputation <- function(data) {
  # Identify outcome columns to skip during NA removal
  outcome_vars <- c("y_pbo", "y_drg")
  
  # Get all other column names
  other_vars <- setdiff(names(data), outcome_vars)
  
  # Remove rows with missing values in variables other than y_pbo and y_drg
  complete_rows <- complete.cases(data[, other_vars])
  imputed_data <- data[complete_rows, ]
  
  return(imputed_data)
}

# 2. Function for Single Imputation (Fixed Value)
#single_imputation <- function(data, fixed_value = 0) {
  # Replace missing values with a fixed value (default is 0)
#  imputed_data <- data %>% mutate_all(~ifelse(is.na(.), fixed_value, .))
#  return(imputed_data)
#}


# 2. Function for Single Imputation Using Linear Regression via mice package
single_imputation <- function(data) {
  
  library(mice)
  
  p <- length(grep("^X", colnames(data), value = TRUE))
  
  # Specify variables to exclude from imputation
  exclude_vars <- c("subject", "group", "true_group", "outcome", paste0("X", p))
  
  # Variables to be imputed
  impute_vars <- setdiff(names(data), exclude_vars)
  
  # Create a method vector for mice
  method <- make.method(data)
  
  # Set the imputation method for each variable
  for (var in names(data)) {
    if (var %in% exclude_vars) {
      method[var] <- ""  # Skip imputation for excluded variables
    } else if (is.numeric(data[[var]])) {
      method[var] <-  "norm"   # Use linear regression 
    } else {
      method[var] <- ""  # Skip non-numeric variables not in exclude_vars
    }
  }
  
  # Create a predictor matrix
  pred <- make.predictorMatrix(data)
  
  # Exclude variables from being used as predictors
  pred[, exclude_vars] <- 0  # Exclude these variables as predictors
  pred[exclude_vars, ] <- 0  # Exclude these variables from being imputed
  
  # Perform single imputation using mice with m = 1 and maxit = 1
  imputed_data <- mice(data, method = method, predictorMatrix = pred, m = 1, maxit = 1, printFlag = FALSE)
  #browser()
  # Extract the completed dataset
  completed_data <- complete(imputed_data, action = 1)
  
  return(completed_data)
}


# 3. Function for Mean Imputation
mean_imputation <- function(data) {
  # Replace missing values with the mean of each column
  imputed_data <- data %>%
    mutate(across(everything(), ~ifelse(is.na(.), mean(., na.rm = TRUE), .)))
  return(imputed_data)
}

# 4# multiple imputation
mice_imputation <- function(data, m = 10, maxit = 50, seed = 123) {
  library(mice)
  
  # Identify predictor variables (X1 to X(p-1))
  p <- length(grep("^X", colnames(data), value = TRUE))
  vars_to_impute <- c(paste0("X", 1:(p)), "outcome")
  
  # Build predictor matrix
  pred_matrix <- make.predictorMatrix(data)
  pred_matrix[, c("subject", "group", "true_group", paste0("X", p))] <- 0
  pred_matrix[c("subject", "group", "true_group", paste0("X", p)), ] <- 0
  
  # Set imputation method
  meth <- make.method(data)
  meth[c("subject", "group", "true_group", paste0("X", p))] <- ""
  # Adaptive method choice
  if (p <= 20) {
    meth[vars_to_impute] <- "norm.nob"
    ridge_penalty <- 0
  } else {
    meth[vars_to_impute] <- "ridge"
    ridge_penalty <- 1e-4
  }
  
  # Run MICE
  mice_result <- mice(
    data = data,
    m = m,
    maxit = maxit,
    method = meth,
    seed = seed,
    ridge = ridge_penalty,
    predictorMatrix = pred_matrix,
    printFlag = FALSE)  # July 4th change to meth[vars_to_impute] <- "norm.nob", previous is "norm" there will be error   x has missing values; consider using makeX() to impute them
  
  # Extract imputed datasets
  imputed_data_independent <- lapply(1:m, function(i) complete(mice_result, i))
  
  # Final safety check: forcefully fill any residual NAs using mean imputation
  
  # imputed_data_independent <- lapply(imputed_data_independent, function(df) {
  #   df[] <- lapply(df, function(col) {
  #     if (is.numeric(col) && anyNA(col)) {
  #       col[is.na(col)] <- mean(col, na.rm = TRUE)
  #     }
  #     return(col)
  #   })
  #   return(df)
  # })
  
  # Return both long and list formats
  return(list(
    imputed_data = complete(mice_result, action = "long", include = TRUE),
    imputed_data_independent = imputed_data_independent
  ))
}


mice.impute.ridge <- function(y, ry, x, ridge = 1e-4, ...) {
  x <- cbind(1, as.matrix(x))  # Add intercept
  xtx <- t(x[ry, ]) %*% x[ry, ]
  penalty <- ridge * diag(ncol(xtx))
  coef <- solve(xtx + penalty, t(x[ry, ]) %*% y[ry])
  yhat <- x[!ry, ] %*% coef
  yhat + rnorm(sum(!ry), 0, sqrt(var(y[ry] - x[ry, ] %*% coef)))
}


#5.mice imputation lasso norm

mice_imputation_lasso_norm <- function(data, m = 10, maxit = 50, seed = 123) {
  library(mice)
  
  # Get number of X columns (exclude the last one Xp used as latent)
  p <- length(grep("^X", colnames(data), value = TRUE))
  vars_to_impute <- c(paste0("X", 1:(p - 1)), "outcome")  # Exclude Xp
  
  # Build predictor matrix
  pred_matrix <- make.predictorMatrix(data)
  pred_matrix[, c("subject", "group", "true_group", paste0("X", p))] <- 0
  pred_matrix[c("subject", "group", "true_group", paste0("X", p)), ] <- 0
  
  # Set imputation methods
  meth <- make.method(data)
  meth[c("subject", "group", "true_group", paste0("X", p))] <- ""   # No imputation
  meth[paste0("X", 1:(p - 1))] <- "lasso.norm"                       # Lasso for predictors
  meth["outcome"] <- "norm"                                         # Standard for outcome
  
  # Run MICE with lasso imputation
  mice_result <- mice(
    data,
    m = m,
    maxit = maxit,
    method = meth,
    seed = seed,
    predictorMatrix = pred_matrix,
    printFlag = TRUE,
    lasso.norm.nfolds = 5,
    lasso.norm.makeX = TRUE# optional: reduce folds for faster CV
  )
  
  # Extract imputed datasets
  imputed_data_independent <- lapply(1:m, function(i) complete(mice_result, i))
  
  return(list(
    imputed_data = complete(mice_result, action = "long", include = TRUE),
    imputed_data_independent = imputed_data_independent
  ))
}


#####Below method no use########################################################################


# 5. Main Driver Function to Call the Appropriate Imputation Method
impute_missing_data <- function(data, method) {
  
  # Exclude specific columns that are not suitable for imputation (if necessary)
  columns_to_exclude <- c("subject", "group")
  imputation_data <- data
  
  # Select the imputation method
  if (method == "complete_case") {
    imputed_data <- complete_case_imputation(imputation_data)
    
  } else if (method == "single_imputation") {
    imputed_data <- single_imputation(imputation_data)
    
  } else if (method == "mean_imputation") {
    imputed_data <- mean_imputation(imputation_data)
    
  } else if (method == "mice_imputation") {
    imputed_data <- mice_imputation(imputation_data)
    
  } else {
    stop("Invalid method specified. Choose from 'complete_case', 'single_imputation', 'mean_imputation', or 'mice_imputation'.")
  }
  
  return(imputed_data)
}

# 6. MissForest Imputation
missforest_imputation <- function(data, maxiter = 10, ntree = 100, seed = 123) {
  library(missForest)
  
  # Preserve outcome and ID columns
  outcome_vars <- c("y_pbo", "y_drg", "outcome")
  exclude_vars <- c("subject", "group", "true_group")
  
  # Select numeric predictors for imputation
  impute_vars <- setdiff(names(data), exclude_vars)
  
  # Subset data for imputation
  data_impute <- data[, impute_vars, drop = FALSE]
  
  set.seed(seed)
  imputed_result <- missForest(
    data_impute,
    maxiter = maxiter,
    ntree = ntree,
    verbose = TRUE
  )
  
  # Replace the imputed columns in the original dataset
  data[, impute_vars] <- imputed_result$ximp
  
  return(data)
}

# 7. Autoencoder / ReMasker Imputation via PyTorch (Python) -- OPTIONAL.
# Loaded only if reticulate + impute_nn.py + PyTorch are available, so the
# rest of the workflow runs on machines without a Python/PyTorch setup.
.remasker_available <- FALSE
if (requireNamespace("reticulate", quietly = TRUE) && file.exists("impute_nn.py")) {
  .remasker_available <- tryCatch({
    reticulate::source_python("impute_nn.py"); TRUE
  }, error = function(e) {
    message("ReMasker/PyTorch not available; ReMasker imputation will be skipped. (",
            conditionMessage(e), ")"); FALSE
  })
}
autoencoder_pytorch_imputation <- function(data, epochs = 50L) {
  # Convert to data frame to ensure compatibility with Python
  data <- as.data.frame(data)
  
  # Call the Python function
  imputed_df <- impute_with_autoencoder(data, epochs = as.integer(epochs))
  
  # Convert back to R data frame
  imputed_df <- as.data.frame(imputed_df)
  
  return(imputed_df)
}

# 8. ReMasker (Transformer) Imputation via PyTorch
remasker_pytorch_imputation <- function(data, epochs = 100L, mask_prob = 0.15) {
  data <- as.data.frame(data)
  
  imputed_df <- impute_with_remasker(
    data,
    epochs = as.integer(epochs),
    mask_prob = as.numeric(mask_prob)
  )
  
  imputed_df <- as.data.frame(imputed_df)
  return(imputed_df)
}

