# =====================================================================
# reproduce_example.R
# Runs ONE representative replicate of the paper's workflow, using the
# paper's own functions:
#   simulate (AR(1)) -> 50% MAR missingness -> impute -> fit each ITR
#   learner -> assign treatment -> compute PCD and Value.
#
# It reproduces the core methodology, not the full 500-replicate study.
# Run (from this folder):
#   source("reproduce_example.R")
# =====================================================================

## ---- packages ----
pkgs <- c("MASS", "dplyr", "tidyr", "mice", "missForest",
          "glmnet", "randomForest", "e1071")
invisible(lapply(pkgs, function(p) suppressMessages(
  if (!require(p, character.only = TRUE)) stop("Please install package: ", p))))

## ---- source the paper's own functions ----
source("compute_c_star_value_star.R")     # population c* / value* helper
source("dat_gen_new_add_star_28Dec.R")    # generate_simulation_data()
source("imputation_function_new_dat_updated_22Oct_addMAE.R")  # imputers
source("missingness.R")                   # generate_missing_data()
for (f in c("calculate_pcd_value.R",            "calculate_pcd_value_rubins2.R",
            "calculate_pcd_value_lasso.R",      "calculate_pcd_value_rubins2_lasso.R",
            "calculate_pcd_value_ridge.R",      "calculate_pcd_value_rubins2_ridge.R",
            "calculate_pcd_value_naivebayes.R", "calculate_pcd_value_rubins2_Naivebayes.R",
            "calculate_pcd_value_rf.R",         "calculate_pcd_value_rubins2_rf.R")) source(f)

## ---- representative setting (AR(1), 50% MAR) ----
p <- 10; n_train <- 50; n_test <- 500; rho <- 0.3
cov_type <- "ar1"; overall_missing <- 0.5
set.seed(1)
beta1 <- matrix(rnorm(p + 1, mean = -0.3, sd = 5), ncol = 1)
beta2 <- matrix(rnorm(p + 1, mean = -0.3, sd = 5), ncol = 1)

## ---- 1. generate one complete dataset (p+1: last column is the latent) ----
full <- generate_simulation_data(n = n_train + n_test, p = p + 1, rho = rho,
                                 cov_type = cov_type, beta1 = beta1, beta2 = beta2)$data
idx        <- sample(seq_len(nrow(full)), 2 * n_train)
train_data <- full[idx, ]
test_data  <- full[-idx, ]
test_data[[paste0("X", p + 1)]] <- NULL          # drop latent from the test set
test_data$group      <- as.factor(test_data$group)
test_data$true_group <- as.factor(test_data$true_group)

## ---- 2. introduce 50% MAR missingness in the training covariates ----
## the original generate_missing_data() reads the globals `dat` and `p`
dat <- train_data
MAR_data <- generate_missing_data(dat, "MAR", overall_missing)

## helper: drop the latent column and make group a factor
prep <- function(d) { d <- as.data.frame(d); nm <- paste0("X", p + 1)
  if (nm %in% names(d)) d[[nm]] <- NULL; d$group <- as.factor(d$group); d }

## ---- learners and imputers (the paper's functions) ----
model_functions <- list(lm = calculate_pcd_value, lasso = calculate_pcd_value_lasso,
  ridge = calculate_pcd_value_ridge, NB = calculate_pcd_value_naivebayes, rf = calculate_pcd_value_rf)
rubins2_functions <- list(lm = calculate_pcd_value_rubins2, lasso = calculate_pcd_value_rubins2_lasso,
  ridge = calculate_pcd_value_rubins2_ridge, NB = calculate_pcd_value_rubins2_naivebayes,
  rf = calculate_pcd_value_rubins2_rf)
imputers <- list("Complete Case" = complete_case_imputation, "Mean" = mean_imputation,
                 "Single" = single_imputation, "MissForest" = missforest_imputation)
if (isTRUE(.remasker_available)) imputers[["ReMasker"]] <- remasker_pytorch_imputation

## ---- 3-6. impute, fit each learner, assign, compute PCD/Value ----
res <- data.frame()
run <- function(imp_label, train, fns) {
  for (m in names(fns)) {
    r <- tryCatch(fns[[m]](train, test_data, p = p),
                  error = function(e) list(PCD = NA, Value = NA))
    res <<- rbind(res, data.frame(Imputation = imp_label, Learner = m,
                                  PCD = r$PCD, Value = r$Value))
  }
}
run("No missing", prep(train_data), model_functions)                 # benchmark
for (nm in names(imputers)) run(nm, prep(imputers[[nm]](MAR_data)), model_functions)
mice_train <- mice_imputation(MAR_data)$imputed_data_independent     # multiple imputation
run("Multiple", mice_train, rubins2_functions)

## ---- report ----
cat(sprintf("\nITR value, one replicate (cov=%s, p=%d, rho=%.1f, MAR %.0f%% missing)\n",
            cov_type, p, rho, 100 * overall_missing))
res$Value <- round(res$Value, 3); res$PCD <- round(res$PCD, 3)
print(res, row.names = FALSE)
