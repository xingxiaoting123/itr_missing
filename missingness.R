# =====================================================================
# missingness.R
# MCAR / MAR / MNAR missingness generator -- the ORIGINAL function used in
# the simulation study, kept verbatim. As in the study it reads the global
# objects `dat` (the current training dataset) and `p` (number of
# covariates), so set those before calling, e.g.
#     dat <- train_data ; p <- 10
#     MAR_data <- generate_missing_data(dat, "MAR", 0.5)
# =====================================================================
generate_missing_data <- function(data, missing_type, missing_rate) {

  # Create a copy of the data
  data_missing <- dat

  # Exclude outcome columns and only include columns with names that contain "X"
  eligible_cols <- grep("^X", colnames(data_missing), value = TRUE)[1:p]  # if gnerate p dimention, here use p-1
  # 30% of the rows in the dataset have at least one missing value in these specified covariates.
  # Determine number of predictors (columns to introduce missingness)
  p1 <- length(eligible_cols) #last column is excluded for x_latent variables for generate
  n <- nrow(dat)

  # Set up varying missingness rates for each predictor
  # v_missing
  # 0.0000000 0.0000000 0.3333333 0.6666667 1.0000000 (here is 100% missing)

  #v_missing <- c(rep(0, floor(p / 2)), (1:(p - floor(p / 2))) / (p - floor(p / 2)))
  v_missing <- c(rep(0, floor(p1 / 2)), (1:(p1 - floor(p1 / 2))) / ((p1 - floor(p1 / 2)) + 1))

  # Introduce missing data based on missing type
  if (missing_type == "MCAR") {
    # Generate missing IDs for rows based on binomial distribution
    missing_rows <- (rbinom(n, size = 1, prob = missing_rate) == 1)
    random_values <- matrix(runif(sum(missing_rows) * p1), nrow = sum(missing_rows), ncol = p1)
    missing_matrix <- (random_values < matrix(v_missing, nrow = sum(missing_rows), ncol = p1, byrow = TRUE))
    row_indices <- which(missing_rows)
    data_missing[row_indices, eligible_cols][missing_matrix] <- NA
  }

  else if (missing_type == "MAR"){
    # Define the missing rate based on X2 covariate values (e.g., only apply to large than X2 median value) only if both conditions are TRUE
    # missing_rows <- (rbinom(n, size = 1, prob = missing_rate) == 1) &
    #  (data_missing[["X2"]] > quantile(data_missing[["X2"]], 0.2))
    #Sep13 modify the MAR
    # Calibrate MAR probability depending on X2 and scale to target missing_rate
    raw <- plogis(-data_missing$X2)   # symmetric around 0 for X2 ~ N(0, sigma)
    prob_missing <- pmin(raw * (missing_rate / mean(raw)), 0.999)

    # Row is eligible if Bernoulli draw with probability depending on X2
    missing_rows <- (rbinom(n, size = 1, prob = prob_missing) == 1)

    random_values <- matrix(runif(sum(missing_rows) * p1), nrow = sum(missing_rows), ncol = p1)
    missing_tmp <- matrix(v_missing, nrow = sum(missing_rows), ncol = p1, byrow = TRUE) #column missing matrix
    missing_matrix <- (random_values < missing_tmp) # column missing true and false
    row_indices <- which(missing_rows)
    data_missing[row_indices, eligible_cols][missing_matrix] <- NA
  }
  else if (missing_type == "MNAR") {

    # Generate row indices for missing data based on a binomial distribution
    # Identify the latent variable as X(p+1)
    latent_var <- paste0("X", p + 1)

    # Use it in the missing row definition
    missing_rows <- (rbinom(n, size = 1, prob = missing_rate) == 1) &
      (data_missing[[latent_var]] > quantile(data_missing[[latent_var]], 0.2))
    random_values <- matrix(runif(sum(missing_rows) * p1), nrow = sum(missing_rows), ncol = p1)
    missing_tmp <- matrix(v_missing, nrow = sum(missing_rows), ncol = p1, byrow = TRUE)
    missing_matrix <- (random_values < missing_tmp)
    row_indices <- which(missing_rows)
    data_missing[row_indices, eligible_cols][missing_matrix] <- NA
  }

  return(data_missing)
}
