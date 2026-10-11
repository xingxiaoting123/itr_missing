calculate_pcd_value_rubins2_lasso <- function(train_data_list, test_data,p=p) {
  # Initialize lists to store coefficients and variances for each imputed dataset
  coef_pbo_list <- list()
  var_pbo_list <- list()
  coef_drg_list <- list()
  var_drg_list <- list()
  
  # Loop through each imputed dataset and fit linear models for placebo and drug groups
  for (i in seq_along(train_data_list)) {
    current_train_data <- train_data_list[[i]]
    current_train_data$group <- as.factor(current_train_data$group)
    
    # Separate the training data into placebo (pbo) and drug (drg) groups
    train_data_pbo <- current_train_data[current_train_data$group == 'pbo', ]
    train_data_drg <- current_train_data[current_train_data$group == 'drg', ]
    
    # Define the formula for the linear model
    # lm_formula <- outcome ~ X1 + X2 + X3 + X4 + X5
    # p <- length(grep("X", names(current_train_data)))
    lm_formula <- as.formula(paste("outcome ~", paste(paste0("X", 1:(p)), collapse = " + ")))
    
    # Define predictor matrix and outcome
    X_train_pbo <- as.matrix(train_data_pbo[, paste0("X", 1:(p))])
    y_train_pbo <- train_data_pbo$outcome
    X_train_drg <- as.matrix(train_data_drg[, paste0("X", 1:(p))])
    y_train_drg <- train_data_drg$outcome
    
    # Cross-validation folds (handle small sample sizes)
    nfolds_pbo <- min(nrow(X_train_pbo), 5)
    nfolds_drg <- min(nrow(X_train_drg), 5)
    
    # Fit LASSO using cross-validation
    #browser()
    cv_pbo <- cv.glmnet(X_train_pbo, y_train_pbo, alpha = 1, nfolds = nfolds_pbo)
    cv_drg <- cv.glmnet(X_train_drg, y_train_drg, alpha = 1, nfolds = nfolds_drg)
    
    # Get best lambda and coefficients
    coef_pbo <- as.vector(coef(cv_pbo, s = "lambda.min"))
    coef_drg <- as.vector(coef(cv_drg, s = "lambda.min"))
    
    # Assign names to match model.matrix
    coef_names <- c("(Intercept)", paste0("X", 1:(p)))
    names(coef_pbo) <- coef_names
    names(coef_drg) <- coef_names
    
    # Store coefficients; variances can be set as NA or approximated
    coef_pbo_list[[i]] <- coef_pbo
    coef_drg_list[[i]] <- coef_drg
    
    # Variance estimation for LASSO is complex, so set NA for now
    var_pbo_list[[i]] <- rep(NA, length(coef_pbo))
    var_drg_list[[i]] <- rep(NA, length(coef_drg))
    
  }
  
  # Apply Rubin's rule to combine the coefficients from multiple imputed datasets
  pooled_pbo <- rubins_rule_combine2(coef_pbo_list, var_pbo_list)
  pooled_drg <- rubins_rule_combine2(coef_drg_list, var_drg_list)
  
  # Extract the pooled coefficients
  pooled_coef_pbo <- pooled_pbo$estimate
  pooled_coef_drg <- pooled_drg$estimate
  
  # Initialize vectors to store PCD and Value for each test dataset
  #pcd_list <- numeric(length(test_data_list))
  #value_list <- numeric(length(test_data_list))
  
  #for (i in seq_along(test_data_list)) {
    # Use the pooled coefficients to generate predictions for test data
    # test_data <- test_data_list
    test_data$group <- as.factor(test_data$group)
    test_data$true_group <- as.factor(test_data$true_group)
    
    # Prepare the design matrix for prediction
    #model_matrix <- model.matrix(~ X1 + X2 + X3 + X4 + X5, data = test_data)
    predictors <- paste0("X", 1:(p))  # X1, X2, X3, ..., X5
    # Create the formula dynamically using as.formula
    formula <- as.formula(paste("~", paste(predictors, collapse = " + ")))
    # Apply model.matrix with the correct formula and dataset
    model_matrix <- model.matrix(formula, data = test_data)
    
    # Ensure the coefficients are aligned with the model matrix
    coef_names <- colnames(model_matrix)
    pooled_coef_pbo_ordered <- pooled_coef_pbo[coef_names]
    pooled_coef_drg_ordered <- pooled_coef_drg[coef_names]
    
    # Compute predicted outcomes for both placebo and drug groups
    pred_pbo <- as.vector(model_matrix %*% pooled_coef_pbo_ordered)
    pred_drg <- as.vector(model_matrix %*% pooled_coef_drg_ordered)
    
    # Assign each test subject to the group with the better predicted outcome
    predicted_group <- ifelse(pred_drg > pred_pbo, 'drg', 'pbo')
    
    # Calculate PCD (Percent Correct Decisions)
    pcd <- sum(predicted_group == test_data$true_group) / nrow(test_data)
    
    # Calculate Value (mean outcome for correctly classified subjects)
    correct_classifications <- test_data[predicted_group == test_data$group, ]
    if (nrow(correct_classifications) > 0) {
      value <- mean(correct_classifications$outcome, na.rm = TRUE)
    } else {
      value <- NA
    }
  #}
  
  # Return the final PCD and Value based on Rubin's pooled coefficients
  return(list(PCD = pcd, Value = value))
}

rubins_rule_combine2 <- function(coef_list, var_list) {
  m <- length(coef_list)
  
  # Convert lists to matrices
  coef_matrix <- do.call(cbind, coef_list)  # Each column is an imputation
  var_matrix <- do.call(cbind, var_list)    # Each column is an imputation
  
  # Calculate the pooled estimates (mean of estimates)
  Q_bar <- rowMeans(coef_matrix)
  
  # Calculate the within-imputation variances (mean of variances)
  U_bar <- rowMeans(var_matrix)
  
  # Calculate the between-imputation variances
  coef_diff <- coef_matrix - Q_bar
  B <- rowSums(coef_diff^2) / (m - 1)
  
  # Total variance
  T <- U_bar + (1 + 1/m) * B
  
  # Standard errors
  pooled_se <- sqrt(T)
  
  # Return a list containing the pooled estimates and standard errors
  return(list(estimate = Q_bar, std_error = pooled_se))
}
