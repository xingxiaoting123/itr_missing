library(glmnet)

calculate_pcd_value_lasso <- function(train_data, test_data, p=p) {
  # Ensure that the group variable is a factor
  train_data$group <- as.factor(train_data$group)
  test_data$group <- as.factor(test_data$group)
  
  # Separate the training data into placebo (pbo) and drug (drg) groups
  train_data_pbo <- train_data[train_data$group == 'pbo', ]
  train_data_drg <- train_data[train_data$group == 'drg', ]
  
  #browser()
  #X_train_pbo <- as.matrix(train_data_pbo[, paste0("X", 1:(p-1))])  # Predictor matrix for placebo
  #X_train_drg <- as.matrix(train_data_drg[, paste0("X", 1:(p-1))])  # Predictor matrix for drug
  X_train_pbo <- as.matrix(train_data_pbo[, paste0("X", 1:(p))])  # This is for non-missing generation
  X_train_drg <- as.matrix(train_data_drg[, paste0("X", 1:(p))])  # This is for non-missing generation
  y_train_pbo <- train_data_pbo$outcome
  y_train_drg <- train_data_drg$outcome
  
  # Fit LASSO regression using cross-validation
  #cv_pbo <- cv.glmnet(X_train_pbo, y_train_pbo, alpha = 1,nfolds = 3)  # LASSO (alpha = 1)
  #cv_drg <- cv.glmnet(X_train_drg, y_train_drg, alpha = 1,nfolds = 3)
  
  # Handle small sample size: set nfolds = min(n, 5)
  nfolds_pbo <- min(nrow(X_train_pbo), 5)
  nfolds_drg <- min(nrow(X_train_drg), 5)
  
  # Fit LASSO regression using cross-validation with adjusted nfolds
  # browser()
  cv_pbo <- cv.glmnet(X_train_pbo, y_train_pbo, alpha = 1, nfolds = nfolds_pbo)  
  cv_drg <- cv.glmnet(X_train_drg, y_train_drg, alpha = 1, nfolds = nfolds_drg)
  
  # Select optimal lambda
  best_lambda_pbo <- cv_pbo$lambda.min
  best_lambda_drg <- cv_drg$lambda.min
  
  # Predict outcomes for the test data
  # X_test <- as.matrix(test_data[, paste0("X", 1:(p-1))])  # Test predictors
  X_test <- as.matrix(test_data[, paste0("X", 1:(p))])  # Test predictors
  pred_pbo <- predict(cv_pbo$glmnet.fit, s = best_lambda_pbo, newx = X_test)
  pred_drg <- predict(cv_drg$glmnet.fit, s = best_lambda_drg, newx = X_test)
  
  # Assign each test subject to the group with the better predicted outcome
  predicted_group <- ifelse(pred_drg > pred_pbo, 'drg', 'pbo')
  
  # Calculate PCD (Percent Correct Decisions)
  pcd <- sum(predicted_group == test_data$true_group) / nrow(test_data)
  
  # Calculate Value (mean outcome for correctly classified subjects)
  correct_classifications <- test_data[predicted_group == test_data$group, ]
  value <- mean(correct_classifications$outcome, na.rm = TRUE)
  
  # Return a list with PCD and Value
  return(list(PCD = pcd, Value = value))
}

