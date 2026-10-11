
#library(e1071)

calculate_pcd_value_naivebayes<- function(train_data,test_data, p) {

  features <- paste0("X", 1:p)
  
  get_nb_coef <- function(data) {
    # Estimate means and variances under NB assumptions
    mu_x <- colMeans(data[, features])
    mu_y <- mean(data$outcome)
    # Naive Bayes: diagonal covariance matrix：diagnol of psix 
    diag_var <- apply(data[, features], 2, var)
    D_x_inv <- diag(1 / diag_var)
    
    cov_xy <- cov(data[, features], data$outcome)
    beta <- D_x_inv %*% cov_xy #vector p*1
    beta0 <- mu_y - t(beta) %*% mu_x #scalar
    return(list(beta0 = as.numeric(beta0), beta = beta))
  }
  
  # Fit models for each group
  drg_coef <- get_nb_coef(subset(train_data, group == "drg"))
  pbo_coef <- get_nb_coef(subset(train_data, group == "pbo"))
  
  # Predict outcomes for each test subject under both treatments
  x_mat <- as.matrix(test_data[, features])
  pred_drg <- as.numeric(drg_coef$beta0 + x_mat %*% drg_coef$beta)
  pred_pbo <- as.numeric(pbo_coef$beta0 + x_mat %*% pbo_coef$beta)
  
  
  # Assign each test subject to the group with the better predicted outcome
  predicted_group <- ifelse(pred_drg > pred_pbo, 'drg', 'pbo')
  ##predicted_value <- ifelse(pred_drg > pred_pbo, pred_drg, pred_pbo)
  
  ####################################################################################
  
  # Calculate PCD (Percent Correct Decisions) If test have the missing, return to NA 
  pcd <- sum(predicted_group == test_data$true_group) / nrow(test_data)
  
  # Calculate Value (mean outcome for correctly classified subjects) using the 
  # figure out this value calculation value , higher can be better. 
  correct_classifications <- test_data[predicted_group == test_data$group, ]# only rows who make the correct decision
  ## modify on 15Apr2025 test_data$group
  value <- mean(correct_classifications$outcome, na.rm = TRUE) 
  
  # Return a list with PCD and Value
  return(list(PCD = pcd, Value = value))
}
