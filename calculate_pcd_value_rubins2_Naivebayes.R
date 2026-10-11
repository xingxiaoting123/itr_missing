
# Apply Rubin's rules for Naive Bayes estimation
calculate_pcd_value_rubins2_naivebayes <- function(train_data_list, test_data, p) {
  features <- paste0("X", 1:p)
  
  get_nb_coef <- function(data) {
    mu_x <- colMeans(data[, features])
    mu_y <- mean(data$outcome)
    diag_var <- apply(data[, features], 2, var)
    D_x_inv <- diag(1 / diag_var)
    cov_xy <- cov(data[, features], data$outcome)
    beta <- D_x_inv %*% cov_xy
    beta0 <- mu_y - t(beta) %*% mu_x
    return(list(beta0 = as.numeric(beta0), beta = beta))
  }
  
  # Store betas and variances for Rubin's combination
  coef_pbo_list <- list()
  var_pbo_list <- list()
  coef_drg_list <- list()
  var_drg_list <- list()
  
  for (i in seq_along(train_data_list)) {
    current_train <- train_data_list[[i]]
    
    coef_pbo <- get_nb_coef(current_train[current_train$group == "pbo", ])
    coef_drg <- get_nb_coef(current_train[current_train$group == "drg", ])
    
    coef_pbo_vec <- c(coef_pbo$beta0, coef_pbo$beta)
    coef_drg_vec <- c(coef_drg$beta0, coef_drg$beta)
    
    coef_pbo_list[[i]] <- coef_pbo_vec
    coef_drg_list[[i]] <- coef_drg_vec
    
    var_pbo_list[[i]] <- rep(NA, length(coef_pbo_vec))
    var_drg_list[[i]] <- rep(NA, length(coef_drg_vec))
  }
  
  pooled_pbo <- rubins_rule_combine2(coef_pbo_list, var_pbo_list)
  pooled_drg <- rubins_rule_combine2(coef_drg_list, var_drg_list)
  
  test_data$group <- as.factor(test_data$group)
  test_data$true_group <- as.factor(test_data$true_group)
  
  x_mat <- as.matrix(test_data[, features])
  pred_pbo <- pooled_pbo$estimate[1] + x_mat %*% pooled_pbo$estimate[-1]
  pred_drg <- pooled_drg$estimate[1] + x_mat %*% pooled_drg$estimate[-1]
  
  predicted_group <- ifelse(pred_drg > pred_pbo, 'drg', 'pbo')
  pcd <- sum(predicted_group == test_data$true_group) / nrow(test_data)
  
  correct_classifications <- test_data[predicted_group == test_data$group, ]
  value <- if (nrow(correct_classifications) > 0) mean(correct_classifications$outcome, na.rm = TRUE) else NA
  
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
