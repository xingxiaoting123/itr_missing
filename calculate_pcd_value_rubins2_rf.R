
library(randomForest)

calculate_pcd_value_rubins2_rf <- function(train_data_list, test_data, p,
                                           ntree = 500, nodesize = 5, seed = 1) {
  features <- paste0("X", seq_len(p))
  
  test_data$group <- as.factor(test_data$group)
  test_data$true_group <- as.factor(test_data$true_group)
  
  m <- length(train_data_list)
  n <- nrow(test_data)
  pred_drg_mat <- matrix(NA_real_, nrow = n, ncol = m)
  pred_pbo_mat <- matrix(NA_real_, nrow = n, ncol = m)
  
  xtest <- test_data[, features, drop = FALSE]
  
  for (i in seq_len(m)) {
    d <- train_data_list[[i]]
    d$group <- as.factor(d$group)
    
    d_drg <- d[d$group == "drg", c(features, "outcome"), drop = FALSE]
    d_pbo <- d[d$group == "pbo", c(features, "outcome"), drop = FALSE]
    
    set.seed(seed + i - 1)
    rf_drg <- randomForest::randomForest(outcome ~ ., data = d_drg, ntree = ntree, nodesize = nodesize)
    rf_pbo <- randomForest::randomForest(outcome ~ ., data = d_pbo, ntree = ntree, nodesize = nodesize)
    
    pred_drg_mat[, i] <- predict(rf_drg, newdata = xtest)
    pred_pbo_mat[, i] <- predict(rf_pbo, newdata = xtest)
  }
  
  pred_drg <- rowMeans(pred_drg_mat, na.rm = TRUE)
  pred_pbo <- rowMeans(pred_pbo_mat, na.rm = TRUE)
  
  predicted_group <- ifelse(pred_drg > pred_pbo, "drg", "pbo")
  
  pcd <- mean(predicted_group == test_data$true_group)
  
  correct_idx <- predicted_group == as.character(test_data$group)
  value <- if (any(correct_idx)) mean(test_data$outcome[correct_idx], na.rm = TRUE) else NA_real_
  
  list(PCD = pcd, Value = value)
}
