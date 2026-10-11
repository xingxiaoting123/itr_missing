library(randomForest)

calculate_pcd_value_rf <- function(train_data, test_data, p) {
  features <- paste0("X", 1:p)
  
  # Ensure group and true_group are factors
  train_data$group <- as.factor(train_data$group)
  test_data$group <- as.factor(test_data$group)
  test_data$true_group <- as.factor(test_data$true_group)
  
  # Fit outcome models for each treatment group
  model_drg <- randomForest(outcome ~ ., data = train_data[train_data$group == "drg", c(features, "outcome")])
  model_pbo <- randomForest(outcome ~ ., data = train_data[train_data$group == "pbo", c(features, "outcome")])
  
  # Predict outcomes under each treatment
  xtest <- test_data[, features]
  pred_drg <- predict(model_drg, newdata = xtest)
  pred_pbo <- predict(model_pbo, newdata = xtest)
  
  # Assign treatment with higher predicted outcome
  predicted_group <- ifelse(pred_drg > pred_pbo, "drg", "pbo")
  
  # Calculate PCD: agreement with optimal group
  pcd <- mean(predicted_group == test_data$true_group)
  
  # Calculate Value: observed outcome where predicted == actual assignment
  correct_classifications <- test_data[predicted_group == test_data$group, ]
  value <- mean(correct_classifications$outcome, na.rm = TRUE)
  
  return(list(PCD = pcd, Value = value))
}

