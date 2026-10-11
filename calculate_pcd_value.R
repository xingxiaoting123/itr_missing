calculate_pcd_value <- function(train_data, test_data,p) {
  # Ensure that the group variable is a factor
  train_data$group <- as.factor(train_data$group)
  test_data$group <- as.factor(test_data$group)
  
  # Separate the training data into placebo (pbo) and drug (drg) groups
  train_data_pbo <- train_data[train_data$group == 'pbo', ]
  train_data_drg <- train_data[train_data$group == 'drg', ]
  
  # Define the formula for the linear model (outcome ~ X1 + X2 + X3 + X4 + X5)
  lm_formula <- as.formula(paste("outcome ~", paste(paste0("X", 1:p), collapse = " + ")))
  #lm_formula <- "outcome ~ X1+X2" 
  
  ## Modify the simulation data to add the potential outcome for the group people are not assigned to. 
  # Fit linear models for both placebo and drug groups using the training data
  lm_pbo <- lm(lm_formula, data = train_data_pbo)
  lm_drg <- lm(lm_formula, data = train_data_drg)
  
  #browser()
  # Predict outcomes for the test data # actual outcome
  pred_pbo <- predict(lm_pbo, newdata = test_data) ## Not doing the imputation for test data.
  pred_drg <- predict(lm_drg, newdata = test_data)
  
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
  #value <- mean(c(test_data[predicted_group=='drg',]$y_drg,test_data[predicted_group=='pbo',]$y_pbo))#100 rows 
  
  #mean(test_data[predicted_group=='pbo',]$y_pbo)
  #mean(test_data[predicted_group=='drg',]$y_drg)
  ##value <- mean(predicted_value, na.rm=TRUE)
  #plot(density(test_data$y_drg))lines(density(test_data$y_pbo),col="red")
  
  # Return a list with PCD and Value
  return(list(PCD = pcd, Value = value))
}
