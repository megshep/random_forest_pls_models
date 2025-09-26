#RANDOM FOREST
  
# Install randomForest 
  install.packages("randomForest")

# Load necessary libraries
library(randomForest)

# Impute missing values with column means (for numerical variables)
exp_low_imputed <- exp_low

# Loop through each column and replace NAs with the mean of that column
for (col in colnames(exp_low_imputed)) {
  exp_low_imputed[[col]][is.na(exp_low_imputed[[col]])] <- mean(exp_low_imputed[[col]], na.rm = TRUE)
}

# Now fit the random forest model for imputed data 
rf_model <- randomForest(els_low ~ ., data = exp_low_imputed)


#real rf model

# Replace NaNs with 0 in the `els_low` variable
# Remove rows with any missing values in the dataset
exp_low_clean <- exp_low[complete.cases(exp_low), ]

rf_model_nonimpute <- randomForest(els_low ~ ., data = exp_low_clean)

# Make predictions
rf_predictions <- predict(rf_model_nonimpute, data=exp_low_clean)

# Evaluate performance (e.g., RMSE)
rf_rmse <- sqrt(mean((els_low - rf_predictions)^2))
print(rf_rmse)

importance(rf_model_nonimpute)

varImpPlot(rf_model_nonimpute)

# Get the variable importance
importance_values <- rf_model_nonimpute$importance

# Display importance values
print(importance_values)

# Permutation-based significance test for each feature - EDIT THIS CODE NEED TO RUN PERMUTATION TESTS FOR EACH FACTOR 
n_permutations <- 1000
importance_permuted <- matrix(NA, ncol = ncol(exp_low_clean) - 1, nrow = n_permutations)

for (i in 1:n_permutations) {
  # Randomize each feature and measure importance
  perm_data <- exp_low_clean
  for (j in 2:ncol(exp_low_clean)) {  # Skipping the target column (1st column)
    perm_data[, j] <- sample(perm_data[, j])
  }
  perm_model <- randomForest(els_low ~ ., data = perm_data)
  importance_permuted[i, ] <- perm_model$importance[, 1]
}

# Compute p-values for each feature's importance
p_values <- apply(importance_permuted, 2, function(x) mean(x >= importance_values[, 1]))
print(p_values)


#for parameters
rf_model_nonimpute$mtry
rf_model_nonimpute$oob.error
rf_model_nonimpute$maxnodes


#code for feature importance plot
# Load necessary packages
library(randomForest)
library(ggplot2)

# Get the importance of each feature
importance_data <- importance(rf_model_nonimpute)

# Convert importance data to a data frame for visualization
importance_df <- data.frame(Feature = rownames(importance_data), Importance = importance_data[, 1])

# Plot feature importance using ggplot2
ggplot(importance_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +
  labs(title = "Feature Importance", x = "Feature", y = "Importance") +
  theme_minimal()

#apartial dependence plot A Partial Dependence Plot shows the effect of a single feature on the predicted outcome, holding other features constant. T
# Load necessary package for Partial Dependence Plots
library(pdp)

pdp_plot <- partial(rf_model_nonimpute, pred.var = "HT3A_map")

plot(pdp_plot, main = "Partial Dependence Plot for HT3a")

# Create a partial dependence plot for both HT3a and CHRNA5 and HT2B
pdp_plot_multi <- partial(rf_model_nonimpute, pred.var = c("HT3A_map", "CHRNA5_map", "HT2B_map"))

# Plot the PDP for both features
plot(pdp_plot_multi, main = "Partial Dependence Plot for HT3a, CHRNA5 and HT2b")



  
