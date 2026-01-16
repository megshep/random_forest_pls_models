#RANDOM FOREST

# Install packages 
# install.packages(c("randomForest", "ggplot2", "pdp"))

# Load necessary libraries
library(randomForest)
library(ggplot2)
library(pdp)

#this is done for reproducibility; therefore, the randomness is still random but fixed for anyone trying to rerun this analysis
set.seed(123)

# Remove rows with missing data (complete cases only)
exp_low_clean <- exp_low[complete.cases(exp_low), ]

# Fit random forest regression model
rf_model <- randomForest(
  els_low ~ .,
  data = exp_low_clean,
  importance = TRUE
)

# Predicts using the dataset.
rf_predictions <- predict(rf_model, newdata = exp_low_clean)

# Calculate RMSE
rf_rmse <- sqrt(mean((exp_low_clean$els_low - rf_predictions)^2))
print(rf_rmse)

# Extract variable importance
importance_values <- importance(rf_model)
print(importance_values)

# Plot default random forest importance - it's default as it's pre-permutation testing
varImpPlot(rf_model)

# Permutation-based feature-wise significance testing
n_permutations <- 1000
feature_names <- setdiff(colnames(exp_low_clean), "els_low")
obs_importance <- importance_values[, 1]

#create a matrix to store the permuted importances
perm_importance <- matrix(
  NA,
  nrow = n_permutations,
  ncol = length(feature_names),
  dimnames = list(NULL, feature_names)
)

#loop over features (predictors) and permutations, generates the 1000 null models for each predictor
for (f in feature_names) {
  for (i in 1:n_permutations) {
#shuffles the features
    perm_data <- exp_low_clean
    perm_data[[f]] <- sample(perm_data[[f]])


#refits the random forest model on the permuted data
    perm_model <- randomForest(
      els_low ~ .,
      data = perm_data,
      importance = TRUE
    )
    #records the permuted importance (the mean decrease in accuracy)
    perm_importance[i, f] <- importance(perm_model)[f, 1]
  }
}
#calculates p values
p_values <- sapply(feature_names, function(f) {
  mean(perm_importance[, f] >= obs_importance[f])
})
#creates a dataframe that includes the name of the predictor, the observed importance and the p value 
importance_results <- data.frame(
  Feature = feature_names,
  Importance = obs_importance[feature_names],
  P_value = p_values
)

importance_results <- importance_results[order(-importance_results$Importance), ]
print(importance_results)

# Plot feature importance using ggplot2
ggplot(importance_results,
       aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +
  labs(
    title = "Random Forest Feature Importance",
    x = "Feature",
    y = "Mean Decrease in Accuracy"
  ) +
  theme_minimal()

# Partial dependence plot for HT3A_map
pdp_ht3a <- partial(
  rf_model,
  pred.var = "HT3A_map",
  train = exp_low_clean
)

plot(pdp_ht3a, main = "Partial Dependence Plot: HT3A_map")

# Partial dependence plot for multiple predictors
pdp_multi <- partial(
  rf_model,
  pred.var = c("HT3A_map", "CHRNA5_map", "HT2B_map"),
  train = exp_low_clean
)

plot(
  pdp_multi,
  main = "Partial Dependence Plot: HT3A_map, CHRNA5_map, HT2B_map"
)

# Model parameters
rf_model$mtry
rf_model$ntree
