RANDOM FOREST

# Install packages 
# install.packages(c("randomForest", "ggplot2", "pdp", "oro.nifti"))

# Load necessary libraries
library(randomForest)
library(ggplot2)
library(oro.nifti)
library(pdp)

setwd("/scratch/j90161ms")

# Define patterns for the files that I want loaded in as input files
patterns <- c("^HT.*\\.nii(\\.gz)?$",   # all HT*.nii or HT*.nii.gz
              "^D.*\\.nii(\\.gz)?$",    # all D*.nii
              "^CHRM.*\\.nii(\\.gz)?$", # all CHRM*.nii
              "^CHRNA.*\\.nii(\\.gz)?$",# all CHRNA*.nii
              "^ADR.*\\.nii(\\.gz)?$")  # all ADR*.nii

# Get all files in scratch folder that match any of these patterns
scratch_dir <- "/scratch/j90161ms/RandomForest/"
nifti_files <- list.files(scratch_dir, full.names = TRUE)
nifti_files <- nifti_files[grepl(paste(patterns, collapse="|"), basename(nifti_files))]

# Check which files are included
print(nifti_files)

# 
# Load NIfTI data and construct dataframe


# Read all predictor NIfTI files and vectorise them
predictor_data <- lapply(nifti_files, function(f) {
  img <- readNIfTI(f, reorient = FALSE)
  as.vector(img)
})

# Name predictors based on file names
names(predictor_data) <- gsub("\\.nii(\\.gz)?$", "", basename(nifti_files))

# Combine into a dataframe
exp_low <- as.data.frame(predictor_data)

# Load ELS map and vectorise
els_img <- readNIfTI("/scratch/j90161ms/RandomForest/els_low.nii", reorient = FALSE)
exp_low$els_low <- as.vector(els_img)

# Replace NA values with 0 (the logic is that biologically an NA represents no gene expression)

exp_low[is.na(exp_low)] <- 0

#this is done for reproducibility; therefore, the randomness is still random but fixed for anyone trying to rerun >
set.seed(123)

# Fit random forest regression model
rf_model <- randomForest(
  els_low ~ .,
  data = exp_low,
  importance = TRUE
)

# Predicts using the dataset.
rf_predictions <- predict(rf_model, newdata = exp_low)

# Calculate RMSE
rf_rmse <- sqrt(mean((exp_low$els_low - rf_predictions)^2))
print(rf_rmse)

# Extract variable importance
importance_values <- importance(rf_model)
print(importance_values)

# Plot default random forest importance - it's default as it's pre-permutation testing
varImpPlot(rf_model)

# Permutation-based feature-wise significance testing

n_permutations <- 1000
feature_names <- setdiff(colnames(exp_low), "els_low")
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
    perm_data <- exp_low
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

# Model parameters
rf_model$mtry
rf_model$ntree
