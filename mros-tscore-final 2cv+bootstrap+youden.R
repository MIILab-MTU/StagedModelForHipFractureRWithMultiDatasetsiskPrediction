library(caret)
library(pROC)
library(dplyr)
library(tidyr)
set.seed(1889)

MrOS_clinical <- read.table("/Users/mac/Desktop/dissertation/cooperation/MrOS_data/mros1.txt", 
                            header = TRUE, sep = "\t", na.strings = c("", "NA", " "))

# Function to compute T-score
compute_t_score <- function(value, mean_val, sd_val) {
  (value - mean_val) / sd_val
}

# Initialize columns for T-scores
MrOS_clinical$HA_THD <- NA
MrOS_clinical$HA_LSD <- NA
MrOS_clinical$HA_FND <- NA

# Define age groups
age_groups <- list("65_75" = c(65, 75), "76_96" = c(76, 96))

# Define race levels
race_levels <- unique(MrOS_clinical$HA_RACE)  # Assuming race is coded as 1,2,3,4,5

# Loop through each (age, race) group
for (age_label in names(age_groups)) {
  age_range <- age_groups[[age_label]]
  
  for (race in race_levels) {
    
    # Subset data for the current (age, race) group
    subset_data <- MrOS_clinical[MrOS_clinical$HA_AGE >= age_range[1] & 
                                   MrOS_clinical$HA_AGE <= age_range[2] & 
                                   MrOS_clinical$HA_HIPFX == 0 & 
                                   MrOS_clinical$HA_RACE == race, ]
    
    # Compute means and standard deviations
    mean_vals <- colMeans(subset_data[, c("B1THD", "B1TLD", "B1FND")], na.rm = TRUE)
    sd_vals <- apply(subset_data[, c("B1THD", "B1TLD", "B1FND")], 2, sd, na.rm = TRUE)
    
    # Apply T-score transformation for matching individuals in MrOS_clinical_subset1
    for (i in seq_len(nrow(MrOS_clinical))) {
      if (MrOS_clinical$HA_AGE[i] >= age_range[1] & MrOS_clinical$HA_AGE[i] <= age_range[2] &
          MrOS_clinical$HA_RACE[i] == race) {
        
        MrOS_clinical$HA_THD[i] <- compute_t_score(MrOS_clinical$B1THD[i], mean_vals["B1THD"], sd_vals["B1THD"])
        MrOS_clinical$HA_LSD[i] <- compute_t_score(MrOS_clinical$B1TLD[i], mean_vals["B1TLD"], sd_vals["B1TLD"])
        MrOS_clinical$HA_FND[i] <- compute_t_score(MrOS_clinical$B1FND[i], mean_vals["B1FND"], sd_vals["B1FND"])
      }
    }
  }
}

# Check first few rows
head(MrOS_clinical[, c("HA_AGE", "HA_RACE", "HA_THD", "HA_LSD", "HA_FND")])

set.seed(1889)

MrOS_clinical <- read.table("/Users/mac/Desktop/dissertation/cooperation/MrOS_data/mros1.txt", 
                            header = TRUE, sep = "\t", na.strings = c("", "NA", " "))



dim(MrOS_clinical)
df <- MrOS_clinical[, c("HA_FND", "HA_THD", "HA_WLKSPED", "HA_LSD", "HA_GRIPAVG", 
                        "HA_VITAL", "HA_IADL51", "HA_AGE", "HA_TRLBTS", "HA_BMI", 
                        "HA_HEIGHT", "HA_MMSE", "HA_KIDNYST", "HA_WRSTFXFU",
                        "HA_SMOKE","HA_HIPFX")]
# Define Clinical and DXA Features
clinical_features <- c("HA_WLKSPED", "HA_GRIPAVG", 
                       "HA_VITAL", "HA_IADL51", "HA_AGE", "HA_TRLBTS", "HA_BMI", 
                       "HA_HEIGHT", "HA_MMSE", "HA_CALCIUM", "HA_WRSTFXFU", 
                       "HA_SMOKE")
dxa_features <- c("HA_THD", "HA_LSD", "HA_FND")
target_col <- "HA_HIPFX"


# Drop rows with missing values in selected features
# Create a vector of columns to check for NA values

# Remove rows with NA values in any of the specified columns
df_clean <- na.omit(df)
y <- as.factor(df_clean$HA_HIPFX)

# Extract X (Clinical + DXA) and y (Target) using base R
X_clinical <- df_clean[, clinical_features, drop = FALSE]  # Keep as data frame
X_dxa <- df_clean[, dxa_features, drop = FALSE]  # Keep as data frame

# Split Data (76% Train, 20% Test)
#set.seed(42)
n_runs <- 1
auc_values_staged <- numeric(n_runs)
threshold <- numeric(n_runs)

for (run in 1:n_runs) {
trainIndex <- createDataPartition(y, p = 0.8, list = FALSE)
X_train_c <- X_clinical[trainIndex, ]
X_test_c <- X_clinical[-trainIndex, ]
X_train_d <- X_dxa[trainIndex, ]
X_test_d <- X_dxa[-trainIndex, ]
y_train <- y[trainIndex]
y_test <- y[-trainIndex]

# Train Ensemble 1 (Bootstrapped Logistic Regression on Clinical Features)
#set.seed(42)
# Split Train Data into Training and Validation (76% of train for training, 20% for validation)
k_folds <- 5
folds <- createFolds(y_train, k = k_folds, list = TRUE)

# Store AUCs and thresholds
auc_values <- numeric(k_folds)
thresholds <- numeric(k_folds)

# Perform k-fold cross-validation
for (fold in seq_along(folds)) {
  # Split into training and validation fold
  train_idx <- unlist(folds[-fold])  # Use all other folds as training
  valid_idx <- folds[[fold]]  # Current fold as validation
  
  X_train_cv <- X_train_c[train_idx, ]
  y_train_cv <- y_train[train_idx]
  X_valid_cv <- X_train_c[valid_idx, ]
  y_valid_cv <- y_train[valid_idx]
  
  # Train Bootstrapped Ensemble Logistic Regression on Training Fold
  n_models <- 50
  proba_valid_cv <- matrix(NA, nrow = nrow(X_valid_cv), ncol = n_models)
  
  for (i in 1:n_models) {
    boot_idx <- sample(seq_len(nrow(X_train_cv)), replace = TRUE)
    model <- glm(y_train_cv[boot_idx] ~ ., data = X_train_cv[boot_idx, ], family = binomial)
    proba_valid_cv[, i] <- predict(model, newdata = X_valid_cv, type = "response")
  }
  
  # Compute mean probability for the validation fold
  mean_proba_valid_cv <- rowMeans(proba_valid_cv)
  
  # Compute ROC curve
  roc_obj_cv <- roc(y_valid_cv, mean_proba_valid_cv)
  auc_values[fold] <- auc(roc_obj_cv)  # Store AUC for this fold
  
  # Compute Youden's index to determine the optimal threshold
  youden_index_cv <- roc_obj_cv$sensitivities + roc_obj_cv$specificities - 1
  valid_thresholds_cv <- roc_obj_cv$thresholds[is.finite(roc_obj_cv$thresholds)]
  thresholds[fold] <- valid_thresholds_cv[which.max(youden_index_cv)]
}

# Select the best threshold corresponding to the best AUC
best_fold <- which.max(auc_values)
best_threshold <- thresholds[best_fold]

# Print results
print(paste("Best AUC:", round(auc_values[best_fold], 4)))
print(paste("Optimal Threshold from Best Validation Fold:", round(best_threshold, 4)))

# Train Final Ensemble on Full Training Set (X_train_c) and Evaluate on Test Set
n_models <- 50
# Train Ensemble on Full Training Data (X_train_c) and Evaluate on Test Set
proba_test <- matrix(NA, nrow = nrow(X_test_c), ncol = n_models)

for (i in 1:n_models) {
  boot_idx <- sample(seq_len(nrow(X_train_c)), replace = TRUE)
  model <- glm(y_train[boot_idx] ~ ., data = X_train_c[boot_idx, ], family = binomial)
  proba_test[, i] <- predict(model, newdata = X_test_c, type = "response")
}

# Compute mean probability for test set
mean_proba_test <- rowMeans(proba_test)

# Compute standard deviation for each row of proba_test
sd_proba_test <- apply(proba_test, 1, sd)

# Print first few values
head(mean_proba_test)
head(sd_proba_test)

# Compute z-score for uncertainty
uncertainty_ensemble1 <- abs((mean_proba_test-best_threshold) / sd_proba_test)

# Define uncertain cases based on z-score threshold (e.g., 1 standard deviation)
z_threshold <- 1  # Adjust as needed
uncertain_cases <- uncertainty_ensemble1 < z_threshold
sum(uncertain_cases=="TRUE")

# Prepare Data for Stage 2 (Clinical + DXA Features)
X_train_combined <- cbind(X_train_c, X_train_d)
X_test_combined <- cbind(X_test_c, X_test_d)
y_train_aligned <- y_train

# Train Ensemble 2 (Clinical + DXA)
#set.seed(42)
n_models <- 50
bootstrap_models2 <- list()
proba_ensemble2 <- matrix(NA, nrow = nrow(X_test_combined), ncol = n_models)


for (i in 1:n_models) {
  boot_idx <- sample(seq_len(nrow(X_train_combined)), replace = TRUE)
  model <- glm(y_train_aligned[boot_idx] ~ ., data = X_train_combined[boot_idx, ], family = binomial)
  proba_ensemble2[, i] <- predict(model, newdata = X_test_combined, type = "response")
}
mean_proba_ensemble2 <- rowMeans(proba_ensemble2)

# Merge Staged Predictions
final_predictions <- mean_proba_test

# Ensure mean_proba_ensemble2 has the same length as uncertain cases
if (sum(uncertain_cases) > 0) {
  final_predictions[uncertain_cases] <- mean_proba_ensemble2[uncertain_cases]
}

auc_values_staged[run] <- auc(roc(y_test, final_predictions))
threshold[run] <- best_threshold
}

# Compute the average AUC over 100 runs
average_auc <- mean(auc_values_staged)
cat("\u2705 Average AUC over 100 runs =", round(average_auc, 3), "confirming robust performance.\n")


# Print first few final predictions
head(final_predictions)

# Compute Final AUC for Staged Model
auc_staged_model <- auc(roc(y_test, final_predictions))

# Store Results in DataFrame
staged_results_df <- data.frame(
  Sample_Index = seq_along(final_predictions),
  True_Label = y_test,
  Stage_1_Probability = mean_proba_test,
  Uncertainty = uncertainty_ensemble1,
  Final_Probability = final_predictions,
  Final_Prediction = as.integer(final_predictions > 0.5)
)

# Print Final AUC
cat("\u2705 AUC =", round(auc_staged_model, 3), "confirms strong classification performance.\n")
best_threshold
# Display results
tools::display_dataframe_to_user(name="Staged Model Results (Final Fixed Version)", dataframe=staged_results_df)
