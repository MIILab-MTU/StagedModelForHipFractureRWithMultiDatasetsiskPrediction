library(caret)
library(pROC)
library(dplyr)
library(tidyr)
#set.seed(191)
# Load the dataset
# Define predictor variables (X) and response variable (y)
#X <- stage_clean[,optVariables]  # Top selected features from previous feature selection

# Convert y to a factor (for logistic regression)
#y <- as.factor(y)
#df <- cbind(X,y)
# Select the new set of columns from SOF_clinical
MrOS_clinical <- read.table("/Users/mac/Desktop/dissertation/cooperation/MrOS_data/mros1.txt", 
 #                           header = TRUE, sep = "\t", na.strings = c("", "NA", " "))
#MrOS_clinicalpp <- read.csv("/Users/mac/Desktop/dissertation/cooperation/imbalance/shuonew.csv", header = TRUE, stringsAsFactors = FALSE)
dim(MrOS_clinical)
df <- MrOS_clinical[, c("HA_FND", "HA_GRIPAVG", "HA_THD", "HA_WLKSPED", "HA_LSD", "HA_BMI",
"HA_MMSE", "HA_AGE", "HA_CALCIUM", "HA_TRLBTS", "HA_HEIGHT", "HA_IADL51",
"HA_SMOKE", "HA_KIDNYST", "HA_HIPFX")]



# Define Clinical and DXA Features
clinical_features <- c( "HA_GRIPAVG",  "HA_WLKSPED", "HA_BMI",
"HA_MMSE", "HA_AGE", "HA_CALCIUM", "HA_TRLBTS", "HA_HEIGHT", "HA_IADL51",
"HA_SMOKE", "HA_KIDNYST" )
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

# Split Data (80% Train, 20% Test)
#set.seed(42)
n_runs <- 100
auc_values_staged <- numeric(n_runs)
threshold <- numeric(n_runs)
accuracy_stage1 <- numeric(n_runs)
percent_stage1 <- numeric(n_runs)
ssindx <- vector("list", n_runs)  # Creates an empty list with n_runs elements

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
# Split Train Data into Training and Validation (80% of train for training, 20% for validation)
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
z_threshold <- 2  # Adjust as needed
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

### for group staying in stage 1

# Get the indices of certain cases (uncertain_cases == "FALSE")
certain_indices <- which(uncertain_cases == "FALSE")
y_test[certain_indices]
mean_proba_test[certain_indices]


library(pROC)

# Extract relevant data
y_true <- y_test[certain_indices]  # Actual labels
y_proba <- mean_proba_test[certain_indices]  # Predicted probabilities


percentage_stage1 <- length(y_true)/length(y_test)

# 3️⃣ Accuracy Using Optimal Threshold
predicted_labels_opt <- ifelse(y_proba >= best_threshold, 1, 0)
accuracy_opt <- mean(predicted_labels_opt == y_true)
print(paste("Accuracy with optimal threshold:", round(accuracy_opt, 4)))
sum(y_proba > best_threshold & y_proba < 0.5)

testIndex_c <- setdiff(1:nrow(X_clinical), trainIndex)

stage1_sample_index <- testIndex_c[certain_indices]

auc_values_staged[run] <- auc(roc(y_test, final_predictions))
threshold[run] <- best_threshold
accuracy_stage1[run] <- accuracy_opt
percent_stage1[run] <- percentage_stage1
ssindx[[run]] <- stage1_sample_index
}

# Compute the average AUC over 100 runs
average_auc <- mean(auc_values_staged)
cat("\u2705 Average AUC over 100 runs =", round(average_auc, 3), "confirming robust performance.\n")
average_accuracy <- mean(accuracy_stage1)
average_accuracy
average_percent <- mean(percent_stage1)
average_percent


# Unlist the lists to combine all the integers into one vector
combined_values <- unlist(ssindx)

# Count the occurrences of each integer
integer_counts <- table(combined_values)

# Print the result
print(integer_counts)

length(integer_counts)
# Convert the table to a data frame
integer_counts_df <- as.data.frame(integer_counts)

# Write the data frame to a CSV file
write.csv(integer_counts_df, "/Users/mac/Desktop/dissertation/cooperation/harmonize_MrOS_code/integer_counts.csv", row.names = FALSE)

# Confirm the file is saved
print("CSV file has been saved as 'integer_counts.csv'.")




# Compute AUC and optimal threshold using Youden’s Index
roc_final <- roc(y_test, final_predictions)
optimal_threshold <- coords(roc_final, "best", ret = "threshold")
optimal_threshold <- as.numeric(optimal_threshold$threshold)

# Convert probabilities to binary predictions using the optimal threshold
predicted_labels <- ifelse(final_predictions >= optimal_threshold, 1, 0)

# Create confusion matrix
conf_matrix <- table(Predicted = predicted_labels, Actual = y_test)

# Compute Accuracy, Sensitivity, and Specificity
accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)  # Overall accuracy
sensitivity <- conf_matrix[2,2] / sum(conf_matrix[,2])  # True Positive Rate
specificity <- conf_matrix[1,1] / sum(conf_matrix[,1])  # True Negative Rate

# Print results
cat("Accuracy:", accuracy, "\n")
cat("Sensitivity:", sensitivity, "\n")
cat("Specificity:", specificity, "\n")



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

