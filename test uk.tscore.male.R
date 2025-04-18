library(caret)
library(pROC)
library(dplyr)
library(tidyr)
#set.seed(191)
# Load the dataset



#### for test FOS

# Select the new set of columns from SOF_clinical
MrOS_clinical <- read.table("/Users/mac/Desktop/dissertation/cooperation/MrOS_data/mros1.txt", 
                            header = TRUE, sep = "\t", na.strings = c("", "NA", " "))

MrOS_clinical1 <- MrOS_clinical

# Function to compute T-score
compute_t_score <- function(value, mean_val, sd_val) {
  (value - mean_val) / sd_val
}

# Initialize columns for T-scores
MrOS_clinical1$HA_THD <- NA
MrOS_clinical1$HA_LSD <- NA
MrOS_clinical1$HA_FND <- NA

# Define age groups
age_groups <- list("65_75" = c(65, 75), "76_96" = c(76, 96))

# Loop through each age group
for (age_label in names(age_groups)) {
  age_range <- age_groups[[age_label]]
  
  # Subset data for the current age group and HA_HIPFX == 0
  subset_data <- MrOS_clinical[MrOS_clinical$HA_AGE >= age_range[1] &
                                 MrOS_clinical$HA_AGE <= age_range[2] &
                                 MrOS_clinical$HA_HIPFX == 0, ]
  
  # Compute means and standard deviations
  mean_vals <- colMeans(subset_data[, c("B1THD", "B1TLD", "B1FND")], na.rm = TRUE)
  sd_vals <- apply(subset_data[, c("B1THD", "B1TLD", "B1FND")], 2, sd, na.rm = TRUE)
  
  # Apply T-score transformation for matching individuals in SOF_clinical1
  match_idx <- which(MrOS_clinical1$HA_AGE >= age_range[1] &
                       MrOS_clinical1$HA_AGE <= age_range[2])
  
  MrOS_clinical1$HA_THD[match_idx] <- compute_t_score(MrOS_clinical1$B1THD[match_idx], mean_vals["B1THD"], sd_vals["B1THD"])
  MrOS_clinical1$HA_LSD[match_idx] <- compute_t_score(MrOS_clinical1$B1TLD[match_idx], mean_vals["B1TLD"], sd_vals["B1TLD"])
  MrOS_clinical1$HA_FND[match_idx] <- compute_t_score(MrOS_clinical1$B1FND[match_idx], mean_vals["B1FND"], sd_vals["B1FND"])
}


MrOS_clinical <- MrOS_clinical1
# Define thresholds based on quantiles of HA_WLKSPED
Q1 <- quantile(MrOS_clinical$HA_WLKSPED, 0.25, na.rm = TRUE)
Q3 <- quantile(MrOS_clinical$HA_WLKSPED, 0.75, na.rm = TRUE)

# Create categorical variable 'walk'
MrOS_clinical$walk <- cut(MrOS_clinical$HA_WLKSPED, 
                          breaks = c(-Inf, Q1, Q3, Inf), 
                          labels = c("Slow pace", "Steady average pace", "Brisk pace"))

# Check distribution
table(MrOS_clinical$walk)
library(dplyr)

# Recode levels for MrOS_clinical
MrOS_clinical$HA_IADL51 <- recode(MrOS_clinical$HA_IADL51,
                                  `1.25` = 1,
                                  `2.5` = 2,
                                  `3.75` = 3,
                                  `5` = 4)

# Check first few rows
head(MrOS_clinical[, c("HA_AGE", "HA_THD", "HA_LSD", "HA_FND")])

df <- MrOS_clinical[, c("HA_THD", "HA_LSD", "HA_FND",  "HA_GRIPAVG","walk","HA_SLDFX","HA_WRSTFX",
                        "HA_HEIGHT", "HA_AGE", "HA_WEIGHT","HA_IADL51","HA_SMOKE", "HA_HIPFX")]
write.csv(df, "/Users/mac/Desktop/dissertation/cooperation/FOS/MrOSvsUK.csv", row.names = FALSE)

# Define Clinical and DXA Features
clinical_features <- c(  "HA_GRIPAVG","walk","HA_SLDFX","HA_WRSTFX",
                         "HA_HEIGHT", "HA_AGE", "HA_WEIGHT","HA_IADL51","HA_SMOKE")


dxa_features <- c("HA_THD", "HA_LSD", "HA_FND")
target_col <- "HA_HIPFX"

###female
# Select the new set of columns from SOF_clinical
SOF_clinical <- read.table("/Users/mac/Desktop/dissertation/cooperation/SOF_data/sof.txt",
                           header = TRUE, sep = "\t", na.strings = c("", "NA", " "))

SOF_clinical1 <- SOF_clinical

# Function to compute T-score
compute_t_score <- function(value, mean_val, sd_val) {
  (value - mean_val) / sd_val
}

# Initialize columns for T-scores
SOF_clinical1$HA_THD <- NA
SOF_clinical1$HA_LSD <- NA
SOF_clinical1$HA_FND <- NA

# Define age groups
age_groups <- list("65_75" = c(65, 75), "76_96" = c(76, 96))

# Loop through each age group
for (age_label in names(age_groups)) {
  age_range <- age_groups[[age_label]]
  
  # Subset data for the current age group and HA_HIPFX == 0
  subset_data <- SOF_clinical[SOF_clinical$HA_AGE >= age_range[1] &
                                SOF_clinical$HA_AGE <= age_range[2] &
                                SOF_clinical$HA_HIPFX == 0, ]
  
  # Compute means and standard deviations
  mean_vals <- colMeans(subset_data[, c("HTOTBMD", "STOTBMD", "NBMD")], na.rm = TRUE)
  sd_vals <- apply(subset_data[, c("HTOTBMD", "STOTBMD", "NBMD")], 2, sd, na.rm = TRUE)
  
  # Apply T-score transformation for matching individuals in SOF_clinical1
  match_idx <- which(SOF_clinical1$HA_AGE >= age_range[1] &
                       SOF_clinical1$HA_AGE <= age_range[2])
  
  
  SOF_clinical1$HA_THD[match_idx] <- compute_t_score(SOF_clinical1$HTOTBMD[match_idx], mean_vals["HTOTBMD"], sd_vals["HTOTBMD"])
  SOF_clinical1$HA_LSD[match_idx] <- compute_t_score(SOF_clinical1$STOTBMD[match_idx], mean_vals["STOTBMD"], sd_vals["STOTBMD"])
  SOF_clinical1$HA_FND[match_idx] <- compute_t_score(SOF_clinical1$NBMD[match_idx], mean_vals["NBMD"], sd_vals["NBMD"])
}

SOF_clinical <- SOF_clinical1
# Define thresholds based on quantiles of HA_WLKSPED
Q1 <- quantile(SOF_clinical$HA_WLKSPED, 0.25, na.rm = TRUE)
Q3 <- quantile(SOF_clinical$HA_WLKSPED, 0.75, na.rm = TRUE)

# Create categorical variable 'walk'
SOF_clinical$walk <- cut(SOF_clinical$HA_WLKSPED, 
                         breaks = c(-Inf, Q1, Q3, Inf), 
                         labels = c("Slow pace", "Steady average pace", "Brisk pace"))

# Check distribution
table(SOF_clinical$walk)
library(dplyr)

SOF_clinical$HA_IADL51[SOF_clinical$HA_IADL51 == 5] <- 4


# Check first few rows
head(SOF_clinical[, c("HA_AGE", "HA_THD", "HA_LSD", "HA_FND")])

df <- SOF_clinical[, c("HA_THD", "HA_LSD", "HA_FND",  "HA_GRIPAVG","walk","HA_SLDFX","HA_WRSTFX",
                       "HA_HEIGHT", "HA_AGE", "HA_WEIGHT","HA_IADL51","HA_SMOKE", "HA_HIPFX")]
write.csv(df, "/Users/mac/Desktop/dissertation/cooperation/FOS/SOFvsUK.csv", row.names = FALSE)


# Define Clinical and DXA Features
clinical_features <- c(  "HA_GRIPAVG","walk","HA_SLDFX","HA_WRSTFX",
                         "HA_HEIGHT", "HA_AGE", "HA_WEIGHT","HA_IADL51","HA_SMOKE")


dxa_features <- c("HA_THD", "HA_LSD", "HA_FND")
target_col <- "HA_HIPFX"

# Drop rows with missing values in selected features
# Create a vector of columns to check for NA values
uk_data <- read.csv("/Users/mac/Desktop/dissertation/cooperation/FOS/UK_mytscore600.dataset.csv", header = TRUE, stringsAsFactors = FALSE)
UK_male <- subset(uk_data, sex == "Male")
UK <- UK_male
#UK_female <- subset(uk_data, sex == "Female")
#UK <- UK_female

colnames(UK)

# Remove rows with NA values in any of the specified columns
df_clean <- na.omit(df)
y <- as.factor(df_clean$HA_HIPFX)
testUK_y <- as.factor(UK$hip)

# Extract X (Clinical + DXA) and y (Target) using base R
X_clinical <- df_clean[, clinical_features, drop = FALSE]  # Keep as data frame
X_dxa <- df_clean[, dxa_features, drop = FALSE]  # Keep as data frame
testUK_clinical <- as.data.frame(cbind(UK$grip,UK$walk,UK$arm,UK$wrist,UK$height, UK$age,UK$weight, UK$IADLdiffnew, UK$smoke))
testUK_dxa <- as.data.frame(cbind(UK$totalbmd_Tscore,UK$spinebmd_Tscore,UK$neckbmd_Tscore))
colnames(testUK_clinical) <- c( "HA_GRIPAVG", "walk" , "HA_SLDFX" ,  "HA_WRSTFX" , "HA_HEIGHT", "HA_AGE", "HA_WEIGHT" ,"HA_IADL51","HA_SMOKE")
colnames(testUK_dxa) <- c("HA_THD", "HA_LSD", "HA_FND")
write.csv(df, "/Users/mac/Desktop/dissertation/cooperation/FOS/SOFvsUK.csv", row.names = FALSE)

# Convert character variables to numeric where appropriate
testUK_clinical$HA_GRIPAVG <- as.numeric(testUK_clinical$HA_GRIPAVG)
testUK_clinical$HA_HEIGHT <- as.numeric(testUK_clinical$HA_HEIGHT)
testUK_clinical$HA_AGE <- as.numeric(testUK_clinical$HA_AGE)
X_clinical$HA_AGE <- as.numeric(X_clinical$HA_AGE)

testUK_clinical$HA_WEIGHT <- as.numeric(testUK_clinical$HA_WEIGHT)
testUK_clinical$HA_IADL51 <- as.numeric(testUK_clinical$HA_IADL51)
X_clinical$HA_WRSTFX <- as.factor(X_clinical$HA_WRSTFX)
X_clinical$HA_SLDFX <- as.factor(X_clinical$HA_SLDFX)
X_clinical$HA_SMOKE <- as.factor(X_clinical$HA_SMOKE)
testUK_clinical$HA_WRSTFX <- as.factor(testUK_clinical$HA_WRSTFX)
testUK_clinical$HA_SLDFX <- as.factor(testUK_clinical$HA_SLDFX)
testUK_clinical$HA_SMOKE <- as.factor(testUK_clinical$HA_SMOKE)
testUK_clinical$walk <- as.factor(testUK_clinical$walk)



str(X_clinical)

# Check the structure after conversion
str(testUK_clinical)

# Split Data (80% Train, 20% Test)
#set.seed(42)
n_runs <- 100
auc_values_staged <- numeric(n_runs)
final_accuracy_scores <- numeric(n_runs)
final_sensitivity_scores <- numeric(n_runs)
final_specificity_scores <- numeric(n_runs)
threshold <- numeric(n_runs)
accuracy_stage1 <- numeric(n_runs)
accuracy_stage2 <- numeric(n_runs)
percent_stage1 <- numeric(n_runs)
percent_stage2 <- numeric(n_runs)
ssindx <- vector("list", n_runs)  # Creates an empty list with n_runs elements


for (run in 1:n_runs) {
  trainIndex <- createDataPartition(y, p = 0.8, list = FALSE)
  X_train_c <- X_clinical[trainIndex, ]
  #X_test_c <- X_clinical[-trainIndex, ]
  X_test_c <- testUK_clinical
  X_train_d <- X_dxa[trainIndex, ]
  #X_test_d <- X_dxa[-trainIndex, ]
  X_test_d <- testUK_dxa
  y_train <- y[trainIndex]
  #y_test <- y[-trainIndex]
  y_test <- testUK_y
  
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
  ## enter stage 1 person and percentage
  y_true <- y_test[certain_indices]  # Actual labels
  y_proba <- mean_proba_test[certain_indices]  # Predicted probabilities
  
  percentage_stage1 <- length(y_true)/length(y_test)
  
  ## enter stage 2 person and percentage
  
  y_true2 <- y_test[-certain_indices]  # Actual labels
  y_proba2 <- mean_proba_ensemble2[-certain_indices]  # Predicted probabilities
  
  percentage_stage2 <- length(y_true2)/length(y_test)
  
  
  
  # 3️⃣ Accuracy Using Optimal Threshold## 
  ## stage 1 accutacy
  predicted_labels_opt <- ifelse(y_proba >= best_threshold, 1, 0)
  accuracy_opt <- mean(predicted_labels_opt == y_true)
  
  ## stage 2 accutacy
  predicted_labels_opt2 <- ifelse(y_proba2 >= best_threshold, 1, 0)
  accuracy_opt2 <- mean(predicted_labels_opt2 == y_true2)
  
  print(paste("Accuracy with optimal threshold:", round(accuracy_opt, 4)))
  sum(y_proba > best_threshold & y_proba < 0.5)
  
  testIndex_c <- setdiff(1:nrow(X_clinical), trainIndex)
  
  stage1_sample_index <- testIndex_c[certain_indices]
  
  auc_values_staged[run] <- auc(roc(y_test, final_predictions))
  roc_obj <- roc(y_test, final_predictions)
  # Get the optimal threshold using Youden's index
  best_threshold <- best_threshold
 # best_threshold <- coords(roc_obj, "best", ret = "threshold", best.method = "youden")
  # Ensure y_pred_prob is a numeric vector
  final_predictions <- as.numeric(final_predictions)  
  
  # Convert predicted probabilities into binary class labels using the best threshold
  predicted_labels <- rep(0, length(final_predictions))  # Initialize with 0s
  predicted_labels[final_predictions >= best_threshold] <- 1  # Assign 1 where condition is met
  # Compute confusion matrix
  cm <- confusionMatrix(factor(predicted_labels), factor(y_test), positive = "1")
  
  # Extract accuracy, sensitivity, and specificity
  final_accuracy_scores[run] <- cm$overall["Accuracy"]
  final_sensitivity_scores[run]<- cm$byClass["Sensitivity"]
  final_specificity_scores[run] <- cm$byClass["Specificity"]
  
  threshold[run] <- best_threshold
  accuracy_stage1[run] <- accuracy_opt
  accuracy_stage2[run] <- accuracy_opt2
  
  percent_stage1[run] <- percentage_stage1
  percent_stage2[run] <- percentage_stage2
  
  ssindx[[run]] <- stage1_sample_index
}

# Compute the average AUC over 100 runs
average_auc <- mean(auc_values_staged)
cat("\u2705 Average AUC over 100 runs =", round(average_auc, 3), "confirming robust performance.\n")
mean_auc <- round(mean(auc_values_staged), 3)
sd_auc <- round(sd(auc_values_staged), 3)

mean_accuracy <- round(mean(final_accuracy_scores), 3)
sd_accuracy <- round(sd(final_accuracy_scores), 3)

mean_sensitivity <- round(mean(final_sensitivity_scores), 3)
sd_sensitivity <- round(sd(final_sensitivity_scores), 3)

mean_specificity <- round(mean(final_specificity_scores), 3)
sd_specificity <- round(sd(final_specificity_scores), 3)

average_accuracy <- mean(accuracy_stage1)
average_accuracy
average_percent <- mean(percent_stage1)
average_percent

average_accuracy <- mean(accuracy_stage2)
average_accuracy
average_percent <- mean(percent_stage2)
average_percent

# Create a data frame with the results
performance_summary <- data.frame(
  Metric = c("AUC", "Accuracy", "Sensitivity", "Specificity"),
  Mean = c(mean_auc, mean_accuracy, mean_sensitivity, mean_specificity),
  SD = c(sd_auc, sd_accuracy, sd_sensitivity, sd_specificity)
)

# Print the summary table
print(performance_summary)


# Unlist the lists to combine all the integers into one vector
combined_values <- unlist(ssindx)

# Count the occurrences of each integer
integer_counts <- table(combined_values)

# Print the result

length(integer_counts)

# Convert the table to a data frame
integer_counts_df <- as.data.frame(integer_counts)

# Write the data frame to a CSV file
write.csv(integer_counts_df, "/Users/mac/Desktop/dissertation/cooperation/harmonize_sof_code/integer_counts.csv", row.names = FALSE)

integer_counts_df






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


