##project

{
  #install.packages(c("caret", "corrplot", "e1071", "lattice", "AppliedPredictiveModeling"))
  #install.packages("mlbench")
  #install.packages( "earth")
  #install.packages( "kernlab")
  #install.packages( "nnet")
  #install.packages( "caret")
  #install.packages("klaR")
  #install.packages("MASS")
  #install.packages("xgboost")
  #install.packages("pls")
  library(AppliedPredictiveModeling)
  library(corrplot)
  library(e1071)
  library(caret) 
  library(VIM) 
  library(dplyr)
  library(mda)
  library(MASS)
  library(klaR)
  library(nnet)
  library(kernlab)
  library(mlbench)
  library(earth)
  library(nnet)   
  library(fastDummies)
  library(glmnet)
  library(ROCR)
  library(pROC)
  library(pls)
  library(xgboost)
}
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

# Define race levels
race_levels <- unique(SOF_clinical$HA_RACE)  # Assuming race is coded as 1,2,3,4,5

# Loop through each (age, race) group
for (age_label in names(age_groups)) {
  age_range <- age_groups[[age_label]]
  
  for (race in race_levels) {
    
    # Subset data for the current (age, race) group
    subset_data <-SOF_clinical[SOF_clinical$HA_AGE >= age_range[1] & 
                                 SOF_clinical$HA_AGE <= age_range[2] & 
                                 SOF_clinical$HA_HIPFX == 0 & 
                                 SOF_clinical$HA_RACE == race, ]
    
    # Compute means and standard deviations
    mean_vals <- colMeans(subset_data[, c("HTOTBMD", "STOTBMD", "NBMD")], na.rm = TRUE)
    sd_vals <- apply(subset_data[, c("HTOTBMD", "STOTBMD", "NBMD")], 2, sd, na.rm = TRUE)
    
    # Apply T-score transformation for matching individuals inSOF_clinical1
    for (i in seq_len(nrow(SOF_clinical1))) {
      if (SOF_clinical1$HA_AGE[i] >= age_range[1] &SOF_clinical1$HA_AGE[i] <= age_range[2] &
          SOF_clinical1$HA_RACE[i] == race) {
        
        SOF_clinical1$HA_THD[i] <- compute_t_score(SOF_clinical$HTOTBMD[i], mean_vals["HTOTBMD"], sd_vals["HTOTBMD"])
        SOF_clinical1$HA_LSD[i] <- compute_t_score(SOF_clinical$STOTBMD[i], mean_vals["STOTBMD"], sd_vals["STOTBMD"])
        SOF_clinical1$HA_FND[i] <- compute_t_score(SOF_clinical$NBMD[i], mean_vals["NBMD"], sd_vals["NBMD"])
      }
    }
  }
}

# Check first few rows
head(SOF_clinical1[, c("HA_AGE", "HA_RACE", "HA_THD", "HA_LSD", "HA_FND")])
SOF_clinical <- SOF_clinical1
ha_sof <- colnames(SOF_clinical)[grepl("^HA_", colnames(SOF_clinical))]
SOF_clinical_subset <- SOF_clinical[, ha_sof]
SOF_clinical_subset1 <- SOF_clinical_subset[, -c(1,2)]
# Remove same varaibles "HA_PHYS"+"HA_IADL"="HA_IADL51"
SOF_clinical_subset1 <- SOF_clinical_subset1[, !names(SOF_clinical_subset1) %in% c("HA_PHYS", "HA_IADL")]
dim(SOF_clinical_subset1)
colnames(SOF_clinical_subset1)
SOF_clinical2 <- SOF_clinical_subset1
SOF_clinical2 <- SOF_clinical_subset1[!is.na(SOF_clinical_subset1$HA_HIPFX), ]
## if missing greater than 20%, remove those variables
SOF_clinical21 <- SOF_clinical2[, colSums(is.na(SOF_clinical2)) <= 720] 
dim(SOF_clinical21)
str(SOF_clinical21)

SOF_clinical21 <- SOF_clinical21 %>%
  mutate(
    HA_VITAL = as.factor(HA_VITAL),
    HA_RACE = as.factor(HA_RACE),
    HA_SMOKE = as.factor(HA_SMOKE),
    # HA_IADL51 = as.factor(HA_IADL51),
    HA_HIPFX = as.factor(HA_HIPFX),
    HA_WRSTFX = as.factor(HA_WRSTFX),
    HA_SLDFX = as.factor(HA_SLDFX),
    HA_KIDNYST = as.factor(HA_KIDNYST),
    HA_CHAIRSTD = as.factor(HA_CHAIRSTD),
  )

#stage1 <- SOF_clinical21[, !names(SOF_clinical21) %in% c("HA_THD", "HA_LSD", "HA_FND","HA_HIPFXFU","HA_WRSTFX", "HA_SLDFX","HA_WRSTFXFU","HA_SLDFXFU","HA_CALCIUM","HA_CROFTOA","HA_VITAL")]
#stage <- stage1
stage2 <- SOF_clinical21[, !names(SOF_clinical21) %in% c("HA_VITAL")]
stage <- stage2



cat_var <- stage[, sapply(stage, is.factor)]
# Identify columns with near-zero variance
degen_index <- nearZeroVar(cat_var)

# Get the names of the columns with near-zero variance
degen_names <- colnames(cat_var)[degen_index]

# Remove degenerate columns from the dataset
df <- stage[, !names(stage) %in% degen_names]

#con_var <- df[, sapply(df, is.numeric)]
# Identify numeric (continuous) variables
con_var <- names(df)[sapply(df, is.numeric)]

# Exclude specific variables
exclude_vars <- c("HA_THD", "HA_LSD", "HA_FND")
con_var <- setdiff(con_var, exclude_vars)

# Extract the filtered continuous variables
df_continuous <- df[, con_var]

# Display the first few rows
head(df_continuous)
# Compute the correlation matrix
cor_matrix <- cor(df_continuous, use = "pairwise.complete.obs")

# Find the indices of highly correlated variables (absolute correlation > 0.75)
high_cor <- findCorrelation(cor_matrix, cutoff = 0.83)
highcor_names <- colnames(df_continuous)[high_cor]


# Remove highcor columns from the dataset
dfs <- df[, !names(df) %in% highcor_names]

library(caret)

# Extract only numeric variables
dfs_numeric <- dfs[sapply(dfs, is.numeric)]

# Apply spatial sign transformation
dfs_spatial_sign <- spatialSign(dfs_numeric)

# Check transformed data
head(dfs_spatial_sign)

dfs_categorical <- dfs[sapply(dfs, is.factor)]  # Extract categorical variables

# Combine transformed numerical data with categorical variables
dfs_final <- cbind(dfs_spatial_sign, dfs_categorical)

# Check final dataset
str(dfs_final)

stage_clean <- as.data.frame(na.omit(dfs_final))
dim(stage_clean)
# Load necessary package
library(VIM)

library(caret)
library(randomForest)
# Example: Define predictor variables (X) and outcome variable (y)
X <- stage_clean[, -which(names(stage_clean) == "HA_HIPFX")]
y <- stage_clean$HA_HIPFX  # Assuming it's binary (0/1)
y <- as.factor(y)  # Convert to factor for classification
# Define the control function
#control <- rfeControl(functions = rfFuncs,  # Use random forest for feature selection
#                    method = "cv",        # Cross-validation
#                     number = 10)          # 10-fold CV
#rfe_result <- rfe(X, y,
#                sizes = c(ncol(X)),  # Include all features
#                 rfeControl = control)

# Print selected features
#print(rfe_result)
# Get the most important selected features
#optVariables <- head(rfe_result$optVariables,15)
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


#ensemble2.optVariables <- c("HA_FND", "HA_GRIPAVG", "HA_THD", "HA_WLKSPED", "HA_LSD", "HA_BMI",
 #                           "HA_MMSE", "HA_AGE", "HA_CALCIUM", "HA_TRLBTS", "HA_HEIGHT", "HA_IADL51",
  #                          "HA_SMOKE", "HA_KIDNYST","HA_HIPFXFU","HA_HIPFX")


ensemble1.optVariables <- c("HA_BMI","HA_AGE","HA_TRLBTS","HA_MMSE", "HA_IADL51","HA_HEIGHT","HA_WLKSPED","HA_HIPFX")
optVariables <- ensemble1.optVariables
stage_clean <- na.omit(MrOS_clinical[,optVariables])
X <- stage_clean[, -which(names(stage_clean) == "HA_HIPFX")]
y <- stage_clean$HA_HIPFX  # Assuming it's binary (0/1)
y <- as.factor(y) 
#ensemble1.optVariables <- c("HA_BMI","HA_GRIPAVG","HA_WLKSPED","HA_MMSE","HA_AGE")
#optVariables <- ensemble1.optVariables

## feature selection

# Load required libraries
library(caret) 
library(pROC)
# Set up parameters
n_runs <- 100  # Run the entire process 100 times
final_auc_scores <- numeric(n_runs)  # Store final AUC scores for each run
final_accuracy_scores <- numeric(n_runs)
final_sensitivity_scores <- numeric(n_runs)
final_specificity_scores <- numeric(n_runs)

for (run in 1:n_runs) {
  
  # Define predictor variables (X) and response variable (y)
  #X <- stage_clean[,optVariables]  # Top selected features from previous feature selection
  
  # Convert y to a factor (for logistic regression)
  #y <- as.factor(y)
  
  # Set up parameters
  n_models <- 100  # Number of bootstrap samples (ensemble size)
  n_folds <- 5    # Number of cross-validation folds
  
  # Store AUC scores for each model
  auc_scores <- numeric(n_models)
  accuracy_scores <- numeric(n_models)
  sensitivity_scores <- numeric(n_models)
  specificity_scores <- numeric(n_models)
  
  # Perform Bootstrap Sampling and Logistic Regression Training
  #set.seed(42)  # For reproducibility
  for (i in 1:n_models) {
    
    # Bootstrap sample with replacement
    boot_idx <- sample(1:nrow(X), replace = TRUE)
    X_boot <- X[boot_idx, ]
    y_boot <- y[boot_idx]
    
    # Perform stratified cross-validation
    cv_folds <- createFolds(y_boot, k = n_folds, list = TRUE)
    auc_values <- numeric(n_folds)
    accuracy_values <- numeric(n_folds)
    sensitivity_values <- numeric(n_folds)
    specificity_values <- numeric(n_folds)
    
    for (j in 1:n_folds) {
      # Split data into training and validation sets
      train_idx <- unlist(cv_folds[-j])
      test_idx <- unlist(cv_folds[j])
      
      X_train <- X_boot[train_idx, ]
      y_train <- y_boot[train_idx]
      X_test <- X_boot[test_idx, ]
      y_test <- y_boot[test_idx]
      
      # Train Logistic Regression Model
      model <- glm(y_train ~ ., data = data.frame(y_train, X_train), family = "binomial")
      
      # Get predictions (probabilities)
      y_pred_prob <- predict(model, newdata = data.frame(X_test), type = "response")
      
      # Compute AUC
      roc_obj <- roc(y_test, y_pred_prob, quiet = TRUE)
      auc_values[j] <- roc_obj$auc
      # Get the optimal threshold using Youden's index
      best_threshold <- coords(roc_obj, "best", ret = "threshold", best.method = "youden")
      # Ensure y_pred_prob is a numeric vector
      y_pred_prob <- as.numeric(y_pred_prob)  
      
      # Convert predicted probabilities into binary class labels using the best threshold
      predicted_labels <- rep(0, length(y_pred_prob))  # Initialize with 0s
      predicted_labels[y_pred_prob >= best_threshold] <- 1  # Assign 1 where condition is met
      # Compute confusion matrix
      cm <- confusionMatrix(factor(predicted_labels), factor(y_test), positive = "1")
      
      # Extract accuracy, sensitivity, and specificity
      accuracy_values[j] <- cm$overall["Accuracy"]
      sensitivity_values[j] <- cm$byClass["Sensitivity"]
      specificity_values[j] <- cm$byClass["Specificity"]
    }
    
    # Store the mean AUC across folds for this bootstrap sample
    auc_scores[i] <- mean(auc_values)
    accuracy_scores[i] <- mean(accuracy_values)
    sensitivity_scores[i] <- mean(sensitivity_values)
    specificity_scores[i] <- mean(specificity_values)
  }
  
  # Compute the final average AUC for this run
  final_auc_scores[run] <- mean(auc_scores)
  final_accuracy_scores[run] <- mean(accuracy_scores)
  final_sensitivity_scores[run] <- mean(sensitivity_scores)
  final_specificity_scores[run] <- mean(specificity_scores)
  
}

# Compute overall statistics for AUC, Accuracy, Sensitivity, and Specificity
mean_auc <- round(mean(final_auc_scores), 3)
sd_auc <- round(sd(final_auc_scores), 3)

mean_accuracy <- round(mean(final_accuracy_scores), 3)
sd_accuracy <- round(sd(final_accuracy_scores), 3)

mean_sensitivity <- round(mean(final_sensitivity_scores), 3)
sd_sensitivity <- round(sd(final_sensitivity_scores), 3)

mean_specificity <- round(mean(final_specificity_scores), 3)
sd_specificity <- round(sd(final_specificity_scores), 3)

# Create a data frame with the results
performance_summary <- data.frame(
  Metric = c("AUC", "Accuracy", "Sensitivity", "Specificity"),
  Mean = c(mean_auc, mean_accuracy, mean_sensitivity, mean_specificity),
  SD = c(sd_auc, sd_accuracy, sd_sensitivity, sd_specificity)
)

# Print the summary table
print(performance_summary)

