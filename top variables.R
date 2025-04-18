# Initialize an empty list to store rankings for each run
variable_ranks_list <- vector("list", 100)


# Run the process 100 times
for (run in 1:100) {
  
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
  
 # stage1 <- SOF_clinical21[, !names(SOF_clinical21) %in% c("HA_THD", "HA_LSD", "HA_FND","HA_VITAL","HA_WRSTFXFU","HA_SLDFXFU","HA_CALCIUM","HA_CROFTOA")]
  # stage <- stage1
  stage2 <- SOF_clinical21[, !names(SOF_clinical21) %in% c("HA_WRSTFXFU","HA_SLDFXFU","HA_VITAL")]
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
  control <- rfeControl(functions = rfFuncs,  # Use random forest for feature selection
                        method = "cv",        # Cross-validation
                        number = 10)          # 10-fold CV
  rfe_result <- rfe(X, y,
                    sizes = c(ncol(X)),  # Include all features
                    rfeControl = control)
  
  # Print selected features
  print(rfe_result)
  # Get top 15 variables and assign ranks (1 to 15)
  variable_ranks_list[[run]] <- setNames(1:15, head(rfe_result$optVariables, 15))
}

# Convert list of named vectors into a data frame
ranking_df <- do.call(rbind, lapply(variable_ranks_list, function(x) {
  data.frame(Variable = names(x), Rank = x)
}))

# Sum ranks for each variable across all runs
ranking_summary <- aggregate(Rank ~ Variable, data = ranking_df, sum)

# Calculate the proportion (average rank per run)
ranking_summary$Proportion <- ranking_summary$Rank / 100  # Since 100 runs

# Order by total rank sum
ranking_summary <- ranking_summary[order(ranking_summary$Rank), ]

# Print the final ranking
print(ranking_summary)
