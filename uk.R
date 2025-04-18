UK <- read.csv("/Users/mac/Desktop/dissertation/cooperation/ukdata.csv", header = TRUE, stringsAsFactors = FALSE)
str(UK)
## mros use right hip: MrOS DXA Quality Assurance Manual for Hologic QDR-4500 Bone Densitometers
## sof use left hip:https://pmc.ncbi.nlm.nih.gov/articles/PMC4388249/
SOF_clinical <- read.table("/Users/mac/Desktop/dissertation/cooperation/SOF_data/sof.txt",
                           header = TRUE, sep = "\t", na.strings = c("", "NA", " "))

MrOS_clinical <- read.table("/Users/mac/Desktop/dissertation/cooperation/MrOS_data/mros1.txt",
                            header = TRUE, sep = "\t", na.strings = c("", "NA", " "))
summary(SOF_clinical$HA_WLKSPED)
summary(MrOS_clinical$HA_WLKSPED)
colnames(UK) <- c("ID", "age", "totalbmd", "smoke", "spinebmd", "leftgrip", "rightgrip", 
                  "walk", "sex", "height", "fracture", "trailtest", "totalbmdtscore",  
                  "rightneckbmd", "leftneckbmd", "BMI", "IADdifficulty",  
                  "leftneckbmdtscore", "rightneckbmdtscore")

# Remove the least frequent level for 'walk'
UK <- UK[UK$walk != "None of the above", ]

# Remove the least frequent level for 'IADdifficulty'
UK <- UK[UK$IADdifficulty != "Do not know", ]

# Remove the least frequent level for 'smoke'
UK <- UK[UK$smoke != "Prefer not to answer", ]

# Check the updated tables
table(UK$walk)
table(UK$IADdifficulty)
table(UK$smoke)

# Create a new variable 'hip' (for Hip fractures)
UK$hip <- ifelse(grepl("Hip", UK$fracture), 1, 0)

# Create a new variable 'arm' (for Arm fractures)
UK$arm <- ifelse(grepl("Arm", UK$fracture), 1, 0)

# Create a new variable 'wrist' (for Wrist fractures)
UK$wrist <- ifelse(grepl("Wrist", UK$fracture), 1, 0)

# Check the first few rows
head(UK[, c("fracture", "hip", "arm", "wrist")])

### Walking speed Approach 1: Based on quantiles

# categorize HA_WLKSPED into three groups:
  #•	Slow pace: Below the 1st quartile (Q1)
#•	Steady average pace: Between the 1st quartile (Q1) and the 3rd quartile (Q3)
#•	Brisk pace: Above the 3rd quartile (Q3)

# Define thresholds based on quantiles of HA_WLKSPED
Q1 <- quantile(SOF_clinical$HA_WLKSPED, 0.25, na.rm = TRUE)
Q3 <- quantile(SOF_clinical$HA_WLKSPED, 0.75, na.rm = TRUE)

# Create categorical variable 'walk'
SOF_clinical$walk <- cut(SOF_clinical$HA_WLKSPED, 
                         breaks = c(-Inf, Q1, Q3, Inf), 
                         labels = c("Slow pace", "Steady average pace", "Brisk pace"))

# Check distribution
table(SOF_clinical$walk)

# Define thresholds based on quantiles of HA_WLKSPED
Q1 <- quantile(MrOS_clinical$HA_WLKSPED, 0.25, na.rm = TRUE)
Q3 <- quantile(MrOS_clinical$HA_WLKSPED, 0.75, na.rm = TRUE)

# Create categorical variable 'walk'
MrOS_clinical$walk <- cut(MrOS_clinical$HA_WLKSPED, 
                         breaks = c(-Inf, Q1, Q3, Inf), 
                         labels = c("Slow pace", "Steady average pace", "Brisk pace"))

# Check distribution
table(MrOS_clinical$walk)

UK$newnecktscore <- ifelse(UK$sex == "Male", UK$rightneckbmdtscore, UK$leftneckbmdtscore)
UK$IADLdiffnew <- factor(UK$IADdifficulty, 
                         levels = c("No difficulty", "Mild difficulty", "Moderate difficulty", 
                                    "Severe difficulty", "Extreme difficulty / unable to do this"), 
                         labels = c(0, 1, 2, 3, 4))

## in MrOS_clinical, change level 1.25 to 1,  change level 2.5 to 2, change level 3.75 to 3, change levle 5 to 4.   in SOF_clinical, change level 5 to 4. 
# Load dplyr if not already loaded
library(dplyr)

# Recode levels for MrOS_clinical
MrOS_clinical$HA_IADL51 <- recode(MrOS_clinical$HA_IADL51,
                                  `1.25` = 1,
                                  `2.5` = 2,
                                  `3.75` = 3,
                                  `5` = 4)

# Recode level 5 to 4 for SOF_clinical
library(dplyr)

SOF_clinical$HA_IADL51[SOF_clinical$HA_IADL51 == 5] <- 4

### T score for UK
# Function to compute T-score
compute_t_score <- function(value, mean_val, sd_val) {
  (value - mean_val) / sd_val
}

# Initialize column for T-scores
UK$spinebmd_Tscore <- NA

# Define age groups
age_groups <- list("54_65" = c(54, 65), "65_70" = c(65, 70))

# Loop through each age group
for (age_label in names(age_groups)) {
  age_range <- age_groups[[age_label]]
  
  # Subset data for the current age group and HA_HIPFX == 0
  subset_data <- UK[UK$age >= age_range[1] &
                      UK$age <= age_range[2] &
                      UK$hip == 0, ]
  
  # Compute mean and standard deviation for spinebmd
  mean_val <- mean(subset_data$spinebmd, na.rm = TRUE)
  sd_val <- sd(subset_data$spinebmd, na.rm = TRUE)
  
  # Apply T-score transformation for matching individuals in the UK dataset
  match_idx <- which(UK$age >= age_range[1] &
                       UK$age <= age_range[2])
  
  UK$spinebmd_Tscore[match_idx] <- compute_t_score(UK$spinebmd[match_idx], mean_val, sd_val)
}

# Check first few rows
head(UK[, c("age", "spinebmd_Tscore")])

colnames(UK)
colnames(MrOS_clinical)

##GSGRPAVG: The average grip strength from the 4 trials on both hands. Please note that at the baseline visit, 2 trials were completed on each hand.
UK$grip <- rowMeans(UK[, c("leftgrip", "rightgrip")], na.rm = TRUE)
UK$smoke[UK$smoke == "Never"] <- 0
UK$smoke[UK$smoke == "Previous"] <- 1
UK$smoke[UK$smoke == "Current"] <- 2

write.csv(UK, "/Users/mac/Desktop/dissertation/cooperation/FOS/UK_dataset.csv", row.names = FALSE)
