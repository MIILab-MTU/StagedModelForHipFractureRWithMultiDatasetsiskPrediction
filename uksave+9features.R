uk_data <- read.csv("/Users/mac/Desktop/dissertation/cooperation/FOS/UK_mytscore600.dataset.csv", header = TRUE, stringsAsFactors = FALSE)
UK_male <- subset(uk_data, sex == "Male")
UK <- UK_male
#UK_female <- subset(uk_data, sex == "Female")
#UK <- UK_female

colnames(UK)

# Remove rows with NA values in any of the specified columns

testUK_y <- as.factor(UK$hip)

# Extract X (Clinical + DXA) and y (Target) using base R
testUK_clinical <- as.data.frame(cbind(UK$grip,UK$walk,UK$arm,UK$wrist,UK$height, UK$age,UK$weight, UK$IADLdiffnew, UK$smoke))
testUK_dxa <- as.data.frame(cbind(UK$totalbmd_Tscore,UK$spinebmd_Tscore,UK$neckbmd_Tscore))
colnames(testUK_clinical) <- c( "HA_GRIPAVG", "walk" , "HA_SLDFX" ,  "HA_WRSTFX" , "HA_HEIGHT", "HA_AGE", "HA_WEIGHT" ,"HA_IADL51","HA_SMOKE")
colnames(testUK_dxa) <- c("HA_THD", "HA_LSD", "HA_FND")
df <- cbind(testUK_clinical,testUK_dxa,testUK_y)
colnames(df) <- c("HA_GRIPAVG", "walk" , "HA_SLDFX" ,  "HA_WRSTFX" , "HA_HEIGHT", "HA_AGE", "HA_WEIGHT" ,"HA_IADL51","HA_SMOKE","HA_THD", "HA_LSD", "HA_FND","HA_HIPFX")

write.csv(df, "/Users/mac/Desktop/dissertation/cooperation/FOS/UK600male.csv", row.names = TRUE)
