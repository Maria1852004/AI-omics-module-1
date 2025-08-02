# --------------------------------------------------
# Step 1: Load dataset from raw_data folder
# --------------------------------------------------

patient_data <- read.csv("raw_data/patient_info.csv")

# Check number of rows
nrow(patient_data)

# --------------------------------------------------
# Step 2: Explore and clean the data
# --------------------------------------------------

# View structure
str(patient_data)

# View unique values in smoker column
unique(patient_data$smoker)

# --------------------------------------------------
# Step 3: Clean smoker column (handle lowercase, extra spaces)
# --------------------------------------------------

# Remove extra spaces and make lowercase
patient_data$smoker <- trimws(patient_data$smoker)
patient_data$smoker <- tolower(patient_data$smoker)

# Convert to standard Yes/No format
patient_data$smoker <- ifelse(patient_data$smoker == "yes", "Yes",
                              ifelse(patient_data$smoker == "no", "No", NA))

# Convert to factor
patient_data$smoker <- as.factor(patient_data$smoker)

# --------------------------------------------------
# Step 4: Fix other data types
# --------------------------------------------------

# Convert gender to factor
patient_data$gender <- as.factor(patient_data$gender)

# (Optional) Convert age to numeric if needed
# patient_data$age <- as.numeric(patient_data$age)

# --------------------------------------------------
# Step 5: Create binary smoking status column
# --------------------------------------------------

# 1 for "Yes", 0 for "No"
patient_data$smoking_status_binary <- ifelse(patient_data$smoker == "Yes", 1, 0)

# --------------------------------------------------
# Step 6: Save cleaned dataset to clean_data folder
# --------------------------------------------------

write.csv(patient_data, "clean_data/patient_info_clean.csv", row.names = FALSE)

# --------------------------------------------------
# Done! ✅
# --------------------------------------------------







