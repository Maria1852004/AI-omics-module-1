#### PRACTICE EXERCISES ####

# ========================================================================
# 1. Check Cholesterol level (using if)
# ========================================================================

cholesterol <- 230  

if (cholesterol > 240) {
  print("High Cholesterol")
}

# ========================================================================
# 2. Blood Pressure Status (using if...else)
# ========================================================================

Systolic_bp <- 130  

if (Systolic_bp < 120) {
  print("Blood Pressure is normal")
} else {
  print("Blood Pressure is high")
}

# ========================================================================
# 3. Automating Data Type Conversion with for loop
# ========================================================================

# Create example datasets
data <- data.frame(
  patient_id = 1:6,
  age = c(45, 30, 55, 26, 40, 38),
  gender = c("Male","Female","Female","Male","Female","Male"),
  diagnosis = c("cancer","normal","normal","cancer","normal","cancer"),
  BMI = c(27.5, 22.4, 31.1, 24.3, 29.0, 26.7),
  smoking_status = c("Yes","No","No","Yes","No","Yes"),
  stringsAsFactors = FALSE
)

metadata <- data.frame(
  column = names(data),
  type = sapply(data, class),
  description = c("ID","Age","Gender","Diagnosis","BMI","Smoker"),
  stringsAsFactors = FALSE
)

# Create copies
data_copy <- data
metadata_copy <- metadata

# Factor columns
factor_cols <- c("gender", "smoking_status")

# Convert in data_copy
for (col in factor_cols) {
  data_copy[[col]] <- as.factor(data_copy[[col]])
}

# Convert in metadata_copy
for (col in factor_cols) {
  if (col %in% names(metadata_copy)) {
    metadata_copy[[col]] <- as.factor(metadata_copy[[col]])
  }
}

# ========================================================================
# 4. Converting Factors to Numeric Codes
# ========================================================================

binary_cols <- c("smoking_status")

for (col in binary_cols) {
  data_copy[[col]] <- ifelse(data_copy[[col]] == "Yes", 1, 0)
}

# ========================================================================
# 5. Verification
# ========================================================================

cat("\n--- Original Data Structure ---\n")
str(data)

cat("\n--- Modified Data Structure ---\n")
str(data_copy

    