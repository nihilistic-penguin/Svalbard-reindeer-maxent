library(dplyr)

# Base output directory
base_dir <- "C:/Users/marie/Documents/R/Barrett's R/Svalbard Reindeer Maxent Paper/output"

# Define weekly folders explicitly
weekly_folders <- c("weeks_09_10", "weeks_11_12", "weeks_13_15a", "weeks_15b_16",
                    "weeks_18a_18b", "weeks_20a_20b", "weeks_21_22", "weeks_23_24",
                    "weeks_25a_25b", "weeks_27_28")

# Initialize empty list
all_ruggedness <- list()

for (week in weekly_folders) {
  
  response_file <- file.path(base_dir, week, "response_curves", "response_ruggedness.csv")
  
  if (file.exists(response_file)) {
    df <- read.csv(response_file)
    df$period <- week
    all_ruggedness[[week]] <- df
    cat("✅ Loaded:", week, "\n")
  } else {
    cat("❌ Missing:", week, "\n")
  }
}

# Combine all into one data frame
combined_ruggedness <- bind_rows(all_ruggedness)

# Save to CSV
write.csv(combined_ruggedness, file.path(base_dir, "ruggedness_response_curves_combined.csv"), row.names = FALSE)

cat("\nTotal rows:", nrow(combined_ruggedness), "\n")
cat("Periods included:", paste(unique(combined_ruggedness$period), collapse = ", "), "\n")