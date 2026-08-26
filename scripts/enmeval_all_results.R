# Load required library
library(dplyr)

# Get all period folders
period_folders <- list.dirs("output", recursive = FALSE)

# Initialize empty list
all_enmeval <- list()

for (folder in period_folders) {
  period_name <- basename(folder)
  
  # Path to ENMeval results
  enmeval_file <- file.path(folder, "enmeval_all_results.csv")
  
  if (file.exists(enmeval_file)) {
    # Read file and add period column
    df <- read.csv(enmeval_file)
    df$period <- period_name
    all_enmeval[[period_name]] <- df
    cat("✅ Loaded:", period_name, "-", nrow(df), "rows\n")
  } else {
    cat("❌ Missing:", period_name, "\n")
  }
}

# Combine all into one data frame
combined_enmeval <- bind_rows(all_enmeval)

# Save to CSV
write.csv(combined_enmeval, "output/all_enmeval_results_combined.csv", row.names = FALSE)

# Check result
cat("\nTotal rows:", nrow(combined_enmeval), "\n")
cat("Periods included:", paste(unique(combined_enmeval$period), collapse = ", "), "\n")
