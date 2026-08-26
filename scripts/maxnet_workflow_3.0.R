# ============================================================================
# ADAPTIVE MAXENT WORKFLOW WITH SAMPLE-SIZE-DRIVEN COMPLEXITY
# maxnet backend (no Java required)
# Expert implementation following Morales et al. 2017,
# Radosavljevic & Anderson 2014
# ============================================================================

library(maxnet)
library(raster)
library(terra)
library(pROC)
library(ecospat)

# ENMeval 2.0.5.x has a bug in block partitioning ("test_bg not found").
# We use 2.0.4 instead. To install it on a new machine:
#   download.file(
#     url      = "https://cran.r-project.org/src/contrib/Archive/ENMeval/ENMeval_2.0.4.tar.gz",
#     destfile = "C:/Users/marie/Downloads/ENMeval_2.0.4.tar.gz",
#     mode     = "wb"
#   )
#   install.packages(
#     "C:/Users/marie/Downloads/ENMeval_2.0.4.tar.gz",
#     repos = NULL,
#     type  = "source"
#   )
library(ENMeval)
library(dplyr)
library(ggplot2)
library(viridis)
library(gridExtra)
library(cowplot)
library(sf)
library(readxl)
library(tidyr)

setwd("C:/Users/marie/Documents/R/Barrett's R/Svalbard Reindeer Maxent Paper")

# ============================================================================
# 1. LOAD PREDICTOR DATA
# ============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("LOADING PREDICTOR DATA\n")
cat(rep("=", 70), "\n\n", sep = "")

predictors <- stack("data/prepared for analyses/predictors_masked2.grd")
cat("Predictor CRS:", projection(predictors), "\n")
cat("Predictor layers:", names(predictors), "\n")
cat("Raster dimensions:", nrow(predictors), "x", ncol(predictors), "cells\n")

categorical_vars <- intersect(c("landcover_fc", "landforms_fc"), names(predictors))
cat("Categorical variables found in raster:",
    paste(categorical_vars, collapse = ", "), "\n")

predictors_terra <- terra::rast(predictors)
cat("Converted to terra SpatRaster for ENMeval\n")

# ============================================================================
# 2. LOAD AND PROCESS OCCURRENCE DATA
# ============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("LOADING OCCURRENCE DATA - INDIVIDUAL REINDEER LEVEL\n")
cat(rep("=", 70), "\n\n", sep = "")

occurrence_raw <- read_excel("data/BIG_reindeer_counts_2023_upto15.7.xlsx")

cat("RAW DATA SUMMARY:\n")
cat("  Total observations (rows):", nrow(occurrence_raw), "\n")
cat("  Date range:", as.character(min(occurrence_raw$date)),
    "to", as.character(max(occurrence_raw$date)), "\n")
cat("  Total reindeer counted:", sum(occurrence_raw$total_n, na.rm = TRUE), "\n")
cat("  Weeks present:",
    paste(sort(unique(occurrence_raw$week)), collapse = ", "), "\n\n")

occurrence_clean <- occurrence_raw %>%
  dplyr::select(
    date,
    week,
    easting  = utm_easting,
    northing = utm_northing,
    count    = total_n
  ) %>%
  filter(
    !is.na(easting), !is.na(northing),
    !is.na(count), count > 0,
    !is.na(date), !is.na(week)
  ) %>%
  mutate(date = as.Date(date))

cat("AFTER CLEANING:\n")
cat("  Observations with valid data:", nrow(occurrence_clean), "\n")
cat("  Total reindeer counted:", sum(occurrence_clean$count), "\n")
cat("  Mean group size:", round(mean(occurrence_clean$count), 1), "\n")
cat("  Max group size:", max(occurrence_clean$count), "\n\n")

occurrence_individuals <- occurrence_clean %>%
  tidyr::uncount(weights = count) %>%
  mutate(occurrence = 1)

cat("AFTER EXPANDING TO INDIVIDUALS:\n")
cat("  Total individual reindeer:", nrow(occurrence_individuals), "\n")
cat("  Unique locations:",
    dplyr::n_distinct(occurrence_individuals[, c("easting", "northing")]), "\n")
cat("  Mean individuals per unique location:",
    round(nrow(occurrence_individuals) /
            dplyr::n_distinct(occurrence_individuals[, c("easting", "northing")]),
          1), "\n\n")

# ============================================================================
# 3. CHECK POINTS AGAINST RASTER EXTENT
# ============================================================================

cat("CHECKING POINTS AGAINST RASTER EXTENT:\n")

coords_sp <- sp::SpatialPoints(
  occurrence_individuals[, c("easting", "northing")],
  proj4string = sp::CRS(raster::projection(predictors))
)

test_extract  <- raster::extract(predictors[[1]], coords_sp)
points_inside <- !is.na(test_extract)

cat("  Points inside raster extent:",  sum(points_inside), "\n")
cat("  Points outside raster extent:", sum(!points_inside), "\n")

if (sum(!points_inside) > 0) {
  cat("  Removing", sum(!points_inside), "individual reindeer outside raster\n")
  occurrence_individuals <- occurrence_individuals[points_inside, ]
}

cat("  Final individual reindeer for modeling:",
    nrow(occurrence_individuals), "\n\n")

# ============================================================================
# 4. DEFINE PERIODS
# ============================================================================

periods <- list(
  weeks_09_10   = as.Date(c("2023-03-05", "2023-03-12")),
  weeks_11_12   = as.Date(c("2023-03-18", "2023-03-26")),
  weeks_13_15a  = as.Date(c("2023-04-01", "2023-04-10")),
  weeks_15b_16  = as.Date(c("2023-04-15", "2023-04-20")),
  weeks_18a_18b = as.Date(c("2023-05-01", "2023-05-07")),
  weeks_20a_20b = as.Date(c("2023-05-15", "2023-05-19")),
  weeks_21_22   = as.Date(c("2023-05-27", "2023-06-04")),
  weeks_23_24   = as.Date(c("2023-06-10", "2023-06-16")),
  weeks_25a_25b = as.Date(c("2023-06-24", "2023-07-01")),
  weeks_27_28   = as.Date(c("2023-07-09", "2023-07-15"))
)

# ============================================================================
# 5. CREATE PERIOD DATA
# ============================================================================

period_data <- lapply(periods, function(date_range) {
  occurrence_individuals %>%
    filter(date >= date_range[1], date <= date_range[2]) %>%
    dplyr::select(easting, northing) %>%
    rename(lon = easting, lat = northing) %>%
    mutate(lon = as.numeric(lon), lat = as.numeric(lat))
})

# ============================================================================
# 6. SAMPLE SIZE SUMMARY
# ============================================================================

sample_sizes <- data.frame(
  period        = names(periods),
  n_individuals = sapply(period_data, nrow),
  n_unique_cells = sapply(period_data, function(x)
    dplyr::n_distinct(x[, c("lon", "lat")])),
  start_date    = sapply(periods, function(d) as.character(d[1])),
  end_date      = sapply(periods, function(d) as.character(d[2]))
) %>%
  mutate(
    individuals_per_cell  = round(n_individuals / n_unique_cells, 1),
    sample_size_class = case_when(
      n_individuals < 50 ~ "Very small (<50)",
      n_individuals < 80 ~ "Small (50-79)",
      TRUE               ~ "Adequate (>=80)"
    )
  )

cat("\n", rep("=", 70), "\n", sep = "")
cat("SAMPLE SIZES BY PERIOD (INDIVIDUAL REINDEER)\n")
cat(rep("=", 70), "\n\n", sep = "")
print(sample_sizes)

dir.create("output", showWarnings = FALSE)
write.csv(sample_sizes, "output/sample_sizes.csv", row.names = FALSE)

# ============================================================================
# 7. ADAPTIVE COMPLEXITY CONFIGURATION
# ============================================================================

get_adaptive_tuning_params <- function(n_presence) {
  if (n_presence < 50) {
    list(fc = "L", rm = seq(3.0, 6.0, 0.5),
         description = "Linear-only, extreme regularization (n < 50)",
         bg_ratio = 5)
  } else if (n_presence < 80) {
    list(fc = c("L", "LQ"), rm = seq(2.0, 5.0, 0.5),
         description = "Linear/quadratic, high regularization (n = 50-80)",
         bg_ratio = 5)
  } else {
    list(fc = c("L", "LQ", "LQH"), rm = seq(1.0, 4.0, 0.5),
         description = "Full complexity range (n >= 80)",
         bg_ratio = 10)
  }
}

# ============================================================================
# 8. ADAPTIVE BACKGROUND GENERATION
# ============================================================================

generate_adaptive_background <- function(n_presence, predictors, bg_ratio) {
  n_background <- bg_ratio * n_presence
  cat("  Generating", n_background, "background points (",
      bg_ratio, "x ratio)\n")
  set.seed(42)
  bg_points <- dismo::randomPoints(predictors, n = n_background)
  bg_env    <- raster::extract(predictors, bg_points)
  valid     <- complete.cases(bg_env)
  if (sum(!valid) > 0) {
    cat("  Removing", sum(!valid), "background points with NAs\n")
    bg_points <- bg_points[valid, ]
  }
  bg_df <- as.data.frame(bg_points)
  colnames(bg_df) <- c("lon", "lat")
  bg_df$lon <- as.numeric(as.character(bg_df$lon))
  bg_df$lat <- as.numeric(as.character(bg_df$lat))
  bg_df <- na.omit(bg_df)
  cat("  Final background points:", nrow(bg_df), "\n")
  return(bg_df)
}

# ============================================================================
# 9. RESPONSE CURVE EXTRACTION
# ============================================================================

extract_maxnet_response <- function(model, env_data, var_name,
                                    categorical_vars = character(0),
                                    n_points = 100) {
  means <- env_data %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
  for (cv in categorical_vars) {
    if (cv %in% colnames(env_data)) {
      means[[cv]] <- as.integer(
        names(sort(table(env_data[[cv]]), decreasing = TRUE))[1])
    }
  }
  if (var_name %in% categorical_vars) {
    levels_present <- sort(unique(env_data[[var_name]]))
    newdata <- means[rep(1, length(levels_present)), ]
    newdata[[var_name]] <- levels_present
  } else {
    var_range <- seq(min(env_data[[var_name]], na.rm = TRUE),
                     max(env_data[[var_name]], na.rm = TRUE),
                     length.out = n_points)
    newdata <- means[rep(1, n_points), ]
    newdata[[var_name]] <- var_range
  }
  for (cv in categorical_vars) {
    if (cv %in% colnames(newdata)) newdata[[cv]] <- as.integer(newdata[[cv]])
  }
  preds <- predict(model, newdata, type = "logistic")
  data.frame(value = newdata[[var_name]], suitability = as.numeric(preds))
}

# ============================================================================
# 10. BOOTSTRAP FOR MAXNET
# ============================================================================

run_maxnet_bootstrap <- function(model_formula_args, p_env, a_env,
                                 predictors, n_boot = 100,
                                 categorical_vars = character(0)) {
  cat("  Running", n_boot, "bootstrap replicates...\n")
  pred_stack  <- list()
  pred_df     <- as.data.frame(raster::values(predictors))
  valid_cells <- complete.cases(pred_df)
  for (cv in categorical_vars) {
    if (cv %in% colnames(pred_df)) pred_df[[cv]] <- as.integer(pred_df[[cv]])
  }
  for (b in seq_len(n_boot)) {
    boot_idx <- sample(nrow(p_env), nrow(p_env), replace = TRUE)
    boot_p   <- p_env[boot_idx, ]
    combined <- rbind(cbind(boot_p, pa = 1), cbind(a_env, pa = 0))
    boot_model <- tryCatch({
      maxnet::maxnet(
        p                      = combined$pa,
        data                   = combined[, !colnames(combined) %in% "pa",
                                          drop = FALSE],
        f                      = model_formula_args$f,
        regmult                = model_formula_args$regmult,
        addsamplestobackground = FALSE
      )
    }, error = function(e) NULL)
    if (is.null(boot_model)) next
    preds_vec <- rep(NA_real_, nrow(pred_df))
    preds_vec[valid_cells] <- as.numeric(
      predict(boot_model, pred_df[valid_cells, ], type = "logistic"))
    pred_stack[[b]] <- preds_vec
    if (b %% 20 == 0) cat("    Completed", b, "replicates\n")
  }
  pred_matrix <- do.call(cbind, pred_stack)
  mean_vals   <- rowMeans(pred_matrix, na.rm = TRUE)
  sd_vals     <- apply(pred_matrix, 1, sd, na.rm = TRUE)
  mean_raster <- raster(predictors[[1]])
  sd_raster   <- raster(predictors[[1]])
  raster::values(mean_raster) <- mean_vals
  raster::values(sd_raster)   <- sd_vals
  list(mean = mean_raster, sd = sd_raster)
}

# ============================================================================
# 11. PERMUTATION VARIABLE IMPORTANCE
# ============================================================================

calculate_permutation_importance <- function(final_model, p_env, a_env,
                                             n_permutations = 10) {
  env_all <- rbind(p_env, a_env)
  labels  <- c(rep(1, nrow(p_env)), rep(0, nrow(a_env)))
  
  pred_full    <- as.numeric(predict(final_model, env_all, type = "logistic"))
  baseline_auc <- as.numeric(pROC::auc(pROC::roc(labels, pred_full,
                                                 quiet = TRUE)))
  
  vars       <- colnames(p_env)
  importance <- data.frame(variable = vars,
                           importance = NA_real_,
                           sd        = NA_real_)
  
  for (v in seq_along(vars)) {
    perm_aucs <- numeric(n_permutations)
    for (pp in seq_len(n_permutations)) {
      env_perm       <- env_all
      env_perm[, v]  <- sample(env_perm[, v])
      pred_perm      <- as.numeric(predict(final_model, env_perm,
                                           type = "logistic"))
      perm_aucs[pp]  <- as.numeric(pROC::auc(pROC::roc(labels, pred_perm,
                                                       quiet = TRUE)))
    }
    drop               <- baseline_auc - perm_aucs
    importance[v, "importance"] <- max(0, mean(drop))
    importance[v, "sd"]         <- sd(drop)
  }
  
  total <- sum(importance$importance, na.rm = TRUE)
  if (total > 0) importance$importance <- 100 * importance$importance / total
  
  importance[order(-importance$importance), ]
}

# ============================================================================
# 12. CORE ADAPTIVE MAXNET FUNCTION
# ============================================================================

run_adaptive_maxnet <- function(occurrence_coords,
                                predictors,
                                predictors_terra,
                                period_name,
                                output_dir       = "output",
                                categorical_vars = character(0)) {
  
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("PROCESSING:", period_name, "\n")
  cat("Sample size:", nrow(occurrence_coords), "\n")
  cat(rep("=", 70), "\n\n", sep = "")
  
  period_dir <- file.path(output_dir, period_name)
  if (dir.exists(period_dir)) {
    cat("  Cleaning existing directory...\n")
    unlink(period_dir, recursive = TRUE, force = TRUE)
    Sys.sleep(0.5)
  }
  dir.create(period_dir, showWarnings = FALSE, recursive = TRUE)
  
  ## CLEAN OCCURRENCE COORDINATES ----
  occurrence_coords <- as.data.frame(occurrence_coords)
  if ("utm_easting" %in% colnames(occurrence_coords)) {
    occurrence_coords <- occurrence_coords %>%
      rename(lon = utm_easting, lat = utm_northing)
  } else if ("easting" %in% colnames(occurrence_coords)) {
    occurrence_coords <- occurrence_coords %>%
      rename(lon = easting, lat = northing)
  } else if (!all(c("lon", "lat") %in% colnames(occurrence_coords))) {
    colnames(occurrence_coords)[1:2] <- c("lon", "lat")
  }
  occurrence_coords$lon <- as.numeric(as.character(occurrence_coords$lon))
  occurrence_coords$lat <- as.numeric(as.character(occurrence_coords$lat))
  occurrence_coords <- occurrence_coords[complete.cases(occurrence_coords), ]
  
  n_occ        <- nrow(occurrence_coords)
  # FIX: use unique cells for bootstrap decision, not raw individuals
  n_unique     <- nrow(dplyr::distinct(occurrence_coords[, c("lon", "lat")]))
  
  cat("  Cleaned occurrence points:", n_occ, "\n")
  cat("  Unique grid cells:", n_unique, "\n")
  cat("  Individuals per cell:", round(n_occ / n_unique, 1), "\n")
  
  write.csv(occurrence_coords,
            file.path(period_dir, "occurrence_points.csv"), row.names = FALSE)
  
  ## ADAPTIVE PARAMETERS ----
  tuning_params <- get_adaptive_tuning_params(n_occ)
  cat("Complexity strategy:", tuning_params$description, "\n")
  if (n_occ < 50) cat("  WARNING: Very small sample - using extreme regularization\n")
  
  ## ADAPTIVE BACKGROUND ----
  background_coords <- generate_adaptive_background(
    n_presence = n_occ, predictors = predictors,
    bg_ratio   = tuning_params$bg_ratio)
  write.csv(background_coords,
            file.path(period_dir, "background_points.csv"), row.names = FALSE)
  
  ## EXTRACT ENVIRONMENTAL VALUES ----
  p_coords <- occurrence_coords[, c("lon", "lat")]
  a_coords <- background_coords[, c("lon", "lat")]
  
  p_env <- as.data.frame(raster::extract(predictors, p_coords))
  a_env <- as.data.frame(raster::extract(predictors, a_coords))
  
  for (cv in categorical_vars) {
    if (cv %in% colnames(p_env)) {
      p_env[[cv]] <- as.integer(p_env[[cv]])
      a_env[[cv]] <- as.integer(a_env[[cv]])
    }
  }
  p_env <- p_env[complete.cases(p_env), ]
  a_env <- a_env[complete.cases(a_env), ]
  
  cat("  Presence env rows after NA removal:", nrow(p_env), "\n")
  cat("  Background env rows after NA removal:", nrow(a_env), "\n")
  
  ## STAGE 1: ENMeval TUNING ----
  cat("\nSTAGE 1: Model tuning with spatial cross-validation\n")
  
  p_coords_en <- data.frame(x = p_coords$lon, y = p_coords$lat)
  a_coords_en <- data.frame(x = a_coords$lon, y = a_coords$lat)
  
  enmeval_result <- tryCatch({
    ENMevaluate(
      occs         = p_coords_en,
      envs         = predictors_terra,
      bg           = a_coords_en,
      algorithm    = "maxnet",
      partitions   = "block",
      tune.args    = list(fc = tuning_params$fc, rm = tuning_params$rm),
      categoricals = categorical_vars,
      parallel     = FALSE
    )
  }, error = function(e) {
    cat("ERROR in ENMeval:", conditionMessage(e), "\n")
    stop(e)
  })
  
  results <- eval.results(enmeval_result)
  
  # THINNING DIAGNOSTIC — report effective sample size after cell thinning
  n_after_thin <- nrow(enmeval_result@occs)
  cat("\nThinning report:\n")
  cat("  Raw individuals:       ", n_occ, "\n")
  cat("  Unique cells (pre-CV): ", n_unique, "\n")
  cat("  After ENMeval thinning:", n_after_thin, "\n")
  cat("  Retained:              ",
      round(100 * n_after_thin / n_occ, 1), "% of individuals\n")
  
  cat("\nTuning complete. Testing", nrow(results), "model configurations\n")
  write.csv(results,
            file.path(period_dir, "enmeval_all_results.csv"), row.names = FALSE)
  
  ## STAGE 2: MODEL SELECTION ----
  cat("\nSTAGE 2: Model screening and selection\n")
  
  screened_models <- results %>% filter(or.10p.avg <= 0.50)
  cat("  Screening: Removed",
      nrow(results) - nrow(screened_models),
      "pathological models with OR10 > 0.50\n")
  
  all_models_failed <- FALSE
  
  if (nrow(screened_models) == 0) {
    all_models_failed <- TRUE
    cat("\nCRITICAL: ALL MODELS PATHOLOGICAL\n")
    cat("  Sample size (n =", n_occ, ") insufficient for reliable modeling\n")
    optimal <- data.frame(fc = "L", rm = max(tuning_params$rm),
                          delta.AICc = NA_real_, or.10p.avg = NA_real_,
                          auc.val.avg = NA_real_)
    model_quality <- "INSUFFICIENT - sample inadequate for modeling"
    
  } else {
    best_aicc       <- min(screened_models$delta.AICc, na.rm = TRUE)
    aicc_candidates <- screened_models %>%
      filter(delta.AICc <= best_aicc + 2)
    cat("  AICc selection:", nrow(aicc_candidates),
        "models within 2 delta AICc of best\n")
    
    if (nrow(aicc_candidates) > 1) {
      fc_complexity <- c("L" = 1, "LQ" = 2, "H" = 3,
                         "LQH" = 4, "LQHP" = 5, "LQHPT" = 6)
      aicc_candidates$complexity <-
        fc_complexity[as.character(aicc_candidates$fc)]
      optimal <- aicc_candidates %>%
        filter(complexity == min(complexity)) %>% slice(1)
      cat("  Multiple models within delta AICc < 2;",
          "selected simplest feature class\n")
    } else {
      optimal <- aicc_candidates[1, ]
    }
  }
  
  ## COERCE OPTIMAL VALUES TO PLAIN NUMERIC ----
  optimal_fc <- as.character(optimal$fc)
  if (length(optimal_fc) == 0 || is.na(optimal_fc)) optimal_fc <- "L"
  
  optimal_rm <- as.numeric(as.character(optimal$rm))
  if (length(optimal_rm) == 0 || is.na(optimal_rm)) optimal_rm <- 4
  
  or10_val <- NA_real_
  auc_val  <- NA_real_
  aicc_val <- NA_real_
  
  if (!all_models_failed) {
    or10_val <- as.numeric(as.character(optimal$or.10p.avg))
    auc_val  <- as.numeric(as.character(optimal$auc.val.avg))
    aicc_val <- as.numeric(as.character(optimal$delta.AICc))
    
    model_quality <- case_when(
      or10_val > 0.30  ~ "POOR - severe overfitting",
      or10_val > 0.20  ~ "MARGINAL - moderate overfitting",
      or10_val > 0.15  ~ "ACCEPTABLE - mild overfitting",
      or10_val <= 0.15 ~ "GOOD - minimal overfitting",
      TRUE             ~ "INSUFFICIENT - sample inadequate for modeling"
    )
  }
  
  cat("\nOptimal model:\n")
  cat("  Feature classes:", optimal_fc, "\n")
  cat("  Regularization:", optimal_rm, "\n")
  if (!all_models_failed) {
    cat("  delta AICc:", round(aicc_val, 3), "\n")
    cat("  Mean OR10:", round(or10_val, 3), "\n")
    cat("  Mean test AUC:", round(auc_val, 3), "\n")
  }
  cat("  Quality assessment:", model_quality, "\n")
  
  optimal_save <- optimal[, intersect(
    c("fc", "rm", "delta.AICc", "or.10p.avg", "auc.val.avg"),
    colnames(optimal))]
  write.csv(optimal_save,
            file.path(period_dir, "optimal_model_settings.csv"),
            row.names = FALSE)
  
  ## STAGE 3: BUILD MAXNET FORMULA ----
  fc_formula <- maxnet::maxnet.formula(
    p       = rep(1, nrow(p_env)),
    data    = rbind(p_env, a_env),
    classes = tolower(gsub("LQH", "lqh",
                           gsub("LQ", "lq",
                                gsub("^L$", "l", optimal_fc))))
  )
  model_formula_args <- list(f = fc_formula, regmult = optimal_rm)
  
  ## STAGE 4: FINAL MODEL FITTING ----
  cat("\nSTAGE 4: Fitting final model\n")
  
  # Bootstrap threshold: ≥45 unique cells = EPV ≥7.5 with 6 predictor variables
  # Following Vittinghoff & McCulloch (2007) minimum EPV of 5-10
  run_bootstrap <- (n_unique >= 45 &&
                      !all_models_failed &&
                      model_quality %in% c("GOOD - minimal overfitting",
                                           "ACCEPTABLE - mild overfitting"))
  if (n_unique >= 45 && n_unique < 80) {
    cat("  Note: Bootstrap run with n_unique =", n_unique,
        "(EPV =", round(n_unique / length(names(predictors)), 1),
        "with", length(names(predictors)), "predictors)\n")
  }
  
  combined_env <- rbind(cbind(p_env, pa = 1), cbind(a_env, pa = 0))
  
  final_model <- maxnet::maxnet(
    p                      = combined_env$pa,
    data                   = combined_env[, !colnames(combined_env) %in% "pa",
                                          drop = FALSE],
    f                      = fc_formula,
    regmult                = optimal_rm,
    addsamplestobackground = FALSE
  )
  
  if (run_bootstrap) {
    cat("  Running bootstrap for uncertainty estimation...\n")
    boot_results <- run_maxnet_bootstrap(
      model_formula_args = model_formula_args,
      p_env = p_env, a_env = a_env,
      predictors = predictors, n_boot = 100,
      categorical_vars = categorical_vars)
    prediction_map <- boot_results$mean
    cat("  Using bootstrap mean raster\n")
    writeRaster(boot_results$sd,
                file.path(period_dir, "prediction_uncertainty_sd.tif"),
                overwrite = TRUE)
    cv_raster <- boot_results$sd / boot_results$mean
    writeRaster(cv_raster,
                file.path(period_dir, "prediction_uncertainty_cv.tif"),
                overwrite = TRUE)
    cat("  Uncertainty maps saved\n")
    
    # CV reliability check — if mean CV > 0.5 fall back to single model prediction
    mean_cv <- cellStats(cv_raster, mean, na.rm = TRUE)
    cat(sprintf("  Mean CV across prediction area: %.3f", mean_cv))
    if (mean_cv > 0.5) {
      cat(" WARNING: mean CV > 0.5 — falling back to single model prediction\n")
      cat("  Bootstrap uncertainty maps retained but prediction map uses single model\n")
      pred_df     <- as.data.frame(raster::values(predictors))
      valid_cells <- complete.cases(pred_df)
      for (cv in categorical_vars) {
        if (cv %in% colnames(pred_df)) pred_df[[cv]] <- as.integer(pred_df[[cv]])
      }
      pred_vals <- rep(NA_real_, nrow(pred_df))
      pred_vals[valid_cells] <- as.numeric(
        predict(final_model, pred_df[valid_cells, ], type = "logistic"))
      prediction_map <- raster(predictors[[1]])
      raster::values(prediction_map) <- pred_vals
      run_bootstrap  <- FALSE   # record as not bootstrapped in summary
    } else {
      cat(" [OK]\n")
    }
    
  } else {
    cat("  Fitting single model (no bootstrap)...\n")
    pred_df     <- as.data.frame(raster::values(predictors))
    valid_cells <- complete.cases(pred_df)
    for (cv in categorical_vars) {
      if (cv %in% colnames(pred_df)) pred_df[[cv]] <- as.integer(pred_df[[cv]])
    }
    pred_vals <- rep(NA_real_, nrow(pred_df))
    pred_vals[valid_cells] <- as.numeric(
      predict(final_model, pred_df[valid_cells, ], type = "logistic"))
    prediction_map <- raster(predictors[[1]])
    raster::values(prediction_map) <- pred_vals
  }
  
  ## STAGE 5: SAVE PREDICTION MAP ----
  cat("\nSTAGE 5: Saving prediction map\n")
  writeRaster(prediction_map,
              file.path(period_dir, "prediction_continuous.tif"),
              overwrite = TRUE)
  
  ## STAGE 6: THRESHOLDS ----
  cat("\nSTAGE 6: Calculating thresholds\n")
  
  p_preds <- as.numeric(predict(final_model, p_env, type = "logistic"))
  a_preds <- as.numeric(predict(final_model, a_env, type = "logistic"))
  
  thresh_10ptp <- quantile(p_preds, 0.10, na.rm = TRUE)
  
  all_vals  <- sort(unique(c(p_preds, a_preds)))
  sens_spec <- sapply(all_vals, function(thr) {
    sens <- mean(p_preds >= thr, na.rm = TRUE)
    spec <- mean(a_preds <  thr, na.rm = TRUE)
    c(sens = sens, spec = spec,
      sum  = sens + spec, diff = abs(sens - spec))
  })
  thresh_maxsss <- all_vals[which.max(sens_spec["sum",  ])]
  thresh_etss   <- all_vals[which.min(sens_spec["diff", ])]
  
  thresholds_df <- data.frame(
    threshold_name  = c("10PTP", "ETSS", "MaxSSS"),
    threshold_value = c(thresh_10ptp, thresh_etss, thresh_maxsss)
  )
  write.csv(thresholds_df,
            file.path(period_dir, "thresholds.csv"), row.names = FALSE)
  
  cell_size_km2 <- prod(res(prediction_map)) / 1e6
  for (i in seq_len(nrow(thresholds_df))) {
    binary_map <- prediction_map >= thresholds_df$threshold_value[i]
    area_km2   <- cellStats(binary_map, "sum") * cell_size_km2
    cat(sprintf("  Area (%s): %.1f km2\n",
                thresholds_df$threshold_name[i], area_km2))
    writeRaster(binary_map,
                file.path(period_dir,
                          paste0("binary_", thresholds_df$threshold_name[i],
                                 ".tif")),
                overwrite = TRUE)
  }
  
  ## STAGE 7: RESPONSE CURVES ----
  cat("\nSTAGE 7: Extracting response curves\n")
  all_env      <- rbind(p_env, a_env)
  response_dir <- file.path(period_dir, "response_curves")
  dir.create(response_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (var_name in names(predictors)) {
    resp_df <- tryCatch({
      extract_maxnet_response(
        model = final_model, env_data = all_env,
        var_name = var_name, categorical_vars = categorical_vars)
    }, error = function(e) {
      cat("  Could not extract response for", var_name, ":",
          conditionMessage(e), "\n")
      NULL
    })
    if (!is.null(resp_df)) {
      write.csv(resp_df,
                file.path(response_dir,
                          paste0("response_", var_name, ".csv")),
                row.names = FALSE)
    }
  }
  
  ## STAGE 8: PERMUTATION VARIABLE IMPORTANCE ----
  cat("\nSTAGE 8: Calculating permutation variable importance\n")
  
  vi <- tryCatch({
    calculate_permutation_importance(
      final_model   = final_model,
      p_env         = p_env,
      a_env         = a_env,
      n_permutations = 10
    )
  }, error = function(e) {
    cat("  Could not calculate variable importance:", conditionMessage(e), "\n")
    NULL
  })
  
  if (!is.null(vi)) {
    write.csv(vi,
              file.path(period_dir, "variable_importance.csv"),
              row.names = FALSE)
    cat("  Variable importance saved\n")
    print(vi)
  }
  
  ## RETURN RESULTS ----
  performance_summary <- data.frame(
    period         = period_name,
    n_occ          = n_occ,
    n_unique_cells = n_unique,
    n_after_thin   = n_after_thin,
    pct_retained   = round(100 * n_after_thin / n_occ, 1),
    n_background   = nrow(background_coords),
    complexity     = tuning_params$description,
    fc             = optimal_fc,
    rm             = optimal_rm,
    quality        = model_quality,
    bootstrapped   = run_bootstrap,
    mean_test_auc  = auc_val,
    mean_or10      = or10_val,
    mean_cv        = ifelse(run_bootstrap, mean_cv, NA_real_)
  )
  
  return(list(
    period              = period_name,
    n_occurrence        = n_occ,
    n_unique_cells      = n_unique,
    n_after_thin        = n_after_thin,
    n_background        = nrow(background_coords),
    performance_summary = performance_summary,
    prediction_map      = prediction_map,
    model_quality       = model_quality,
    all_models_failed   = all_models_failed,
    bootstrapped        = run_bootstrap,
    variable_importance = vi
  ))
}

# ============================================================================
# 13. RUN ALL PERIODS
# ============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("RUNNING ADAPTIVE MAXNET FOR ALL PERIODS\n")
cat(rep("=", 70), "\n\n", sep = "")

all_results    <- list()
failed_periods <- character()

for (period_name in names(period_data)) {
  result <- tryCatch({
    run_adaptive_maxnet(
      occurrence_coords = period_data[[period_name]],
      predictors        = predictors,
      predictors_terra  = predictors_terra,
      period_name       = period_name,
      output_dir        = "output",
      categorical_vars  = categorical_vars
    )
  }, error = function(e) {
    cat("\nERROR in", period_name, ":", conditionMessage(e), "\n\n")
    failed_periods <<- c(failed_periods, period_name)
    return(NULL)
  })
  if (!is.null(result)) {
    all_results[[period_name]] <- result
    saveRDS(result, file = file.path("output", period_name, "results.rds"))
  }
}

saveRDS(all_results, "output/all_results.rds")

if (length(failed_periods) > 0) {
  cat("\nFailed periods:", paste(failed_periods, collapse = ", "), "\n")
} else {
  cat("\nAll periods completed successfully!\n")
}

# ============================================================================
# 14. COMPILE RESULTS SUMMARY
# ============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("COMPILING RESULTS SUMMARY\n")
cat(rep("=", 70), "\n\n", sep = "")

performance_list <- list()
for (period_name in names(all_results)) {
  res <- all_results[[period_name]]
  if (!is.null(res$performance_summary)) {
    performance_list[[period_name]] <- res$performance_summary
  }
}

if (length(performance_list) > 0) {
  performance_table <- do.call(rbind, performance_list)
  write.csv(performance_table,
            "output/performance_summary_all_periods.csv", row.names = FALSE)
  
  cat("\nPerformance Summary:\n")
  print(performance_table)
  
  # Thinning summary table — important for methods section reporting
  thinning_summary <- performance_table %>%
    dplyr::select(period, n_occ, n_unique_cells, n_after_thin, pct_retained)
  cat("\nThinning Summary (for methods reporting):\n")
  print(thinning_summary)
  write.csv(thinning_summary, "output/thinning_summary.csv", row.names = FALSE)
  
  quality_summary <- performance_table %>%
    group_by(quality) %>%
    summarise(
      n_periods      = dplyr::n(),
      mean_n_occ     = mean(n_occ,        na.rm = TRUE),
      mean_auc       = mean(mean_test_auc, na.rm = TRUE),
      mean_or10      = mean(mean_or10,     na.rm = TRUE),
      n_bootstrapped = sum(bootstrapped,   na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("\nQuality Assessment Summary:\n")
  print(quality_summary)
  write.csv(quality_summary, "output/quality_summary.csv", row.names = FALSE)
  
  # Variable importance summary across all periods
  vi_list <- lapply(names(all_results), function(pn) {
    vi <- all_results[[pn]]$variable_importance
    if (!is.null(vi)) { vi$period <- pn; vi }
  })
  vi_all <- do.call(rbind, Filter(Negate(is.null), vi_list))
  if (!is.null(vi_all)) {
    write.csv(vi_all, "output/variable_importance_all_periods.csv",
              row.names = FALSE)
    cat("\nVariable importance saved across all periods\n")
  }
}

cat("\n", rep("=", 70), "\n", sep = "")
cat("ADAPTIVE MAXNET WORKFLOW COMPLETE\n")
cat(rep("=", 70), "\n\n", sep = "")
cat("Results saved in output/ directory\n")
cat("Key outputs for methods reporting:\n")
cat("  thinning_summary.csv               — raw vs thinned sample sizes\n")
cat("  performance_summary_all_periods.csv — full model performance\n")
cat("  variable_importance_all_periods.csv — permutation importance\n")
