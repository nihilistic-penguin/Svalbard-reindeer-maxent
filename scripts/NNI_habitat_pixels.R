# ============================================================================
# NNI AND AGGREGATION INDEX ON SUITABLE HABITAT PIXELS
# Replaces previous occurrence-based NNI analysis
# NNI treated as descriptive index only (not formal test) due to raster lattice
# constraints. AI added as complementary contiguity metric.
# Follows Beumer et al. (2019) and Pedersen et al. (2023)
# ============================================================================

library(raster)
library(terra)
library(spatstat.geom)
library(spatstat.explore)
library(landscapemetrics)
library(dplyr)
library(readr)
library(ggplot2)
library(rcartocolor)
library(colorspace)
library(patchwork)

setwd("C:/Users/marie/Documents/R/Barrett's R/Svalbard Reindeer Maxent Paper")

# ============================================================================
# PERIOD DEFINITIONS — 8 MARGINAL-TO-GOOD PERIODS
# ============================================================================

period_order <- c(
  "weeks_09_10", "weeks_15b_16",
  "weeks_18a_18b", "weeks_20a_20b", "weeks_21_22",
  "weeks_23_24", "weeks_25a_25b", "weeks_27_28"
)

period_labels <- c(
  weeks_09_10   = "early\nMarch",
  weeks_15b_16  = "late\nApril",
  weeks_18a_18b = "early\nMay",
  weeks_20a_20b = "mid\nMay",
  weeks_21_22   = "late\nMay",
  weeks_23_24   = "early\nJune",
  weeks_25a_25b = "late\nJune",
  weeks_27_28   = "mid\nJuly"
)

cell_size    <- 20          # raster resolution in metres
cell_area_m2 <- cell_size^2 # 400 m²

# ============================================================================
# LOAD STUDY AREA BOUNDARY FROM MASKED RASTER
# Used to define the observation window for NNI expected distance calculation
# ============================================================================

cat("Loading study area boundary...\n")
pred_masked  <- raster("data/prepared for analyses/predictors_masked2.grd")
study_area_m2 <- sum(!is.na(getValues(pred_masked))) * cell_area_m2
cat("  Study area:", round(study_area_m2 / 1e6, 3), "km2\n")

# ============================================================================
# CORE FUNCTION: NNI + AI FOR ONE BINARY RASTER
# NNI calculated using Clark-Evans formula with study-area-corrected expected
# distance: E[d] = 0.5 / sqrt(n / A) where A = study area in same units as coords
# ============================================================================

calculate_habitat_metrics <- function(binary_raster, study_area_m2, cell_size = 20) {

  # Coordinates of suitable cell centroids
  suitable_cells  <- Which(binary_raster == 1, cells = TRUE)
  n_suitable      <- length(suitable_cells)

  if (n_suitable < 2) {
    return(list(nni = NA_real_, ai = NA_real_,
                n_suitable = n_suitable,
                area_km2 = n_suitable * cell_area_m2 / 1e6))
  }

  suitable_coords <- xyFromCell(binary_raster, suitable_cells)

  # NNI: Clark-Evans formula with study-area correction
  # Observed mean NN distance
  nn_dists    <- spatstat.geom::nndist(suitable_coords)
  obs_mean_nn <- mean(nn_dists)

  # Expected mean NN distance under CSR given study area and n
  intensity   <- n_suitable / study_area_m2
  exp_mean_nn <- 0.5 / sqrt(intensity)

  nni_val <- obs_mean_nn / exp_mean_nn

  # Aggregation Index via landscapemetrics
  # Requires binary 0/1 raster with no NAs — set NA to 0 (background)
  binary_copy <- binary_raster
  binary_copy[is.na(binary_copy)] <- 0

  binary_terra <- terra::rast(binary_copy)

  ai_result <- tryCatch({
    ai_df <- landscapemetrics::lsm_c_ai(binary_terra)
    # Filter for class 1 (suitable habitat)
    ai_df$value[ai_df$class == 1]
  }, error = function(e) NA_real_)

  ai_val <- if (length(ai_result) == 0) NA_real_ else ai_result

  list(
    nni        = nni_val,
    obs_mean_nn = obs_mean_nn,
    exp_mean_nn = exp_mean_nn,
    ai         = ai_val,
    n_suitable = n_suitable,
    area_km2   = n_suitable * cell_area_m2 / 1e6
  )
}

# ============================================================================
# CALCULATE METRICS ACROSS ALL PERIODS AND THRESHOLDS
# ============================================================================

cat("\nCalculating NNI and AI across periods and thresholds...\n")

results_list <- list()

for (pn in period_order) {

  cat("  Period:", pn, "\n")
  pdir <- file.path("output", pn)

  # Load continuous suitability raster
  cont_file <- file.path(pdir, "prediction_continuous.tif")
  if (!file.exists(cont_file)) { cat("    Missing raster\n"); next }
  suit_r <- raster(cont_file)

  # Load period-specific thresholds
  thresh_file <- file.path(pdir, "thresholds.csv")
  if (!file.exists(thresh_file)) { cat("    Missing thresholds\n"); next }
  thresholds <- read_csv(thresh_file, show_col_types = FALSE)

  for (i in seq_len(nrow(thresholds))) {
    tname <- thresholds$threshold_name[i]
    tval  <- thresholds$threshold_value[i]

    binary_r <- suit_r >= tval

    metrics <- calculate_habitat_metrics(binary_r, study_area_m2, cell_size)

    results_list[[paste(pn, tname)]] <- data.frame(
      period        = pn,
      period_label  = factor(period_labels[pn],
                             levels = period_labels[period_order]),
      threshold     = tname,
      threshold_val = tval,
      nni           = metrics$nni,
      obs_mean_nn   = metrics$obs_mean_nn,
      exp_mean_nn   = metrics$exp_mean_nn,
      ai            = metrics$ai,
      n_suitable    = metrics$n_suitable,
      area_km2      = metrics$area_km2
    )

    cat(sprintf("    %s: NNI=%.3f  AI=%.1f  area=%.2f km2\n",
                tname, metrics$nni, metrics$ai, metrics$area_km2))
  }
}

nni_df <- bind_rows(results_list)

write.csv(nni_df, "output/habitat_nni_ai_all_periods.csv", row.names = FALSE)
cat("\nResults saved to output/habitat_nni_ai_all_periods.csv\n")

# ============================================================================
# SUPPLEMENTARY TABLE
# ============================================================================

supp_table <- nni_df %>%
  mutate(
    period_name = gsub("\n", " ", as.character(period_label)),  # <- swap newline for space
    nni_r       = round(nni, 3),
    ai_r        = round(ai, 1),
    area_r      = round(area_km2, 2),
    obs_nn_r    = round(obs_mean_nn, 1),
    exp_nn_r    = round(exp_mean_nn, 1)
  ) %>%
  dplyr::select(
    "Period"                      = period_name,
    "Threshold"                   = threshold,
    "Threshold value"             = threshold_val,
    "Suitable cells (n)"          = n_suitable,
    "Area (km2)"                  = area_r,
    "Observed mean NN dist. (m)"  = obs_nn_r,
    "Expected mean NN dist. (m)"  = exp_nn_r,
    "NNI"                         = nni_r,
    "AI (%)"                      = ai_r
  )

cat("\nSupplementary Table:\n")
print(supp_table)

write.csv(supp_table, "output/supplementary_habitat_nni_table.csv",
          row.names = FALSE)

# ============================================================================
# PANEL B: NNI PLOT — coloured by NNI value, faceted by threshold
# ============================================================================

cat("\nBuilding NNI figure...\n")

# Threshold display order
nni_df$threshold <- factor(nni_df$threshold,
                           levels = c("10PTP", "MaxSSS", "ETSS"))

shared_theme <- theme_minimal(base_size = 10) +
  theme(
    axis.text.x      = element_text(size = 8),
    axis.text.y      = element_text(size = 8),
    axis.title       = element_text(size = 11),
    panel.grid.major = element_line(colour = "grey75", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_blank(),
    strip.text       = element_text(face = "bold", size = 9.5),
    strip.background = element_rect(fill = "grey82", colour = NA),
    plot.margin      = margin(10, 10, 10, 10)
  )

thresh_cols <- setNames(
  carto_pal(3, "Teal"),
  c("10PTP", "MaxSSS", "ETSS")
)


p_nni <- ggplot(nni_df, aes(x = period_label, y = nni,
                            fill = threshold,
                            group = threshold)) +
  geom_hline(yintercept = 1, linetype = "dashed",
             colour = "grey40", linewidth = 0.6) +
  geom_col(position = position_dodge(width = 0.75),
           width = 0.7, colour = "grey30", linewidth = 0.2) +
  geom_text(aes(label = round(nni, 2),
                y     = nni + 0.03),
            position = position_dodge(width = 0.75),
            size     = 2.2,
            fontface = "bold",
            colour   = "grey20") +
  scale_fill_manual(
    values = thresh_cols,
    name   = "Threshold"
  ) +
  scale_y_continuous(limits = c(0, 1.6),
                     breaks = seq(0, 1.5, 0.25)) +
  annotate("text", x = 0.5, y = 1.02,
           label = "random", hjust = 0, size = 3,
           colour = "grey40", fontface = "italic") +
  labs(x = "Period", y = "Nearest Neighbour Index (NNI)") +
  shared_theme +
  theme(
    legend.position  = "right",
    legend.title     = element_text(size = 9),
    legend.text      = element_text(size = 8),
    legend.key.width = unit(0.5, "cm")
  )

print(p_nni)
# ============================================================================
# SAVE — uncomment once happy
# ============================================================================

# figures_dir <- "figures"
# dir.create(figures_dir, showWarnings = FALSE)
#
# ggsave(file.path(figures_dir, "habitat_nni_ai.png"),
#        p_nni, width = 12, height = 5, dpi = 300, units = "in")
# ggsave(file.path(figures_dir, "habitat_nni_ai.tiff"),
#        p_nni, width = 12, height = 5, dpi = 300, units = "in",
#        compression = "lzw")