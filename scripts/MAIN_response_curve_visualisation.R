# ============================================================================
# RESPONSE CURVES — 8 PERIODS (MARGINAL TO GOOD), MAIN FIGURE
# 6 panels (one per variable), 8 lines per panel coloured by Geyser palette
# Excludes weeks_11_12 (INSUFFICIENT) and weeks_13_15a (POOR)
# ============================================================================

library(ggplot2)
library(dplyr)
library(readr)
library(rcartocolor)
library(patchwork)

setwd("C:/Users/marie/Documents/R/Barrett's R/Svalbard Reindeer Maxent Paper")

# ============================================================================
# SELECTED PERIODS AND LABELS
# ============================================================================

period_order <- c(
  "weeks_09_10",    # early March  — MARGINAL
  "weeks_15b_16",   # late April   — GOOD
  "weeks_18a_18b",  # early May    — ACCEPTABLE
  "weeks_20a_20b",  # mid May      — GOOD
  "weeks_21_22",    # late May     — ACCEPTABLE
  "weeks_23_24",    # early June   — GOOD
  "weeks_25a_25b",  # late June    — GOOD (best)
  "weeks_27_28"     # mid July     — MARGINAL
)

period_labels <- c(
  weeks_09_10   = "early March",
  weeks_15b_16  = "late April",
  weeks_18a_18b = "early May",
  weeks_20a_20b = "mid May",
  weeks_21_22   = "late May",
  weeks_23_24   = "early June",
  weeks_25a_25b = "late June",
  weeks_27_28   = "mid July"
)

# Variable labels
var_labels <- c(
  elevation    = "Elevation (m)",
  hillshade    = "Hillshade",
  ruggedness   = "Ruggedness",
  landforms_fc = "Landform class",
  landcover_fc = "Land cover class",
  NDVI22.23    = "NDVI"
)

continuous_vars  <- c("elevation", "hillshade", "ruggedness", "NDVI22.23")
categorical_vars <- c("landforms_fc", "landcover_fc")
var_order        <- c(continuous_vars, categorical_vars)

# ============================================================================
# LOAD QUALITY LABELS FROM PERFORMANCE SUMMARY
# ============================================================================

perf <- read.csv("output/performance_summary_all_periods.csv")

quality_short <- c(
  "GOOD - minimal overfitting"                    = "GOOD",
  "ACCEPTABLE - mild overfitting"                 = "ACCEPTABLE",
  "MARGINAL - moderate overfitting"               = "MARGINAL",
  "POOR - severe overfitting"                     = "POOR",
  "INSUFFICIENT - sample inadequate for modeling" = "INSUFFICIENT"
)

# Build period labels without quality tags
quality_tags <- period_labels[period_order]
names(quality_tags) <- period_labels[period_order]

# ============================================================================
# COLOUR PALETTE — Geyser, one colour per period
# ============================================================================

geyser_8       <- carto_pal(8, "Geyser")
season_cols_tagged <- setNames(geyser_8,
                               quality_tags[period_labels[period_order]])

# All solid lines
ltype_map <- setNames(rep("solid", length(period_order)),
                      quality_tags[period_labels[period_order]])

# ============================================================================
# SHARED THEME
# ============================================================================

shared_theme <- theme_bw(base_size = 10) +
  theme(
    strip.text        = element_text(face = "bold", size = 9.5),
    strip.background  = element_rect(fill = "grey82", colour = NA),
    axis.text.x       = element_text(size = 11),
    axis.text.y       = element_text(size = 11),
    axis.title        = element_text(size = 14, face = "bold"),
    panel.grid.major  = element_line(colour = "grey88", linewidth = 0.3),
    panel.grid.minor  = element_blank(),
    panel.border      = element_blank(),
    legend.position   = "bottom",
    legend.title      = element_text(size = 12),
    legend.text       = element_text(size = 11),
    legend.key.width  = unit(0.7, "cm"),
    plot.margin       = margin(10, 10, 5, 10)
  )

# ============================================================================
# LOAD RESPONSE CURVES
# ============================================================================

cat("Loading response curves...\n")

rc_list <- list()

for (pn in period_order) {
  rc_dir <- file.path("output", pn, "response_curves")
  if (!dir.exists(rc_dir)) { cat("  Missing:", pn, "\n"); next }
  
  for (vn in var_order) {
    f <- file.path(rc_dir, paste0("response_", vn, ".csv"))
    if (!file.exists(f)) next
    
    df <- read_csv(f, show_col_types = FALSE)
    df$period       <- pn
    df$period_label <- factor(quality_tags[period_labels[pn]],
                              levels = quality_tags[period_labels[period_order]])
    df$variable     <- vn
    df$var_label    <- factor(var_labels[vn],
                              levels = var_labels[var_order])
    rc_list[[paste(pn, vn)]] <- df
  }
}

rc_all  <- bind_rows(rc_list)
rc_cont <- rc_all %>% filter(variable %in% continuous_vars)
rc_cat  <- rc_all %>% filter(variable %in% categorical_vars)

cat("  Loaded", nrow(rc_all), "rows across",
    length(unique(rc_all$period)), "periods\n")

# ============================================================================
# CONTINUOUS PREDICTORS — lines
# ============================================================================

cat("\nBuilding continuous response curves...\n")

p_cont <- ggplot(rc_cont,
                 aes(x        = value,
                     y        = suitability,
                     colour   = period_label,
                     linetype = period_label,
                     group    = period_label)) +
  geom_line(linewidth = 0.9, alpha = 0.9) +
  scale_colour_manual(values = season_cols_tagged,
                      name   = "Period (model quality)") +
  scale_linetype_manual(values = ltype_map,
                        name   = "Period (model quality)") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  facet_wrap(~ var_label, scales = "free_x", ncol = 4) +
  labs(x = "Variable value", y = "Predicted suitability") +
  shared_theme +
  theme(legend.position = "none", 
        axis.title = element_text(size = 14))

# ============================================================================
# CATEGORICAL PREDICTORS — points
# ============================================================================

cat("Building categorical response curves...\n")

p_cat <- ggplot(rc_cat,
                aes(x      = value,
                    y      = suitability,
                    colour = period_label,
                    shape  = period_label,
                    group  = period_label)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_line(linewidth = 0.6, alpha = 0.9, linetype = "dotted") +
  scale_colour_manual(values = season_cols_tagged,
                      name   = "Period ") +
  scale_shape_manual(values = c(16, 17, 15, 18, 16, 17, 15, 18),
                     name   = "Period ") +
  scale_y_continuous(limits = c(0, 0.5), breaks = c(0, 0.25, 0.5)#, 0.75, 1)
                     ) +
  facet_wrap(~ var_label, scales = "free_x", ncol = 2) +
  labs(x = "Class", y = "Predicted suitability") +
  shared_theme

# ============================================================================
# COMBINE
# ============================================================================

p_rc_main <- p_cont / p_cat +
  plot_layout(heights = c(1, 1)) &
  theme(plot.margin = margin(5, 10, 5, 10))

print(p_rc_main)

# ============================================================================
# SAVE — uncomment once happy
# ============================================================================

figures_dir <- "figures"
dir.create(figures_dir, showWarnings = FALSE)

ggsave(file.path(figures_dir, "MAIN_response_curves.png"),
       p_rc_main, width = 9, height = 6, dpi = 600, units = "in")
ggsave(file.path(figures_dir, "MAIN_response_curves.tiff"),
       p_rc_main, width = 9, height = 6, dpi = 600, units = "in",
       compression = "lzw")
