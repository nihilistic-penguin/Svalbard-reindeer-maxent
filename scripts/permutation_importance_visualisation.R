# ============================================================================
# PERMUTATION VARIABLE IMPORTANCE — BAR PLOTS
# Main figure: 8 periods MARGINAL-GOOD (4x2)
# Appendix figure: all 9 periods (3x3)
# One panel per period, bars per variable, error bars from SD
# ============================================================================

library(ggplot2)
library(dplyr)
library(readr)
library(rcartocolor)
library(colorspace)
library(patchwork)

setwd("C:/Users/marie/Documents/R/Barrett's R/Svalbard Reindeer Maxent Paper")

# ============================================================================
# PERIOD DEFINITIONS
# ============================================================================

# All 9 periods (appendix)
all_period_order <- c(
  "weeks_09_10", "weeks_13_15a", "weeks_15b_16",
  "weeks_18a_18b", "weeks_20a_20b", "weeks_21_22",
  "weeks_23_24", "weeks_25a_25b", "weeks_27_28"
)

all_period_labels <- c(
  weeks_09_10   = "early March",
  weeks_13_15a  = "early April",
  weeks_15b_16  = "late April",
  weeks_18a_18b = "early May",
  weeks_20a_20b = "mid May",
  weeks_21_22   = "late May",
  weeks_23_24   = "early June",
  weeks_25a_25b = "late June",
  weeks_27_28   = "mid July"
)

# 8 periods for main figure (MARGINAL to GOOD, excludes POOR and INSUFFICIENT)
main_period_order <- c(
  "weeks_09_10", "weeks_15b_16",
  "weeks_18a_18b", "weeks_20a_20b", "weeks_21_22",
  "weeks_23_24", "weeks_25a_25b", "weeks_27_28"
)

main_period_labels <- all_period_labels[main_period_order]

# Variable labels
var_labels <- c(
  elevation    = "Elevation",
  hillshade    = "Hillshade",
  ruggedness   = "Ruggedness",
  landforms_fc = "Landform class",
  landcover_fc = "Land cover class",
  NDVI22.23    = "NDVI"
)

# Poor fit periods (appendix only)
poor_periods <- c("weeks_13_15a")

# ============================================================================
# LOAD QUALITY LABELS
# ============================================================================

perf <- read.csv("output/performance_summary_all_periods.csv")

quality_short <- c(
  "GOOD - minimal overfitting"                    = "GOOD",
  "ACCEPTABLE - mild overfitting"                 = "ACCEPTABLE",
  "MARGINAL - moderate overfitting"               = "MARGINAL",
  "POOR - severe overfitting"                     = "POOR",
  "INSUFFICIENT - sample inadequate for modeling" = "INSUFFICIENT"
)

get_quality_tag <- function(pn, labels) {
  labels[pn]
}

# ============================================================================
# LOAD VARIABLE IMPORTANCE
# ============================================================================

cat("Loading variable importance...\n")

vi_list <- list()

for (pn in all_period_order) {
  f <- file.path("output", pn, "variable_importance.csv")
  if (!file.exists(f)) { cat("  Missing:", pn, "\n"); next }
  
  df <- read_csv(f, show_col_types = FALSE)
  df$period       <- pn
  df$period_label <- get_quality_tag(pn, all_period_labels)
  df$var_label    <- factor(var_labels[df$variable],
                            levels = var_labels)
  df$is_poor      <- pn %in% poor_periods
  vi_list[[pn]]   <- df
  cat("  Loaded:", pn, "\n")
}

vi_all  <- bind_rows(vi_list)
vi_main <- vi_all %>% filter(period %in% main_period_order)

cat("  Loaded", nrow(vi_all), "rows\n")

# ============================================================================
# COLOUR PALETTE — darker Earth for bars
# ============================================================================

earth_cols  <- rev(carto_pal(7, "BrwnYl"))
darker_cols <- darken(earth_cols, amount = 0.15)
bar_cols    <- setNames(
  colorRampPalette(darker_cols)(length(var_labels)),
  var_labels
)

# ============================================================================
# SHARED THEME
# ============================================================================

shared_theme <- theme_grey(base_size = 10) +
  theme(
    strip.text        = element_text(face = "bold", size = 10),
    strip.background  = element_rect(fill = "grey82", colour = NA),
    axis.text.x       = element_text(size = 11, angle = 30, hjust = 1),
    axis.text.y       = element_text(size = 11),
    axis.title        = element_text(size = 14),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    panel.grid.major  = element_line(colour = "grey70", linewidth = 0.3),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(colour = "grey60"),
    legend.position   = "none",
    plot.margin       = margin(10, 10, 10, 10), 
    panel.background = element_rect(fill = "white")
  )

# ============================================================================
# HELPER: build one bar plot figure
# ============================================================================

build_importance_fig <- function(vi_df, period_ord, period_lbls, ncols) {
  
  vi_df$period_label <- factor(
    vi_df$period_label,
    levels = sapply(period_ord, get_quality_tag, labels = period_lbls)
  )
  
  poor_labs <- vi_df %>%
    filter(is_poor) %>%
    group_by(period_label) %>%
    summarise(x = var_labels[length(var_labels)],
              y = max(importance + sd, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(label = "*") %>%
    filter(is.finite(y))
  
  ggplot(vi_df, aes(x = var_label, y = importance, fill = var_label)) +
    geom_col(width = 0.7) +
    geom_errorbar(aes(ymin = pmax(importance - sd, 0),
                      ymax = importance + sd),
                  width     = 0.25,
                  linewidth = 0.5,
                  colour    = "grey30") +
    geom_text(data        = poor_labs,
              aes(x = x, y = y, label = label),
              inherit.aes = FALSE,
              colour      = "darkred",
              size        = 5,
              fontface    = "bold",
              hjust       = 1,
              vjust       = 0) +
    scale_fill_manual(values = bar_cols) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    facet_wrap(~ period_label, ncol = ncols) +
    labs(x = NULL, y = "Permutation importance (%)") +
    shared_theme
}

# ============================================================================
# BUILD FIGURES
# ============================================================================

cat("\nBuilding main importance figure (8 periods)...\n")
p_imp_main <- build_importance_fig(vi_main, main_period_order,
                                   main_period_labels, ncols = 4)
print(p_imp_main)

cat("\nBuilding appendix importance figure (9 periods)...\n")
p_imp_all  <- build_importance_fig(vi_all, all_period_order,
                                   all_period_labels, ncols = 3)
print(p_imp_all)

# ============================================================================
# SAVE — uncomment once happy
# ============================================================================

figures_dir <- "figures"
dir.create(figures_dir, showWarnings = FALSE)

ggsave(file.path(figures_dir, "FINAL_MAIN_permutation_importance.png"),
       p_imp_main, width = 10, height = 6, dpi = 600, units = "in")
ggsave(file.path(figures_dir, "FINAL_MAIN_permutation_importance.tiff"),
       p_imp_main, width = 10, height = 6, dpi = 600, units = "in",
       compression = "lzw")

ggsave(file.path(figures_dir, "FINAL_APPENDIX_permutation_importance_all.png"),
       p_imp_all, width = 9, height = 7, dpi = 600, units = "in")
ggsave(file.path(figures_dir, "FINAL_APPENDIX_permutation_importance_all.tiff"),
       p_imp_all, width = 9, height = 7, dpi = 600, units = "in",
       compression = "lzw")
