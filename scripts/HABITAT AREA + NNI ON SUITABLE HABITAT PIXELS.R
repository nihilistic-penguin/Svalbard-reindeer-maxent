# ============================================================================
# DUAL PANEL FIGURE — HABITAT AREA + NNI ON SUITABLE HABITAT PIXELS
# Panel A: Suitable habitat area (km²) by threshold and period
# Panel B: NNI on suitable habitat pixels by threshold and period
# X-axis: all 10 period labels including excluded periods (shown as gaps)
# White background, Teal palette
# ============================================================================

library(ggplot2)
library(dplyr)
library(readr)
library(raster)
library(patchwork)
library(rcartocolor)

setwd("C:/Users/marie/Documents/R/Barrett's R/Svalbard Reindeer Maxent Paper")

# ============================================================================
# 1. PERIOD DEFINITIONS
# ============================================================================

# All 10 periods — including excluded ones for x-axis spacing
all_period_order <- c(
  "weeks_09_10", "weeks_11_12", "weeks_13_15a", "weeks_15b_16",
  "weeks_18a_18b", "weeks_20a_20b", "weeks_21_22",
  "weeks_23_24", "weeks_25a_25b", "weeks_27_28"
)

all_period_labels <- c(
  weeks_09_10   = "early March",
  weeks_11_12   = "late March",
  weeks_13_15a  = "early April",
  weeks_15b_16  = "late April",
  weeks_18a_18b = "early May",
  weeks_20a_20b = "mid May",
  weeks_21_22   = "late May",
  weeks_23_24   = "early June",
  weeks_25a_25b = "late June",
  weeks_27_28   = "mid July"
)

# 8 modelled periods
period_order <- c(
  "weeks_09_10", "weeks_15b_16", "weeks_18a_18b", "weeks_20a_20b",
  "weeks_21_22", "weeks_23_24", "weeks_25a_25b", "weeks_27_28"
)

bookend_periods <- c("weeks_09_10", "weeks_27_28")

# Factor with all 10 levels for x-axis — gaps appear for missing data
x_levels <- all_period_labels[all_period_order]

# ============================================================================
# 2. THRESHOLD SETTINGS — Teal palette
# ============================================================================

teal_cols <- carto_pal(3, "Teal")

thresh_settings <- data.frame(
  threshold = c("10PTP", "MaxSSS", "ETSS"),
  label     = c("10PTP (core)", "MaxSSS (balanced)", "ETSS (equal sens/spec)"),
  colour    = teal_cols,
  stringsAsFactors = FALSE
)

thresh_cols <- setNames(teal_cols, c("10PTP", "MaxSSS", "ETSS"))

# ============================================================================
# 3. SHARED ELEMENTS
# ============================================================================

shared_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.text.x      = element_text(size = 10),
    axis.text.y      = element_text(size = 10),
    axis.title.x = element_text(size = 13, face = "bold"),
    axis.title.y = element_text(size = 13, face = "bold"),
    panel.grid.major = element_line(colour = "grey88", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_blank(),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin      = margin(10, 0, 5, 15), 
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

period_xscale <- scale_x_discrete(
  limits = x_levels,
  labels = function(x) gsub(" ", "\n", x),
  drop   = FALSE
)

# ============================================================================
# 4. CALCULATE HABITAT AREA
# ============================================================================

cat("Calculating suitable habitat areas...\n")

pred_masked   <- raster("data/prepared for analyses/predictors_masked2.grd")
cell_size_km2 <- prod(res(pred_masked)) / 1e6

area_list <- list()

for (pn in period_order) {
  pdir        <- file.path("output", pn)
  thresh_file <- file.path(pdir, "thresholds.csv")
  cont_file   <- file.path(pdir, "prediction_continuous.tif")
  if (!file.exists(thresh_file) | !file.exists(cont_file)) next
  
  thresholds <- read_csv(thresh_file, show_col_types = FALSE)
  pred_r     <- raster(cont_file)
  sd_file    <- file.path(pdir, "prediction_uncertainty_sd.tif")
  has_sd     <- file.exists(sd_file)
  if (has_sd) sd_r <- raster(sd_file)
  
  for (i in seq_len(nrow(thresholds))) {
    thr   <- thresholds$threshold_value[i]
    tname <- thresholds$threshold_name[i]
    
    binary   <- pred_r >= thr
    area_km2 <- cellStats(binary, "sum") * cell_size_km2
    
    area_sd <- NA_real_
    if (has_sd) {
      area_upper <- cellStats((pred_r + sd_r) >= thr, "sum") * cell_size_km2
      area_lower <- cellStats((pred_r - sd_r) >= thr, "sum") * cell_size_km2
      area_sd    <- (area_upper - area_lower) / 2
    }
    
    area_list[[paste(pn, tname)]] <- data.frame(
      period       = pn,
      x_label      = factor(all_period_labels[pn], levels = x_levels),
      threshold    = tname,
      area_km2     = area_km2,
      area_sd      = area_sd,
      bootstrapped = has_sd,
      is_bookend   = pn %in% bookend_periods
    )
  }
  cat("  Done:", pn, "\n")
}

area_df <- bind_rows(area_list) %>%
  left_join(thresh_settings[, c("threshold", "label")], by = "threshold") %>%
  mutate(threshold_label = factor(label, levels = thresh_settings$label))

# ============================================================================
# 5. LOAD NNI DATA
# ============================================================================

cat("\nLoading habitat NNI results...\n")

nni_df <- read.csv("output/habitat_nni_ai_all_periods.csv") %>%
  mutate(
    x_label   = factor(all_period_labels[period], levels = x_levels),
    threshold = factor(threshold, levels = c("10PTP", "MaxSSS", "ETSS"))
  )

cat("  Loaded", nrow(nni_df), "rows\n")

# ============================================================================
# 6. BUILD PANEL A — HABITAT AREA
# ============================================================================

cat("\nBuilding Panel A...\n")

p_area <- ggplot(area_df,
                 aes(x = x_label, y = area_km2,
                     colour = threshold_label,
                     group  = threshold_label)) +
  geom_line(linewidth = 0.9) +
  geom_point(data = area_df %>% filter(!is_bookend), size = 2.5) +
  geom_point(data = area_df %>% filter(is_bookend),
             size = 2.5, shape = 21, fill = "white") +
  geom_errorbar(
    data = area_df %>% filter(bootstrapped & !is.na(area_sd)),
    aes(ymin = area_km2 - area_sd, ymax = area_km2 + area_sd),
    width = 0.2, linewidth = 0.5
  ) +
  scale_colour_manual(
    values = setNames(thresh_settings$colour, thresh_settings$label),
    name   = "Threshold"
  ) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  period_xscale +
  labs(x = NULL,
       y = expression(bold("Suitable habitat area (km"^2*")"))) +
  shared_theme +
  theme(
    legend.position  = "right",
    legend.title     = element_text(size = 12),
    legend.text      = element_text(size = 11),
    legend.key.width = unit(1.2, "cm")
  )

# ============================================================================
# 7. BUILD PANEL B — NNI
# ============================================================================

cat("Building Panel B...\n")

p_nni <- ggplot(nni_df,
                aes(x = x_label, y = nni,
                    fill  = threshold,
                    group = threshold)) +
  geom_hline(yintercept = 1, linetype = "dashed",
             colour = "grey40", linewidth = 0.6) +
  geom_col(position = position_dodge(width = 0.75),
           width = 0.7, colour = "grey30", linewidth = 0.2) +
  # geom_text(aes(label = round(nni, 2), y = nni + 0.04),
  #           position = position_dodge(width = 0.75),
  #           size = 2.2, fontface = "bold", colour = "grey20") +
  scale_fill_manual(values = thresh_cols, name = "Threshold") +
  scale_y_continuous(limits = c(0, 1.6), breaks = seq(0, 1.5, 0.25)) +
  annotate("text", x = 1, y = 1.07,
           label = "random", hjust = 0.5, size = 3,
           colour = "grey40", fontface = "italic") +
  period_xscale +
  labs(x = "Period", y = "Nearest Neighbour Index (NNI)") +
  shared_theme +
  theme(
    legend.position  = "right",
    legend.title     = element_text(size = 12),
    legend.text      = element_text(size = 11),
    legend.key.width = unit(0.5, "cm")
  )

# ============================================================================
# 8. COMBINE
# ============================================================================

# p_combined <- p_area / p_nni +
#   plot_layout(heights = c(1, 1)) +
#   plot_annotation(tag_levels = "A") &
#   theme(plot.tag        = element_text(face = "bold", size = 12),
#         plot.tag.position = c(-0.015, 1),
#         plot.background = element_rect(fill = "white", colour = NA))

p_combined <- p_area / p_nni +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag          = element_text(face = "bold", size = 12),
        plot.tag.position = c(-0.015, 1),
        plot.background   = element_rect(fill = "white", colour = NA))

print(p_combined)


# ============================================================================
# SAVE — uncomment once happy
# ============================================================================

figures_dir <- "figures"
dir.create(figures_dir, showWarnings = FALSE)

ggsave(file.path(figures_dir, "habitat_nni.png"),
       p_combined, width = 8, height = 5.5, dpi = 600, units = "in",
       bg = "white")
ggsave(file.path(figures_dir, "habitat_nni.jpg"),
       p_combined, width = 10, height = 8, dpi = 600, units = "in",
        bg = "white")
