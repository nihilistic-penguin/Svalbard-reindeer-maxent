library(raster)
library(terra)
library(ggplot2)
library(viridis)
library(dplyr)
library(sf)
library(ggnewscale)
library(rcartocolor)
library(colorspace)
library(cowplot)
library(ggh4x)  # for facet_wrap2 with colored strips
library(ggspatial)

setwd("C:/Users/marie/Documents/R/Barrett's R/Svalbard Reindeer Maxent Paper")

# ============================================================================
# PERIOD ORDER AND LABELS
# ============================================================================

period_order <- c(
  "weeks_09_10", "weeks_15b_16",
  "weeks_18a_18b", "weeks_20a_20b", "weeks_21_22",
  "weeks_23_24", "weeks_25a_25b", "weeks_27_28"
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

all_panel_labels <- c(period_labels[period_order], "All observations")

# Earth palette — one colour per period
earth_cols <- carto_pal(8, "Earth")
names(earth_cols) <- period_labels[period_order]

# Named vector for strip colours in right panel
strip_colors <- earth_cols  # same palette for strips

utm_crs <- "+proj=utm +zone=33 +datum=WGS84 +units=m"

# ============================================================================
# HELPERS
# ============================================================================

ll_to_utm_x <- function(lon, lat) {
  pts <- st_as_sf(data.frame(lon = lon, lat = lat),
                  coords = c("lon", "lat"), crs = 4326)
  st_coordinates(st_transform(pts, utm_crs))[, 1]
}

ll_to_utm_y <- function(lon, lat) {
  pts <- st_as_sf(data.frame(lon = lon, lat = lat),
                  coords = c("lon", "lat"), crs = 4326)
  st_coordinates(st_transform(pts, utm_crs))[, 2]
}

# ============================================================================
# LOAD TERRAIN RASTERS
# ============================================================================

cat("Loading elevation and hillshade rasters...\n")
elevation_r <- rast("data/prepared for analyses/predictors_og_1.bil")
hillshade_r <- rast("data/prepared for analyses/predictors_og_2.bil")

to_facet_df <- function(r, value_name, panel_levels) {
  df <- as.data.frame(r, xy = TRUE)
  colnames(df)[3] <- value_name
  df <- df[!is.na(df[[value_name]]), ]
  do.call(rbind, lapply(panel_levels, function(lbl) {
    df$period_label <- factor(lbl, levels = panel_levels)
    df
  }))
}

elev_df  <- to_facet_df(elevation_r, "elevation",  all_panel_labels)
shade_df <- to_facet_df(hillshade_r, "hillshade", all_panel_labels)

# ============================================================================
# LOAD SUITABILITY RASTERS
# ============================================================================

cat("Loading prediction rasters...\n")

all_cont <- list()
for (pn in period_order) {
  f <- file.path("output", pn, "prediction_continuous.tif")
  if (!file.exists(f)) { cat("  Missing:", pn, "\n"); next }
  r  <- raster(f)
  df <- as.data.frame(r, xy = TRUE)
  colnames(df)[3] <- "suitability"
  df <- df[!is.na(df$suitability), ]
  df$period_label <- factor(period_labels[pn], levels = all_panel_labels)
  all_cont[[pn]] <- df
  cat("  Loaded:", pn, "\n")
}
plot_df <- bind_rows(all_cont)

# ============================================================================
# LOAD OCCURRENCE POINTS
# ============================================================================

occ_list <- list()
for (pn in period_order) {
  f <- file.path("output", pn, "occurrence_points.csv")
  if (!file.exists(f)) next
  occ <- read.csv(f)
  colnames(occ)[1:2] <- c("x", "y")
  occ$period_label <- factor(period_labels[pn], levels = all_panel_labels)
  occ$period_name  <- period_labels[pn]  # for colour mapping
  occ_list[[pn]] <- occ[, c("x", "y", "period_label", "period_name")]
}
occ_df <- bind_rows(occ_list)

# All observations pooled — keep period_name for colour
all_occ_df <- occ_df
all_occ_df$period_label <- factor("All observations", levels = all_panel_labels)
all_occ_df <- dplyr::distinct(all_occ_df, x, y, period_name, .keep_all = TRUE)

# ============================================================================
# AXIS LABELS
# ============================================================================

lon_seq <- seq(15.25, 15.40, by = 0.1)
lat_seq <- seq(78.16, 78.24, by = 0.03)

x_breaks <- ll_to_utm_x(lon_seq, rep(mean(lat_seq), length(lon_seq)))
x_labels <- paste0(lon_seq, "\u00b0E")

y_breaks <- ll_to_utm_y(rep(mean(lon_seq), length(lat_seq)), lat_seq)
y_labels <- paste0(lat_seq, "\u00b0N")

# ============================================================================
# SHARED THEME
# ============================================================================

theme_map <- function() {
  theme_grey(base_size = 10) +
    theme(
      panel.background  = element_rect(fill = "white", colour = NA),
      panel.border      = element_rect(fill = NA, colour = "grey60", linewidth = 0.4),
      panel.grid.major  = element_line(colour = "grey70", linewidth = 0.3),
      panel.grid.minor  = element_blank(),
      strip.text        = element_text(face = "bold", size = 12),
      strip.background  = element_rect(fill = "grey82", colour = NA),
      axis.text.x       = element_text(size = 14),
      axis.text.y       = element_text(size = 14),
      axis.title.x      = element_text(size = 14, face = "bold"),
      axis.title.y      = element_text(size = 14, face = "bold"),
      legend.key.height = unit(1.2, "cm"),
      legend.key.width  = unit(0.4, "cm"),
      legend.title      = element_text(size = 12, vjust = 3),
      legend.text       = element_text(size = 11),
      plot.margin       = margin(5, 5, 5, 5)
    )
}

# ============================================================================
# LEFT PANEL — All observations coloured by period
# ============================================================================

shade_df_all <- shade_df %>%
  filter(period_label == "All observations") %>%
  mutate(period_label = factor("All observations", levels = all_panel_labels))

elev_df_all <- elev_df %>%
  filter(period_label == "All observations") %>%
  mutate(period_label = factor("All observations", levels = all_panel_labels))

# Define the correct chronological order explicitly, matching your
# biweekly period naming convention (adjust the exact strings to match
# whatever labels are actually in all_occ_df$period_name)
period_order <- c(
  "early March", "late March", "early April", "late April",
  "early May", "mid May", "late May",
  "early June", "late June",
  "mid July"
)

all_occ_df$period_name <- factor(all_occ_df$period_name, levels = period_order)

p_all_obs <- ggplot() +
  geom_raster(data = shade_df_all,
              aes(x = x, y = y, fill = hillshade),
              show.legend = FALSE) +
  scale_fill_gradient(low = "grey5", high = "grey95", na.value = "white") +
  ggnewscale::new_scale_fill() +
  geom_raster(data = elev_df_all,
              aes(x = x, y = y, fill = elevation),
              alpha = 0.4, show.legend = FALSE) +
  scale_fill_gradient(low = "snow2", high = "grey90", na.value = NA) +
  ggnewscale::new_scale_fill() +
  geom_point(data        = all_occ_df,
             aes(x = x, y = y, fill = period_name),
             inherit.aes = FALSE,
             colour      = "black",
             alpha       = 1,
             size        = 1.7,
             shape       = 21) +
  scale_fill_manual(
    name   = "Survey period",
    values = earth_cols
  ) +
  annotation_scale(
    location      = "br",
    width_hint    = 0.3,
    style         = "ticks",
    line_width    = 0.5,
    text_cex      = 0.8,
    unit_category = "metric"
  ) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  scale_y_continuous(breaks = y_breaks, labels = y_labels) +
  coord_equal(expand = FALSE) +
  labs(x = "Longitude", y = "") +
  theme_map() +
  guides(fill = guide_legend(override.aes = list(size = 4))) +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 16, face = "bold"),
    legend.text     = element_text(size = 14),
    legend.key      = element_rect(fill = "white", colour = NA),
    axis.title.x    = element_text(colour = "white"),
    panel.border    = element_blank()
  ) +
  # guides(fill = guide_legend(override.aes = list(size = 2.2))) + 
  theme(legend.key.height = unit(0.5, "cm"))

print(p_all_obs)



#################adding gpx track ############################################
library(sf)

# ============================================================================
# LOAD GPX TRACK
# ============================================================================
# GPX tracks are usually stored under the "tracks" layer; "track_points" also
# works if you want individual points instead of a continuous line
# Option 1: forward slashes (works fine on Windows despite looking Unix-y)
gpx_track <- st_read(
  "C:/Users/marie/Documents/R/Barrett's R/Svalbard Reindeer Maxent Paper/data/Bj_rndalen_survey_transect.gpx",
  layer = "tracks", quiet = TRUE
)
# Reproject to match your raster/plot CRS (UTM zone 33N, as used elsewhere)
gpx_track_utm <- st_transform(gpx_track, utm_crs)

library(dplyr)

# Extract line coordinates as a plain data frame (no sf/coord_sf involved)
track_coords <- st_coordinates(gpx_track_utm) %>%
  as.data.frame()

# MULTILINESTRING coordinates include an L1/L2 grouping column —
# check what columns you actually got:
head(track_coords)

# ============================================================================
# ADD TRACK TO p_all_obs
# ============================================================================
p_all_obs_track <- p_all_obs +
  geom_path(
    data        = track_coords,
    aes(x = X, y = Y, group = L1),
    inherit.aes = FALSE,
    colour      = "black",
    linewidth   = 0.6,
    alpha       = 0.9
  )

print(p_all_obs_track)

figures_dir <- "figures"
dir.create(figures_dir, showWarnings = FALSE)

ggsave(file.path(figures_dir, "all_observations.png"),
       p_all_obs_track, width = 5, height = 6, dpi = 600, units = "in", bg = "white")


# ============================================================================
# RIGHT PANELS — 8 models with coloured strip backgrounds
# ============================================================================

plot_df_models  <- plot_df  %>% filter(period_label != "All observations")
shade_df_models <- shade_df %>% filter(period_label != "All observations")
elev_df_models  <- elev_df  %>% filter(period_label != "All observations")

p_models <- ggplot() +
  geom_raster(data = shade_df_models,
              aes(x = x, y = y, fill = hillshade),
              show.legend = FALSE) +
  scale_fill_gradient(low = "grey5", high = "grey95", na.value = "white") +
  ggnewscale::new_scale_fill() +
  geom_raster(data = elev_df_models,
              aes(x = x, y = y, fill = elevation),
              alpha = 0.4, show.legend = FALSE) +
  scale_fill_gradient(low = "snow2", high = "grey90", na.value = NA) +
  ggnewscale::new_scale_fill() +
  geom_raster(data = plot_df_models,
              aes(x = x, y = y, fill = suitability),
              alpha = 1) +
  scale_fill_carto_c(palette  = "Teal",
                     name     = "Habitat\nsuitability",
                     limits   = c(0, 1),
                     direction = -1) +
  # Coloured strip backgrounds matching observation point colours
  facet_wrap2(
    ~ period_label,
    ncol   = 4,
    nrow   = 2,
    strip  = strip_themed(
      background_x = lapply(earth_cols, function(col)
        element_rect(fill = col, colour = NA))
    )
  ) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  scale_y_continuous(breaks = y_breaks, labels = y_labels) +
  coord_equal(expand = FALSE) +
  labs(x = "Longitude", y = NULL) +
  theme_map() +
  theme(
    axis.text.y    = element_blank(),
    axis.ticks.y   = element_blank(),
    axis.title.x   = element_text(hjust = 0.1),
    legend.position = "right",
    panel.spacing = unit(0.3, "cm")
  )

p_models

p_models_new <- ggplot() +
  geom_raster(data = shade_df_models,
              aes(x = x, y = y, fill = hillshade),
              show.legend = FALSE) +
  scale_fill_gradient(low = "grey5", high = "grey95", na.value = "white") +
  ggnewscale::new_scale_fill() +
  geom_raster(data = elev_df_models,
              aes(x = x, y = y, fill = elevation),
              alpha = 0.4, show.legend = FALSE) +
  scale_fill_gradient(low = "snow2", high = "grey90", na.value = NA) +
  ggnewscale::new_scale_fill() +
  geom_raster(data = plot_df_models,
              aes(x = x, y = y, fill = suitability),
              alpha = 1) +
  scale_fill_carto_c(palette  = "Teal",
                     name     = "Habitat\nsuitability",
                     limits   = c(0, 1),
                     direction = -1) +
  # Coloured strip backgrounds matching observation point colours
  facet_wrap2(
    ~ period_label,
    ncol   = 4,
    nrow   = 2,
    strip  = strip_themed(
      background_x = lapply(earth_cols, function(col)
        element_rect(fill = col, colour = NA))
    )
  ) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  scale_y_continuous(breaks = y_breaks, labels = y_labels) +
  coord_equal(expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_map() +
  theme(
    legend.position = "right",
    panel.spacing = unit(0.3, "cm"), 
    panel.border = element_blank()
  )

p_models_new


figures_dir <- "figures"
dir.create(figures_dir, showWarnings = FALSE)

ggsave(file.path(figures_dir, "suitability_models_alone.png"),
       p_models_new, width = 8, height = 8, dpi = 600, units = "in", bg = "white")

# ============================================================================
# COMBINE
# ============================================================================

# Fix B label position
p_final <- plot_grid(
  p_all_obs,
  p_models,
  ncol           = 2,
  rel_widths     = c(0.93, 1.6),
  align          = "h",
  axis           = "tb",
  labels         = c("A", "B"),
  label_size     = 14,
  label_fontface = "bold",
  label_x        = c(0, -0.02),  # nudge B slightly right
  label_y        = c(1, 1)
)

print(p_final)

# ============================================================================
# SAVE
# ============================================================================

figures_dir <- "figures"
dir.create(figures_dir, showWarnings = FALSE)

ggsave(file.path(figures_dir, "FINAL_MAIN_suitability_all_periods.png"),
       p_final, width = 9, height = 6, dpi = 600, units = "in", bg = "white")

ggsave(file.path(figures_dir, "FINAL_MAIN_suitability_all_periods.tiff"),
       p_final, width = 9, height = 6, dpi = 600, units = "in", bg = "white",
       compression = "lzw")

cat("Saved to figures/\n")
