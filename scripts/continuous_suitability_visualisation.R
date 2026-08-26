# ============================================================================
# CONTINUOUS SUITABILITY — 3x3 PANEL FIGURE
# Excludes weeks_11_12 (INSUFFICIENT)
# Marks weeks_13_15a (POOR) with label in lower right corner
# Terrain background: hillshade + elevation in greyscale
# Suitability: Geyser (rcartocolor) palette
# ============================================================================

library(raster)
library(terra)
library(ggplot2)
library(viridis)
library(dplyr)
library(sf)
library(ggnewscale)
library(rcartocolor)

setwd("C:/Users/marie/Documents/R/Barrett's R/Svalbard Reindeer Maxent Paper")

# ============================================================================
# PERIOD ORDER AND LABELS
# weeks_11_12 excluded (INSUFFICIENT)
# ============================================================================

period_order <- c(
  "weeks_09_10", "weeks_13_15a", "weeks_15b_16",
  "weeks_18a_18b", "weeks_20a_20b", "weeks_21_22",
  "weeks_23_24", "weeks_25a_25b", "weeks_27_28"
)

period_labels <- c(
  weeks_09_10   = "Early March",
  weeks_13_15a  = "Early April",
  weeks_15b_16  = "Late April",
  weeks_18a_18b = "Early May",
  weeks_20a_20b = "Mid May",
  weeks_21_22   = "Late May",
  weeks_23_24   = "Early June",
  weeks_25a_25b = "Late June",
  weeks_27_28   = "Mid July"
)

quality_labels <- data.frame(
  period_label = factor("Early April", levels = period_labels[period_order]),
  label        = "POOR FIT"
)

utm_crs <- "+proj=utm +zone=33 +datum=WGS84 +units=m"

# ============================================================================
# HELPERS: UTM to lon/lat for axis labels
# ============================================================================

utm_to_lon <- function(x_utm, y_utm_mid) {
  pts <- st_as_sf(data.frame(x = x_utm, y = y_utm_mid),
                  coords = c("x", "y"), crs = utm_crs)
  st_coordinates(st_transform(pts, 4326))[, 1]
}

utm_to_lat <- function(x_utm_mid, y_utm) {
  pts <- st_as_sf(data.frame(x = x_utm_mid, y = y_utm),
                  coords = c("x", "y"), crs = utm_crs)
  st_coordinates(st_transform(pts, 4326))[, 2]
}

# ============================================================================
# LOAD TERRAIN RASTERS (unmasked, for background)
# ============================================================================

cat("Loading elevation and hillshade rasters...\n")
elevation_r <- rast("data/prepared for analyses/predictors_og_1.bil")
hillshade_r <- rast("data/prepared for analyses/predictors_og_2.bil")

to_facet_df <- function(r, value_name) {
  df <- as.data.frame(r, xy = TRUE)
  colnames(df)[3] <- value_name
  df <- df[!is.na(df[[value_name]]), ]
  do.call(rbind, lapply(period_labels[period_order], function(lbl) {
    df$period_label <- factor(lbl, levels = period_labels[period_order])
    df
  }))
}

elev_df  <- to_facet_df(elevation_r, "elevation")
shade_df <- to_facet_df(hillshade_r, "hillshade")
cat("  Terrain layers loaded and tiled across", length(period_order), "panels\n")

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
  df$period_label <- factor(period_labels[pn], levels = period_labels[period_order])
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
  occ$period_label <- factor(period_labels[pn], levels = period_labels[period_order])
  occ_list[[pn]] <- occ[, c("x", "y", "period_label")]
}
occ_df <- bind_rows(occ_list)

# ============================================================================
# AXIS LABELS
# ============================================================================

x_mid    <- mean(range(plot_df$x))
y_mid    <- mean(range(plot_df$y))
x_breaks <- pretty(range(plot_df$x), n = 2)
y_breaks <- pretty(range(plot_df$y), n = 4)
x_labels <- paste0(round(utm_to_lon(x_breaks, rep(y_mid, length(x_breaks))), 2), "\u00b0E")
y_labels <- paste0(round(utm_to_lat(rep(x_mid, length(y_breaks)), y_breaks), 2), "\u00b0N")

# ============================================================================
# BUILD FIGURE
# ============================================================================

cat("\nBuilding figure...\n")

library(colorspace)
earth_cols <- rev(carto_pal(7, "Earth"))
darker_cols <- darken(earth_cols, amount = 0.15)  # 0-1, higher = darker



p <- ggplot() +
  # Layer 1: hillshade as greyscale base
  geom_raster(data = shade_df,
              aes(x = x, y = y, fill = hillshade),
              show.legend = FALSE) +
  scale_fill_gradient(low = "grey5", high = "grey95",
                      na.value = "white") +
  # Layer 2: elevation overlaid with transparency
  ggnewscale::new_scale_fill() +
  geom_raster(data = elev_df,
              aes(x = x, y = y, fill = elevation),
              alpha = 0.4,
              show.legend = FALSE) +
  scale_fill_gradient(low = "snow2", high = "grey90",
                      na.value = NA) +
  # Layer 3: suitability in Geyser
  ggnewscale::new_scale_fill() +
  geom_raster(data = plot_df,
              aes(x = x, y = y, fill = suitability),
              alpha = 0.9) +
  # geom_point(data        = occ_df,
  #            aes(x = x, y = y),
  #            inherit.aes = FALSE,
  #            colour      = "grey",
  #            fill = "white" ,
  #            alpha       = 0.5,
  #            size        = 1,
  #            shape       = 21) +
  geom_text(data        = quality_labels,
            aes(x       = quantile(plot_df$x, 0.97),
                y       = quantile(plot_df$y, 0.03),
                label   = label),
            inherit.aes = FALSE,
            colour      = "darkred",
            size        = 2.8,
            fontface    = "bold",
            hjust       = 1,
            vjust       = 0) +
  scale_fill_gradientn(
    colours = darker_cols,
    name    = "Suitability",
    limits  = c(0, 1)
  ) +
  # scale_fill_carto_c(
  #   palette   = "Earth",
  #   name      = "Suitability",
  #   direction = -1,
  #   limits    = c(0, 1),   # shift lower limit below 0 to darken
  #   breaks    = c(0, 0.25, 0.5, 0.75, 1),
  #   labels    = c("0", "0.25", "0.50", "0.75", "1")
  # ) +
  # scale_fill_viridis_c(
  #   option    = "cividis",
  #   name      = "Suitability",
  #   direction = 1,
  #   limits    = c(0, 1),
  #   breaks    = c(0, 0.25, 0.5, 0.75, 1),
  #   labels    = c("0", "0.25", "0.50", "0.75", "1")
  # ) +
  facet_wrap(~ period_label, ncol = 5
            # , nrow = 3
             ) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  scale_y_continuous(breaks = y_breaks, labels = y_labels) +
  coord_equal(expand = FALSE) +
  labs(
    #title = "Habitat suitability — Svalbard reindeer 2023",
    x     = "Longitude",
    y     = "Latitude"
  ) +
  theme_grey(base_size = 10) +
  theme(
    strip.text        = element_text(face = "bold", size = 9),
    strip.background  = element_rect(colour = "grey"),
    axis.text.x       = element_text(size = 7, angle = 30, hjust = 1),
    axis.text.y       = element_text(size = 7),
    axis.title        = element_text(size = 9),
    panel.grid.major  = element_line(colour = "grey88", linewidth = 0.3),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(colour = "grey70", fill = NA,
                                     linewidth = 0.4),
    legend.position   = "right",
    legend.key.height = unit(1.2, "cm"),
    legend.key.width  = unit(0.4, "cm"),
    legend.title      = element_text(size = 9),
    legend.text       = element_text(size = 8),
    plot.title        = element_text(face = "bold", size = 13, hjust = 0),
    plot.margin       = margin(10, 10, 10, 10)
  )

print(p)

# ============================================================================
# SAVE — uncomment once happy with the figure
# ============================================================================

# figures_dir <- "figures"
# dir.create(figures_dir, showWarnings = FALSE)
#
# ggsave(file.path(figures_dir, "suitability_all_periods.png"),
#        p, width = 12, height = 12, dpi = 300, units = "in")
#
# ggsave(file.path(figures_dir, "suitability_all_periods.tiff"),
#        p, width = 12, height = 12, dpi = 300, units = "in",
#        compression = "lzw")
#
# cat("Saved to figures/\n")