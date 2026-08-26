# ============================================================================
# PREDICTOR VARIABLES — FACET FIGURE (per-panel legends, own scales)
# ============================================================================
library(raster)
library(terra)
library(ggplot2)
library(dplyr)
library(sf)
library(rcartocolor)
library(colorspace)
library(patchwork)

setwd("C:/Users/marie/Documents/R/Barrett's R/Svalbard Reindeer Maxent Paper")

# ============================================================================
# PREDICTOR LABELS
# ============================================================================
predictor_names <- c(
  "elevation", "hillshade", "ruggedness",
  "landforms_fc", "landcover_fc", "NDVI22.23"
)
predictor_labels <- c(
  elevation    = "Elevation",
  hillshade    = "Hillshade",
  ruggedness   = "Ruggedness",
  landforms_fc = "Landform class",
  landcover_fc = "Landcover class",
  NDVI22.23    = "NDVI"
)
categorical_vars <- c("landforms_fc", "landcover_fc")
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
# LOAD PREDICTOR RASTERS — RAW VALUES, NO NORMALISATION
# ============================================================================
bil_files <- c(
  elevation    = "data/prepared for analyses/predictors_og_1.bil",
  hillshade    = "data/prepared for analyses/predictors_og_2.bil",
  ruggedness   = "data/prepared for analyses/predictors_og_3.bil",
  landforms_fc = "data/prepared for analyses/predictors_og_4.bil",
  landcover_fc = "data/prepared for analyses/predictors_og_5.bil",
  NDVI22.23    = "data/prepared for analyses/predictors_og_6.bil"
)

cat("Loading predictor rasters (raw, unnormalised)...\n")
pred_list_raw <- list()
for (pn in predictor_names) {
  r  <- rast(bil_files[pn])
  df <- as.data.frame(r, xy = TRUE)
  colnames(df)[3] <- "value"
  df <- df[!is.na(df$value), ]
  
  if (pn %in% categorical_vars) {
    df$value <- as.factor(df$value)
  }
  
  df$panel_label <- factor(predictor_labels[pn],
                           levels = predictor_labels[predictor_names])
  pred_list_raw[[pn]] <- df
  cat("  Loaded:", pn, "\n")
}

# ============================================================================
# SHARED THEME
# ============================================================================
shared_theme <- theme_grey(base_size = 10) +
  theme(
    strip.text        = element_blank(),
    axis.text.x       = element_text(size = 9),
    axis.text.y       = element_text(size = 9),
    axis.title        = element_blank(),
    panel.grid.major  = element_line(colour = "grey70", linewidth = 0.3),
    panel.grid.minor  = element_blank(),
    panel.border      = element_blank(),
    panel.background  = element_rect(fill = "white"),
    plot.title        = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.margin       = margin(5, 5, 5, 5)
  )

# ============================================================================
# CATEGORY LABEL LOOKUPS
# ============================================================================
landform_labels <- c(
  "1" = "Valley",
  "2" = "Lower slope",
  "3" = "Flat area",
  "4" = "Middle slope",
  "5" = "Upper slope",
  "6" = "Ridge"
)

landcover_labels <- c(
  "1" = "Ocean",
  "2" = "Freshwater",
  "3" = "Shadow",
  "4" = "Barren",
  "5" = "Mesic",
  "6" = "Wetland",
  "7" = "Bird cliff",
  "8" = "Dryas",
  "9" = "Cassiope"
)

category_lookups <- list(
  landforms_fc = landform_labels,
  landcover_fc = landcover_labels
)

# ============================================================================
# PLOT-BUILDING FUNCTION — updated: named categories, no axis ticks/labels
# ============================================================================
make_predictor_plot <- function(df, label, is_categorical = FALSE, var_name = NULL) {
  
  if (is_categorical) {
    df$value <- as.factor(df$value)
    
    # Relabel factor levels using the lookup table, if provided
    if (!is.null(var_name) && var_name %in% names(category_lookups)) {
      lookup <- category_lookups[[var_name]]
      levels(df$value) <- lookup[levels(df$value)]
    }
    
    p <- ggplot(df, aes(x = x, y = y, fill = value)) +
      geom_raster() +
      scale_fill_carto_d(
        name    = NULL,
        palette = "Earth"
      )
  } else {
    df$value <- as.numeric(as.character(df$value))
    
    p <- ggplot(df, aes(x = x, y = y, fill = value)) +
      geom_raster() +
      scale_fill_carto_c(
        name      = NULL,
        palette   = "Earth",
        direction = -1
      )
  }
  
  p +
    coord_fixed(ratio = 1, expand = FALSE) +
    labs(title = label, x = "", y = "") +
    shared_theme +
    theme(
      axis.text.x             = element_blank(),
      axis.text.y             = element_blank(),
      axis.ticks              = element_blank(),
      panel.grid.major        = element_blank(),   # graticule lines were tied to axis breaks; drop too
      legend.position          = "inside",
      legend.position.inside   = c(0.95, 0.05),
      legend.justification     = c("right", "bottom"),
      legend.background        = element_rect(fill = alpha("white", 0.65), colour = NA),
      legend.key.height        = unit(0.35, "cm"),
      legend.key.width         = unit(0.2, "cm"),
      legend.title             = element_blank(),
      legend.text              = element_text(size = 7)
    )
}

# ============================================================================
# BUILD ALL SIX PLOTS — now passing var_name for the lookup
# ============================================================================
predictor_plots <- lapply(predictor_names, function(pn) {
  make_predictor_plot(
    df             = pred_list_raw[[pn]],
    label          = predictor_labels[pn],
    is_categorical = pn %in% categorical_vars,
    var_name       = pn
  )
})

p_all <- wrap_plots(predictor_plots, ncol = 6) +
  plot_layout(axis_titles = "collect")

print(p_all)

# TEMPORARY — test bottom_grid sizing alone
ggsave("figures/predictors_LINE.png",
       p_all,
       width = 9, height = 4,
       dpi = 600, units = "in", bg = "white")

# ============================================================================
# COMBINE WITH LOCATION MAP + ALL-OBSERVATIONS MAP (A / B / C layout)
# ============================================================================
library(cowplot)

# ============================================================================
# COMPUTE ACTUAL DATA ASPECT RATIO (width:height of the valley extent)
# ============================================================================
bbox <- range(pred_list_raw$elevation$x) # xmin/xmax
bbox_y <- range(pred_list_raw$elevation$y)

data_width  <- diff(bbox)
data_height <- diff(bbox_y)
aspect_ratio <- data_width / data_height   # e.g. ~0.35 for a narrow tall valley

cat("Data aspect ratio (w/h):", aspect_ratio, "\n")

# ============================================================================
# BOTTOM GRID — six predictor panels (C-H), 2 rows x 3 cols
# ============================================================================
bottom_grid <- plot_grid(
  plotlist = predictor_plots,   # list of 6 ggplot objects, in C-H order
  ncol     = 3,
  align    = "hv",
  axis     = "tblr"
)

# ============================================================================
# TOP ROW — location map (A) + all-observations map (B)
# ============================================================================
top_row <- plot_grid(
  p_location,
  p_all_obs_track,
  ncol       = 2,
  rel_widths = c(1, 0.7),
  align      = "hv",
  axis       = "tblr"
)

n_col <- 3
n_row <- 2
panel_w <- 2.2          # inches — starting guess, tune to taste
panel_h <- panel_w / aspect_ratio   # 2.2 / 0.4412 ≈ 4.99 in

bottom_grid_width  <- n_col * panel_w    # 6.6 in
bottom_grid_height <- n_row * panel_h    # ~9.98 in

cat("Bottom grid:", bottom_grid_width, "x", bottom_grid_height, "in\n")
# ============================================================================
# FULL FIGURE — top row + bottom grid stacked, no labels for now
# ============================================================================
bottom_grid <- plot_grid(
  plotlist = predictor_plots,
  ncol     = 3
)

top_row <- plot_grid(
  p_location,
  p_all_obs_track,
  ncol       = 2,
  rel_widths = c(1, 0.7),
  align      = "hv",
  axis       = "tblr"
)

p_final <- plot_grid(
  top_row,
  bottom_grid,
  ncol        = 1,
  rel_heights = c(0.8, 1)
)

print(p_final)


# TEMPORARY — test bottom_grid sizing alone
ggsave("figures/predictors_grid.png",
       bottom_grid,
       width = bottom_grid_width, height = bottom_grid_height,
       dpi = 600, units = "in", bg = "white")

# ============================================================================
ggsave(file.path(figures_dir, "NEW_map_predictor_panel.png"),
       p_final,
       width  = bottom_grid_width + 1,
       height = bottom_grid_height + 4,
       dpi = 600, units = "in", bg = "white")


ggsave(file.path(figures_dir, "FINAL_map_predictor_panel.png"),
       p_final,
       width  = bottom_grid_width + 1,     # + buffer for legends/margins
       height = bottom_grid_height + 4,    # + top row (A/B) height
       dpi    = 600, units = "in", bg = "white")

ggsave(file.path(figures_dir, "FINAL_map_predictor_panel.tiff"),
       p_final,
       width  = bottom_grid_width + 1,
       height = bottom_grid_height + 4,
       dpi    = 600, units = "in", bg = "white",
       compression = "lzw")