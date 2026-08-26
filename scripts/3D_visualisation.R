library(terra)
library(raster)
library(plotly)
library(scales)
library(htmlwidgets)
library(rcartocolor)

setwd("C:/Users/marie/Documents/R/Barrett's R/Svalbard Reindeer Maxent Paper")

# ============================================================================
# PERIOD DEFINITIONS
# ============================================================================

period_order <- c(
  "weeks_09_10", "weeks_15b_16",
  "weeks_18a_18b", "weeks_20a_20b", "weeks_21_22",
  "weeks_23_24", "weeks_25a_25b", "weeks_27_28"
)

period_labels <- c(
  weeks_09_10   = "Early March",
  weeks_15b_16  = "Late April",
  weeks_18a_18b = "Early May",
  weeks_20a_20b = "Mid May",
  weeks_21_22   = "Late May",
  weeks_23_24   = "Early June",
  weeks_25a_25b = "Late June",
  weeks_27_28   = "Mid July"
)

# ============================================================================
# LOAD DEM
# ============================================================================

cat("Loading DEM...\n")
dem <- rast("data/prepared for analyses/predictors_og_1.bil")

# Aggregate to reduce rendering load — increase factor if slow
dem_agg  <- aggregate(dem, fact = 3, fun = mean)
dem_mat  <- as.matrix(raster(dem_agg))
dem_mat  <- dem_mat[nrow(dem_mat):1, ]  # flip for plotly orientation

# Get x/y coordinates
dem_ext  <- ext(dem_agg)
x_coords <- seq(dem_ext$xmin, dem_ext$xmax, length.out = ncol(dem_mat))
y_coords <- seq(dem_ext$ymin, dem_ext$ymax, length.out = nrow(dem_mat))

# ============================================================================
# FUNCTION TO BUILD ONE 3D PLOTLY FIGURE
# ============================================================================
earth_cols <- carto_pal(8, "Earth")
names(earth_cols) <- period_labels[period_order]

make_3d_plot <- function(period_name, period_label) {
  
  cat("Processing:", period_label, "\n")
  
  # Load suitability raster
  suit_file <- file.path("output", period_name, "prediction_continuous.tif")
  if (!file.exists(suit_file)) {
    cat("  Missing:", period_name, "\n")
    return(NULL)
  }
  
  suit <- rast(suit_file)
  suit_agg <- resample(suit, dem_agg, method = "bilinear")
  suit_mat <- as.matrix(raster(suit_agg))
  suit_mat <- suit_mat[nrow(suit_mat):1, ]
  suit_mat[is.na(suit_mat)] <- 0
  dem_mat_clean <- dem_mat
  dem_mat_clean[is.na(dem_mat_clean)] <- 0
  ocean_z <- matrix(0, nrow = nrow(dem_mat_clean), ncol = ncol(dem_mat_clean))
  
  # Load occurrence points
  occ_file <- file.path("output", period_name, "occurrence_points.csv")
  occ <- read.csv(occ_file)
  colnames(occ)[1:2] <- c("x", "y")
  
  # Extract elevation at each occurrence point for z coordinate
  dem_raster <- raster(dem_agg)
  occ_z <- raster::extract(dem_raster, occ[, c("x", "y")])
  occ_z[is.na(occ_z)] <- 0
  occ_z <- occ_z + 15  # offset slightly above terrain so points are visible
  
  # Get point colour from Earth palette
  point_col <- earth_cols[period_labels[period_name]]
  
  fig <- plot_ly(height = 800) %>%
    add_surface(
      x            = x_coords,
      y            = y_coords,
      z            = dem_mat_clean,
      surfacecolor = suit_mat,
      colorscale   = list(
        list(0,    "#00334d"),
        list(0.25, "#00667a"),
        list(0.5,  "#4da6a8"),
        list(0.75, "#99d4d4"),
        list(1,    "#e0f5f5")
      ),
      cmin     = 0,
      cmax     = 1,
      colorbar = list(
        title     = "Habitat<br>suitability",
        titlefont = list(size = 14),
        tickfont  = list(size = 12),
        len       = 0.5
      ),
      contours = list(z = list(show = FALSE)),
      lighting = list(
        ambient   = 1.5,   # reduce from 0.9 — more shadow contrast
        diffuse   = 1.2,   # increase — stronger directional lighting
        specular  = 0.5,
        roughness = 0.9,   # increase — more matte surface shows relief better
        fresnel   = 0.1
      ),
      lightposition = list(x = 2, y = 2, z = 3)  # higher light position
    ) %>%
    add_trace(
      type       = "scatter3d",
      mode       = "markers",
      x          = occ$x,
      y          = occ$y,
      z          = occ_z,
      marker     = list(
        size    = 4,
        color   = point_col,
        opacity = 0.85,
        line    = list(color = "black", width = 0.5)
      ),
      name       = period_label,
      showlegend = FALSE
    )  %>%
    layout(
      title = list(
        text = paste0("<b>", period_label, " — Habitat Suitability</b>"),
        font = list(size = 16)
      ),
      scene = list(
        xaxis = list(title = "Easting (m)", titlefont = list(size = 12)),
        yaxis = list(title = "Northing (m)", titlefont = list(size = 12)),
        zaxis = list(title = "Elevation (m)", titlefont = list(size = 12)),
        camera = list(
          #eye = list(x = -1.5, y = -1.5, z = 1.2)
          eye = list(x = -0.8, y = 2.2, z = 0.6)),
        aspectmode  = "manual",
        aspectratio = list(x = 1, y = 2, z = 0.45)
      ),
      margin        = list(l = 0, r = 0, t = 50, b = 0),
      paper_bgcolor = "white",
      height        = 800
    )
  
  return(fig)
}

# ============================================================================
# GENERATE AND SAVE ALL PERIODS
# ============================================================================

output_dir <- "figures/3d_supplementary"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Make plotly available offline
plotly::plotly_IMAGE # not needed, just use:
options(htmltools.dir.version = FALSE)

for (pn in period_order) {
  
  fig <- make_3d_plot(pn, period_labels[pn])
  
  if (!is.null(fig)) {
    # Use period label for filename, replacing spaces with underscores
    file_label <- gsub(" ", "_", tolower(period_labels[pn]))
    out_file   <- file.path(output_dir,
                            paste0(file_label, "_3d_suitability.html"))
    htmltools::save_html(
      htmltools::browsable(
        htmltools::tagList(
          htmltools::tags$head(
            htmltools::tags$style("body { margin: 0; } .plotly { width: 100vw !important; height: 100vh !important; }")
          ),
          plotly::as_widget(fig)
        )
      ),
      file = out_file
    )
    saveWidget(fig, file = out_file, selfcontained = TRUE)
  }
}

cat("\nAll 3D figures saved to", output_dir, "\n")

# Delete all folders in the output directory, keep only HTML files
dirs_to_remove <- list.dirs(output_dir, recursive = FALSE)
unlink(dirs_to_remove, recursive = TRUE)
cat("Cleaned up", length(dirs_to_remove), "dependency folders\n")
