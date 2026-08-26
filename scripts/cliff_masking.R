# ============================================================================
# CLIFF MASKING SCRIPT
# Removes two cliff areas from the masked predictor raster:
#   - Right cliff: everything to the right of line (507480,8684960)-(508380,8684320)
#   - Left cliff:  everything to the left of line  (505880,8683020)-(506360,8682360)
# Saves updated raster as predictors_masked2.grd
# ============================================================================

library(raster)
library(sf)

setwd("C:/Users/marie/Documents/R/Barrett's R/Svalbard Reindeer Maxent Paper")

utm_crs <- "+proj=utm +zone=33 +datum=WGS84 +units=m"

# Raster extent (from your output)
xmin <- 505010; xmax <- 509010
ymin <- 8675970; ymax <- 8685010

# ============================================================================
# 1. LOAD CURRENT MASKED RASTER
# ============================================================================

predictors <- stack("data/prepared for analyses/predictors_masked.grd")
cat("Loaded raster:", nlayers(predictors), "layers\n")

# Preview before masking
plot(predictors[[1]], main = "Before additional masking")

# ============================================================================
# 2. BUILD MASKING POLYGONS
# ============================================================================

# --- RIGHT CLIFF ---
# Line goes from (507480, 8684960) to (508380, 8684320)
# Mask everything to the RIGHT of this line
# Close the polygon using the top-right and right raster edges

right_cliff <- st_polygon(list(matrix(c(
  507480, 8684960,   # start of cut line (top point)
  508380, 8684320,   # end of cut line (bottom point)
  xmax,   8684320,   # extend to right raster edge at same y
  xmax,   ymax,      # top-right corner of raster
  507480, ymax,      # across top edge to above start point
  507480, 8684960    # close polygon
), ncol = 2, byrow = TRUE)))

# --- LEFT CLIFF / ARM ---
# Traces the outline of the arm and closes via the left raster edge

left_cliff <- st_polygon(list(matrix(c(
  505840, 8683020,   # top of arm
  506320, 8682220,   # upper right of arm
  505440, 8681000,   # mid right of arm
  505720, 8679920,   # bottom of arm   
  xmin,   8679440,   # extend to left raster edge at bottom
  xmin,   8683020,   # up the left edge to top
  505840, 8683020    # close polygon
), ncol = 2, byrow = TRUE)))

# Combine into one sf object
mask_polys <- st_sfc(right_cliff, left_cliff, crs = utm_crs)
mask_sf    <- st_sf(geometry = mask_polys)

# ============================================================================
# 3. PREVIEW POLYGONS ON THE RASTER — check before applying
# ============================================================================

plot(predictors[[1]], main = "Check masking polygons (red = areas to remove)")
plot(st_geometry(mask_sf), add = TRUE, col = adjustcolor("red", alpha = 0.4),
     border = "red", lwd = 2)

cat("\nCheck the plot — do the red polygons cover the correct cliff areas?\n")
cat("If yes, run the rest of the script to apply the mask.\n")
cat("If not, adjust the polygon coordinates above and re-run sections 2-3.\n\n")

# ============================================================================
# STOP HERE AND CHECK THE PLOT BEFORE CONTINUING
# Highlight lines below and run only when happy with the polygons
# ============================================================================

# ============================================================================
# 4. APPLY MASK TO ALL RASTER LAYERS
# ============================================================================

# Rasterize the mask polygons (1 = area to remove)
mask_r <- rasterize(as(mask_sf, "Spatial"), predictors[[1]], field = 1)

# Set masked cells to NA in all layers
predictors_masked2 <- mask(predictors, mask_r, inverse = TRUE)

cat("Masking applied.\n")
cat("Original NA cells:  ", cellStats(is.na(predictors[[1]]), sum), "\n")
cat("New NA cells:       ", cellStats(is.na(predictors_masked2[[1]]), sum), "\n")
cat("Additional cells masked:", 
    cellStats(is.na(predictors_masked2[[1]]), sum) - 
      cellStats(is.na(predictors[[1]]), sum), "\n")

# ============================================================================
# 5. PREVIEW RESULT
# ============================================================================

plot(predictors_masked2[[1]], main = "After additional masking")
plot(st_geometry(mask_sf), add = TRUE, border = "red", lwd = 1, lty = 2)

# ============================================================================
# 6. SAVE — run once happy with the result
# ============================================================================

writeRaster(
  predictors_masked2,
  filename  = "data/prepared for analyses/predictors_masked2.grd",
  format    = "raster",
  overwrite = TRUE
)
cat("Saved: data/prepared for analyses/predictors_masked2.grd\n")
cat("Update your main script to load predictors_masked2.grd instead\n")
