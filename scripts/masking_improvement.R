#### this is a HELPER SCRIPT to find the right UTMs for the cliff masking script ###

library(raster)
library(terra)

# Load your raster
r <- raster("data/prepared for analyses/predictors_masked.grd")

# Open interactive plot — click on the corners of the areas to exclude
plot(r[[1]])

# Then run this and click points along your boundary lines
# Press Escape when done
clicked <- click(r[[1]], n = 20, id = TRUE, xy = TRUE, cell = TRUE)
print(clicked[, c("x", "y")])


extent(r)

cat("xmin:", extent(r[[1]])@xmin, "\n")
cat("xmax:", extent(r[[1]])@xmax, "\n")
cat("ymin:", extent(r[[1]])@ymin, "\n")
cat("ymax:", extent(r[[1]])@ymax, "\n")

plot(predictors[[1]])
clicked3 <- click(predictors[[1]], n = 10, id = TRUE, xy = TRUE, cell = TRUE)
print(clicked3[, c("x", "y")])