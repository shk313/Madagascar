######################################################
####### Time to Health Facilities in Madagascar 
####### Using Weiss et al Friction surface model
####################################################

###Time to each HF in madagascar ###

# clear workspace
rm(list = ls())

## Required Packages
require(gdistance)
library(raster)
library(rgdal)

#pull in polygon to get extent
poly <- readOGR(dsn = "Z:/Madagascar-NMCP/Madagascar Shapefiles/Madagascar_Admin_shapefiles/Malareo_District.shp", #
                layer = "Malareo_District")
e <- extent(poly)

# Input Files
friction.surface.filename <- "Z:/mastergrids/Other_Global_Covariates/Accessibility/Weiss/friction_surface_2015_v1.tif"
################### FILL IN BELOW ##############
point.filename <- "G:/Madagascar Research/"# structured as [UNIQUE_ID, X_COORD, Y_COORD] aka [LONG, LAT].  Use a header.

#  Define the spatial information from the friction surface
friction <- raster(friction.surface.filename)
fs1 <- crop(friction, e)
plot(fs1)
## Read in the points table.
points <- read.csv("lat_long_id.csv", header = TRUE)
head(points)

points <- subset(points, ID >2666)
head(points)
#Loop through points

HF_list <- unique(points$ID)
filepath <- ("C:/Users/SuzanneK/Documents/MAP stuff/Madagascar_Project/") ######## FILL IN HERE ######

for (i in seq_along(HF_list)) { 
  
  T.filename <- paste(filepath, HF_list[i], ".HF.T.rds", sep='')
  T.GC.filename <- paste(filepath, HF_list[i], ".HF.T.GC.rds", sep='')
  output.filename <- paste(filepath, HF_list[i], "HF.access.tif", sep='')
  
  T <- transition(fs1, function(x) 1/mean(x), 8) # RAM intensive, can be very slow for large areas
  saveRDS(T, T.filename)
  T.GC <- geoCorrection(T)                    
  saveRDS(T.GC, T.GC.filename)
  
  HF.coords <- c(points$X[i],points$Y[i])
  HF.raster <- accCost(T.GC, HF.coords)
  
  writeRaster(HF.raster, output.filename)
}
