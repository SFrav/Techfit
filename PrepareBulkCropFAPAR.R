library(raster)
#library(sp)
library(sf)
library(rgdal)

setwd("D:/GIS/Copernicus/fAPAR") #To location with fAPAR sub-folders

Maskshp <- readOGR("aoi.shp")

# Getting the spatial extent of the shapefile
e <- extent(Maskshp)

for (i in 1:length(list.files(full.names=TRUE, pattern = ".nc", recursive = TRUE))) {
  
  files <- list.files(full.names=TRUE, pattern = ".nc", recursive = TRUE)
  fAPAR_raster <- raster(files[i], varname = "FAPAR")
  crs(fAPAR_raster) <- crs(Maskshp)

  filename <- (paste("./EA/", tools::file_path_sans_ext(basename(files[i])), ".tif", sep=""))
  
  # Crop the raster
  fAPAR_raster.crop <- crop(fAPAR_raster, extent(Maskshp))
  fAPAR_raster.crop <- mask(fAPAR_raster.crop, Maskshp)

  #Export file 
  writeRaster(Maskshp.masked, filename)
}
