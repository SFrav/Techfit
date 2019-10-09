library(raster)
#library(sp)
library(sf)
library(rgdal)

setwd("D:/GIS/Copernicus/fAPAR") #To location with fAPAR sub-folders
#create EA folder before running loop

Maskshp <- readOGR("aoi.shp")

for (i in 1:length(list.files(full.names=TRUE, pattern = "PROBAV_V1.0.1.nc", recursive = TRUE))) {
  
  files <- list.files(full.names=TRUE, pattern = "PROBAV_V1.0.1.nc", recursive = TRUE)
  fAPAR_raster <- raster(files[i], varname = "FAPAR")
  crs(fAPAR_raster) <- crs(Maskshp) #Copernicus fAPAR and DMP are in regular latitude/longitude grid WGS84

  filename <- (paste("./EA/", tools::file_path_sans_ext(basename(files[i])), ".tif", sep=""))
  
  # Crop the raster
  fAPAR_raster.crop <- crop(fAPAR_raster, extent(Maskshp))
  fAPAR_raster.crop <- mask(fAPAR_raster.crop, Maskshp)

  #Export file 
  writeRaster(fAPAR_raster.crop, filename)
}
