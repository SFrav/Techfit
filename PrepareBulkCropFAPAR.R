library(raster)
library(sf)
library(rgdal)

iDir <-"D:/GIS/Copernicus/fAPAR" #To location with fAPAR sub-folders

Maskshp <- readOGR("aoi.shp")

for (i in 1:length(list.files(paste0(iDir, "/"), full.names = TRUE, pattern = ".nc", recursive = TRUE))) {
  
  EA <- paste0(iDir, "/EA/", sep = "")
  if (!file.exists(EA)) {dir.create(EA, recursive = T)}
  
  files <- list.files(paste0(iDir, "/"), full.names = TRUE, pattern = ".nc", recursive = TRUE)
  
  fAPAR_raster <- raster(files[i], varname = "FAPAR")
  
  crs(fAPAR_raster) <- crs(Maskshp) #Copernicus fAPAR and DMP are in regular latitude/longitude grid WGS84
  
  filename <- paste("EA/", tools::file_path_sans_ext(basename(files[i])), ".tif", sep = "")
  
  # Crop and mask the raster
  fAPAR_raster.masked <- mask(crop(fAPAR_raster, extent(Maskshp)), Maskshp)
  
  # Export file
  writeRaster(fAPAR_raster.masked, paste(iDir, "/", filename, sep = ""), overwrite = TRUE)
  
}

