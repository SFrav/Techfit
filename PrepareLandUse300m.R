library(raster)
#library(sp)
library(sf)
library(rgdal)

setwd("D:/GIS/Copernicus/Landcover/") #To location of Copernicus 100m landcover layers

###Grass
##Import and merge rasters
LUGrass1 <- raster('E020N00/E020N00_ProbaV_LC100_epoch2015_global_v2.0.1_grass-coverfraction-layer_EPSG-4326.tif')
LUGrass2 <- raster('E020N20/E020N20_ProbaV_LC100_epoch2015_global_v2.0.1_grass-coverfraction-layer_EPSG-4326.tif')
LUGrass3 <- raster('E040N20/E040N20_ProbaV_LC100_epoch2015_global_v2.0.1_grass-coverfraction-layer_EPSG-4326.tif')

LUGrass <- Reduce(function(...)merge(...,tolerance=1),c(LUGrass1,LUGrass2, LUGrass3)) #!E040N00 Add in missing corner or north-east Kenya

##Crop raster to same extent
base300m <- raster("D:/GIS/Copernicus/fAPAR/EA/c_gls_FAPAR300_201501100000_GLOBE_PROBAV_V1.0.1.tif")
LUGrass.cropped <- crop(LUGrass, base300m) #LUGrass.masked <- crop(LUGrass.masked, extent(base300m)) #@Either works

##Aggregate to lower resolution
LUGrass300m <- aggregate(LUGrass.cropped, fact = 3, fun = mean)

##Mask out values beyond 300m res layers
LUGrass300m.masked <- mask(LUGrass300m, base300m)

##Export file
writeRaster(LUGrass.cropped, 'LUgrass100.tif')
writeRaster(LUGrass300m.masked, 'LUgrass300.tif')

rm(list=ls())
.rs.restartR() #Clear disk usage

####Crops
##Import and merge rasters
LUcrops1 <- raster('E020N00/E020N00_ProbaV_LC100_epoch2015_global_v2.0.1_crops-coverfraction-layer_EPSG-4326.tif')
LUcrops2 <- raster('E020N20/E020N20_ProbaV_LC100_epoch2015_global_v2.0.1_crops-coverfraction-layer_EPSG-4326.tif')
LUcrops3 <- raster('E040N20/E040N20_ProbaV_LC100_epoch2015_global_v2.0.1_crops-coverfraction-layer_EPSG-4326.tif')

LUcrops <- Reduce(function(...)merge(...,tolerance=1),c(LUcrops1,LUcrops2, LUcrops3))

##Crop raster to same extent
base300m <- raster("D:/GIS/Copernicus/fAPAR/EA/c_gls_FAPAR300_201501100000_GLOBE_PROBAV_V1.0.1.tif")
LUcrops.cropped <- crop(LUcrops, base300m) #LUcrops.masked <- crop(LUcrops.masked, extent(base300m)) #@Either works


##Aggregate to lower resolution
LUcrops300m <- aggregate(LUcrops.cropped, fact = 3, fun = mean)

##Mask out values beyond 300m res layers
LUcrops300m.masked <- mask(LUcrops300m, base300m)

##Export file
writeRaster(LUcrops.cropped, 'LUcrops100.tif')
writeRaster(LUcrops300m.masked, 'LUcrops300.tif')
rm(list=ls())
.rs.restartR() #Clear disk usage

###Trees
##Import and merge rasters
LUtree1 <- raster('E020N00/E020N00_ProbaV_LC100_epoch2015_global_v2.0.1_tree-coverfraction-layer_EPSG-4326.tif')
LUtree2 <- raster('E020N20/E020N20_ProbaV_LC100_epoch2015_global_v2.0.1_tree-coverfraction-layer_EPSG-4326.tif')
LUtree3 <- raster('E040N20/E040N20_ProbaV_LC100_epoch2015_global_v2.0.1_tree-coverfraction-layer_EPSG-4326.tif')

LUtree <- Reduce(function(...)merge(...,tolerance=1),c(LUtree1,LUtree2, LUtree3))

##Crop raster to same extent
base300m <- raster("D:/GIS/Copernicus/fAPAR/EA/c_gls_FAPAR300_201501100000_GLOBE_PROBAV_V1.0.1.tif")
LUtree.cropped <- crop(LUtree, base300m) #LUtree.masked <- crop(LUtree.masked, extent(base300m)) #@Either works


##Aggregate to lower resolution
LUtree300m <- aggregate(LUtree.cropped, fact = 3, fun = mean)

##Mask out values beyond 300m res layers
LUtree300m.masked <- mask(LUtree300m, base300m)

##Export file
#writeRaster(Maskshp.r, 'aoiMask.tif')
writeRaster(LUtree.cropped, 'LUtree100.tif')
writeRaster(LUtree300m.masked, 'LUtree300.tif')

rm(list=ls())
.rs.restartR() #Clear disk usage

###Forest Type
###Trees
##Import and merge rasters
LUforest1 <- raster('E020N00/E020N00_ProbaV_LC100_epoch2015_global_v2.0.1_forest-type-layer_EPSG-4326.tif')
LUforest2 <- raster('E020N20/E020N20_ProbaV_LC100_epoch2015_global_v2.0.1_forest-type-layer_EPSG-4326.tif')
LUforest3 <- raster('E040N20/E040N20_ProbaV_LC100_epoch2015_global_v2.0.1_forest-type-layer_EPSG-4326.tif')

LUforest <- Reduce(function(...)merge(...,tolerance=1),c(LUforest1,LUforest2, LUforest3))

##Crop raster to same extent
base300m <- raster("D:/GIS/Copernicus/fAPAR/EA/c_gls_FAPAR300_201501100000_GLOBE_PROBAV_V1.0.1.tif")
LUforest.cropped <- crop(LUforest, base300m) #LUforest.masked <- crop(LUforest.masked, extent(base300m)) #@Either works


##Aggregate to lower resolution
LUforest300m <- aggregate(LUforest.cropped, fact = 3, fun = max) # take the maximum value of the 9 pixels. Types are either 2 or 4

##Mask out values beyond 300m res layers
LUforest300m.masked <- mask(LUforest300m, base300m)

##Export file
#writeRaster(Maskshp.r, 'aoiMask.tif')
writeRaster(LUforest.cropped, 'LUforestType100.tif')
writeRaster(LUforest300m.masked, 'LUforestType300.tif')
