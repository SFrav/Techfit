library(raster)
library(gdalUtils)

##@Data is 'days since 1970' -> as.Date("1970-01-01") + 16770

memory.limit(size=64000)

LUgrass <- raster('SpatialData/inputs/Feed_DrySeason/LandUse/LUgrass300.tif')

phenPath <- 'SpatialData/inputs/Feed_DrySeason/PhenologyModis/'
filenames <- list.files(path = phenPath,pattern=".hdf$",full.names = T)

#Convert HDF4 files to tif, extracting both bands. netHDF packages don't work in this instance. gdal_translate doesn't make a perfect copy though.
for (filename in filenames)
{
  sds <- get_subdatasets(filename)
  gdal_translate(sds[grep(pattern = ":Greenup", sds)], b = 1, dst_dataset=paste0(substr(filename, 1, nchar(filename)-4), "_greenup1" ,".tif"))
  
  gdal_translate(sds[grep(pattern = "Peak", sds)], b = 1, dst_dataset=paste0(substr(filename, 1, nchar(filename)-4), "_peak1" ,".tif"))
  
  gdal_translate(sds[grep(pattern = "Senescence", sds)], b = 1, dst_dataset=paste0(substr(filename, 1, nchar(filename)-4), "_senescence1" ,".tif"))
  
  gdal_translate(sds[grep(pattern = "Maturity", sds)], b = 1, dst_dataset=paste0(substr(filename, 1, nchar(filename)-4), "_maturity1" ,".tif"))
  
  gdal_translate(sds[grep(pattern = "Dormancy", sds)], b = 1, dst_dataset=paste0(substr(filename, 1, nchar(filename)-4), "_dormancy1" ,".tif"))
  
  gdal_translate(sds[grep(pattern = ":Greenup", sds)], b = 2, dst_dataset=paste0(substr(filename, 1, nchar(filename)-4), "_greenup2" ,".tif"))
  
  gdal_translate(sds[grep(pattern = "Peak", sds)], b = 2, dst_dataset=paste0(substr(filename, 1, nchar(filename)-4), "_peak2" ,".tif"))
  
  gdal_translate(sds[grep(pattern = "Senescence", sds)], b = 2, dst_dataset=paste0(substr(filename, 1, nchar(filename)-4), "_senescence2" ,".tif"))
  
  gdal_translate(sds[grep(pattern = "Maturity", sds)], b = 2, dst_dataset=paste0(substr(filename, 1, nchar(filename)-4), "_maturity2" ,".tif"))
  
  gdal_translate(sds[grep(pattern = "Dormancy", sds)], b = 2, dst_dataset=paste0(substr(filename, 1, nchar(filename)-4), "_dormancy2" ,".tif"))
  
  gdal_translate(sds[grep(pattern = "NumCycles", sds)], dst_dataset=paste0(substr(filename, 1, nchar(filename)-4), "_numcycles" ,".tif"))
  
}


##Merge, interpolate and crop all rasters
width = 19

filenamesTif <- list.files(path = phenPath ,pattern=".tif$",full.names = T)

jloop <- 8 #How many tiles (HDF files included) - 8 for EA

for(j in 1:jloop){ #Number of layers to merge
  assign(paste0("tmp", j), raster(filenamesTif[grep(pattern = "greenup1", filenamesTif)][j]))
 }
  phenoGreenup1 <- Reduce(function(...)merge(...,tolerance=1), c(tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8))
  rm(list=ls(pattern="tmp"))
  save.image("~/ILRI/FeedSurfaces/FEASTdata/SpatialData/inputs/Feed_DrySeason/PhenologyModis/tmpWorkspace.RData")
  .rs.restartR() 
  memory.limit(size=64000)
  phenoGreenup1 <- projectRaster(phenoGreenup1, crs = crs(LUgrass))
  phenoGreenup1 <- crop(phenoGreenup1, extent(LUgrass))
  phenoGreenup1 <- resample(phenoGreenup1, LUgrass, method = "bilinear") #Resample from 500m to 300m using method for continuous data data
  phenoGreenup1 <- focal(phenoGreenup1, w=matrix(1,nrow=15, ncol=7), fun=mean, NAonly=TRUE, na.rm=TRUE) #!!!Very rough interpolation
  phenoGreenup1 <- focal(phenoGreenup1, w=matrix(1,nrow=width, ncol=width), fun=mean, NAonly=TRUE, na.rm=TRUE)
  phenoGreenup1 <- focal(phenoGreenup1, w=matrix(1,nrow=width, ncol=width), fun=mean, NAonly=TRUE, na.rm=TRUE)
  phenoGreenup1 <- mask(phenoGreenup1, LUgrass)
  writeRaster(phenoGreenup1, paste0(phenPath, "/outputTif/phenoGreenup1.tif"))
  rm(phenoGreenup1)
  gc()
  
    
for(j in 1:jloop){
  assign(paste0("tmp", j), raster(filenamesTif[grep(pattern = "peak1", filenamesTif)][j]))
  }
  phenoPeak1 <- Reduce(function(...)merge(...,tolerance=1),c(tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8))
  rm(list=ls(pattern="tmp"))
  save.image("~/ILRI/FeedSurfaces/FEASTdata/SpatialData/inputs/Feed_DrySeason/PhenologyModis/tmpWorkspace.RData")
  .rs.restartR() 
  memory.limit(size=64000)
  phenoPeak1 <- projectRaster(phenoPeak1, crs = crs(LUgrass))
  phenoPeak1 <- crop(phenoPeak1, extent(LUgrass))
  phenoPeak1 <- resample(phenoPeak1, LUgrass, method = "bilinear") #Resample from 500m to 300m using method for continuous data data
  phenoPeak1 <- focal(phenoPeak1, w=matrix(1,nrow=15, ncol=7), fun=mean, NAonly=TRUE, na.rm=TRUE) #!!!Very rough interpolation
  phenoPeak1 <- focal(phenoPeak1, w=matrix(1,nrow=width, ncol=width), fun=mean, NAonly=TRUE, na.rm=TRUE)
  phenoPeak1 <- focal(phenoPeak1, w=matrix(1,nrow=width, ncol=width), fun=mean, NAonly=TRUE, na.rm=TRUE)
  phenoPeak1 <- mask(phenoPeak1, LUgrass)
  writeRaster(phenoPeak1, paste0(phenPath, "/outputTif/phenoPeak1.tif"))
  rm(phenoPeak1)
  gc()
  save.image("~/ILRI/FeedSurfaces/FEASTdata/SpatialData/inputs/Feed_DrySeason/PhenologyModis/tmpWorkspace.RData")
  .rs.restartR() 
  
for(j in 1:jloop){
  assign(paste0("tmp", j), raster(filenamesTif[grep(pattern = "senescence1", filenamesTif)][j]))
  }
  phenoSenescence1 <- Reduce(function(...)merge(...,tolerance=1),c(tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8))
  rm(list=ls(pattern="tmp"))
  save.image("~/ILRI/FeedSurfaces/FEASTdata/SpatialData/inputs/Feed_DrySeason/PhenologyModis/tmpWorkspace.RData")
  .rs.restartR() 
  memory.limit(size=64000)
  phenoSenescence1 <- projectRaster(phenoSenescence1, crs = crs(LUgrass))
  phenoSenescence1 <- crop(phenoSenescence1, extent(LUgrass))
  phenoSenescence1 <- resample(phenoSenescence1, LUgrass, method = "bilinear") #Resample from 500m to 300m using method for continuous data data
  phenoSenescence1 <- focal(phenoSenescence1, w=matrix(1,nrow=15, ncol=7), fun=mean, NAonly=TRUE, na.rm=TRUE) #!!!Very rough interpolation
  phenoSenescence1 <- focal(phenoSenescence1, w=matrix(1,nrow=width, ncol=width), fun=mean, NAonly=TRUE, na.rm=TRUE)
  phenoSenescence1 <- focal(phenoSenescence1, w=matrix(1,nrow=width, ncol=width), fun=mean, NAonly=TRUE, na.rm=TRUE)
  phenoSenescence1 <- mask(phenoSenescence1, LUgrass)
  writeRaster(phenoSenescence1, paste0(phenPath, "/outputTif/phenoSenescence1.tif"))
  rm(phenoSenescence1)
  gc()
  save.image("~/ILRI/FeedSurfaces/FEASTdata/SpatialData/inputs/Feed_DrySeason/PhenologyModis/tmpWorkspace.RData")
  .rs.restartR() 
  
for(j in 1:jloop){
  assign(paste0("tmp", j), raster(filenamesTif[grep(pattern = "maturity1", filenamesTif)][j]))
  }
  phenoMaturity1 <- Reduce(function(...)merge(...,tolerance=1),c(tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8))
  rm(list=ls(pattern="tmp"))
  save.image("~/ILRI/FeedSurfaces/FEASTdata/SpatialData/inputs/Feed_DrySeason/PhenologyModis/tmpWorkspace.RData")
  .rs.restartR() 
  memory.limit(size=64000)
  phenoMaturity1 <- projectRaster(phenoMaturity1, crs = crs(LUgrass))
  phenoMaturity1 <- resample(phenoMaturity1, LUgrass, method = "bilinear") #Resample from 500m to 300m using method for continuous data data
  phenoMaturity1 <- crop(phenoMaturity1, extent(LUgrass))
  phenoMaturity1 <- focal(phenoMaturity1, w=matrix(1,nrow=width, ncol=width), fun=mean, NAonly=TRUE, na.rm=TRUE) #!!!Very rough interpolation
  phenoMaturity1 <- focal(phenoMaturity1, w=matrix(1,nrow=width, ncol=width), fun=mean, NAonly=TRUE, na.rm=TRUE)
  phenoMaturity1 <- focal(phenoMaturity1, w=matrix(1,nrow=width, ncol=width), fun=mean, NAonly=TRUE, na.rm=TRUE)
  phenoMaturity1 <- mask(phenoMaturity1, LUgrass)
  writeRaster(phenoMaturity1, paste0(phenPath, "/outputTif/phenoMaturity1.tif"))
  rm(phenoMaturity1)
  gc()
  save.image("~/ILRI/FeedSurfaces/FEASTdata/SpatialData/inputs/Feed_DrySeason/PhenologyModis/tmpWorkspace.RData")
  .rs.restartR() 
  
for(j in 1:jloop){
  assign(paste0("tmp", j), raster(filenamesTif[grep(pattern = "dormancy1", filenamesTif)][j]))
  }
  phenoDormancy1 <- Reduce(function(...)merge(...,tolerance=1),c(tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8))
  rm(list=ls(pattern="tmp"))
  save.image("~/ILRI/FeedSurfaces/FEASTdata/SpatialData/inputs/Feed_DrySeason/PhenologyModis/tmpWorkspace.RData")
  .rs.restartR() 
  memory.limit(size=64000)
  phenoDormancy1 <- projectRaster(phenoDormancy1, crs = crs(LUgrass))
  phenoDormancy1 <- crop(phenoDormancy1, extent(LUgrass))
  phenoDormancy1 <- resample(phenoDormancy1, LUgrass, method = "bilinear") #Resample from 500m to 300m using method for continuous data data
  phenoDormancy1 <- focal(phenoDormancy1, w=matrix(1,nrow=15, ncol=7), fun=mean, NAonly=TRUE, na.rm=TRUE) #!!!Very rough interpolation
  phenoDormancy1 <- focal(phenoDormancy1, w=matrix(1,nrow=width, ncol=width), fun=mean, NAonly=TRUE, na.rm=TRUE)
  phenoDormancy1 <- focal(phenoDormancy1, w=matrix(1,nrow=width, ncol=width), fun=mean, NAonly=TRUE, na.rm=TRUE)
  phenoDormancy1 <- mask(phenoDormancy1, LUgrass)
  writeRaster(phenoDormancy1, paste0(phenPath, "/outputTif/phenoDormancy1.tif"))
  rm(phenoDormancy1)
  gc()
  save.image("~/ILRI/FeedSurfaces/FEASTdata/SpatialData/inputs/Feed_DrySeason/PhenologyModis/tmpWorkspace.RData")
  .rs.restartR() 
  
for(j in 1:jloop){
  assign(paste0("tmp", j), raster(filenamesTif[grep(pattern = "numcycles", filenamesTif)][j]))
  }
  phenoNumcycles <- Reduce(function(...)merge(...,tolerance=1),c(tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8))
  rm(list=ls(pattern="tmp"))
  save.image("~/ILRI/FeedSurfaces/FEASTdata/SpatialData/inputs/Feed_DrySeason/PhenologyModis/tmpWorkspace.RData")
  .rs.restartR() 
  memory.limit(size=64000)
  phenoNumcycles <- projectRaster(phenoNumcycles, crs = crs(LUgrass))
  phenoNumcycles <- crop(phenoNumcycles, extent(LUgrass))
  phenoNumcycles <- resample(phenoNumcycles, LUgrass, method = "bilinear") #Resample from 500m to 300m using method for continuous data data
  phenoNumcycles <- focal(phenoNumcycles, w=matrix(1,nrow=15, ncol=7), fun=mean, NAonly=TRUE, na.rm=TRUE) #!!!Very rough interpolation
  phenoNumcycles <- mask(phenoNumcycles, LUgrass)
  writeRaster(phenoNumcycles, paste0(phenPath, "/outputTif/phenoNumcycles.tif"))
  rm(phenoNumcycles)
  gc()
  save.image("~/ILRI/FeedSurfaces/FEASTdata/SpatialData/inputs/Feed_DrySeason/PhenologyModis/tmpWorkspace.RData")
  .rs.restartR() 
  


  for(j in 1:jloop){ #Number of layers to merge
    assign(paste0("tmp", j), raster(filenamesTif[grep(pattern = "greenup2", filenamesTif)][j]))
  }
  phenoGreenup2 <- Reduce(function(...)merge(...,tolerance=1), c(tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8))
  rm(list=ls(pattern="tmp"))
  save.image("~/ILRI/FeedSurfaces/FEASTdata/SpatialData/inputs/Feed_DrySeason/PhenologyModis/tmpWorkspace.RData")
  .rs.restartR() 
  memory.limit(size=64000)
  phenoGreenup2 <- projectRaster(phenoGreenup2, crs = crs(LUgrass))
  phenoGreenup2 <- crop(phenoGreenup2, extent(LUgrass))
  phenoGreenup2 <- resample(phenoGreenup2, LUgrass, method = "bilinear") #Resample from 500m to 300m using method for continuous data data
  phenoGreenup2 <- focal(phenoGreenup2, w=matrix(1,nrow=15, ncol=7), fun=mean, NAonly=TRUE, na.rm=TRUE) #!!!Very rough interpolation
  phenoGreenup2 <- focal(phenoGreenup2, w=matrix(1,nrow=width, ncol=width), fun=mean, NAonly=TRUE, na.rm=TRUE)
  phenoGreenup2 <- focal(phenoGreenup2, w=matrix(1,nrow=width, ncol=width), fun=mean, NAonly=TRUE, na.rm=TRUE)
  phenoGreenup2 <- mask(phenoGreenup2, LUgrass)
  writeRaster(phenoGreenup2, paste0(phenPath, "/outputTif/phenoGreenup2.tif"))
  rm(phenoGreenup2)
  gc()
  save.image("~/ILRI/FeedSurfaces/FEASTdata/SpatialData/inputs/Feed_DrySeason/PhenologyModis/tmpWorkspace.RData")
  .rs.restartR() 
  
  
  
  for(j in 1:jloop){
    assign(paste0("tmp", j), raster(filenamesTif[grep(pattern = "peak2", filenamesTif)][j]))
  }
  phenopeak2 <- Reduce(function(...)merge(...,tolerance=1),c(tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8))
  rm(list=ls(pattern="tmp"))
  save.image("~/ILRI/FeedSurfaces/FEASTdata/SpatialData/inputs/Feed_DrySeason/PhenologyModis/tmpWorkspace.RData")
  .rs.restartR() 
  memory.limit(size=64000)
  phenopeak2 <- projectRaster(phenopeak2, crs = crs(LUgrass))
  phenopeak2 <- crop(phenopeak2, extent(LUgrass))
  phenopeak2 <- resample(phenopeak2, LUgrass, method = "bilinear") #Resample from 500m to 300m using method for continuous data data
  phenopeak2 <- focal(phenopeak2, w=matrix(1,nrow=15, ncol=7), fun=mean, NAonly=TRUE, na.rm=TRUE) #!!!Very rough interpolation
  phenopeak2 <- focal(phenopeak2, w=matrix(1,nrow=width, ncol=width), fun=mean, NAonly=TRUE, na.rm=TRUE)
  phenopeak2 <- focal(phenopeak2, w=matrix(1,nrow=width, ncol=width), fun=mean, NAonly=TRUE, na.rm=TRUE)
  phenopeak2 <- mask(phenopeak2, LUgrass)
  writeRaster(phenopeak2, paste0(phenPath, "/outputTif/phenoPeak2.tif"))
  rm(phenopeak2)
  gc()
  save.image("~/ILRI/FeedSurfaces/FEASTdata/SpatialData/inputs/Feed_DrySeason/PhenologyModis/tmpWorkspace.RData")
  .rs.restartR() 
  
  
  for(j in 1:jloop){
    assign(paste0("tmp", j), raster(filenamesTif[grep(pattern = "senescence2", filenamesTif)][j]))
  }
  phenoSenescence2 <- Reduce(function(...)merge(...,tolerance=1),c(tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8))
  rm(list=ls(pattern="tmp"))
  save.image("~/ILRI/FeedSurfaces/FEASTdata/SpatialData/inputs/Feed_DrySeason/PhenologyModis/tmpWorkspace.RData")
  .rs.restartR() 
  memory.limit(size=64000)
  phenoSenescence2 <- projectRaster(phenoSenescence2, crs = crs(LUgrass))
  phenoSenescence2 <- crop(phenoSenescence2, extent(LUgrass))
  phenoSenescence2 <- resample(phenoSenescence2, LUgrass, method = "bilinear") #Resample from 500m to 300m using method for continuous data data
  phenoSenescence2 <- focal(phenoSenescence2, w=matrix(1,nrow=15, ncol=7), fun=mean, NAonly=TRUE, na.rm=TRUE) #!!!Very rough interpolation
  phenoSenescence2 <- focal(phenoSenescence2, w=matrix(1,nrow=width, ncol=width), fun=mean, NAonly=TRUE, na.rm=TRUE)
  phenoSenescence2 <- focal(phenoSenescence2, w=matrix(1,nrow=width, ncol=width), fun=mean, NAonly=TRUE, na.rm=TRUE)
  phenoSenescence2 <- mask(phenoSenescence2, LUgrass)
  writeRaster(phenoSenescence2, paste0(phenPath, "/outputTif/phenoSenescence2.tif"))
  rm(phenoSenescence2)
  gc()
  save.image("~/ILRI/FeedSurfaces/FEASTdata/SpatialData/inputs/Feed_DrySeason/PhenologyModis/tmpWorkspace.RData")
  .rs.restartR() 
  
  
  for(j in 1:jloop){
    assign(paste0("tmp", j), raster(filenamesTif[grep(pattern = "maturity2", filenamesTif)][j]))
  }
  phenoMaturity2 <- Reduce(function(...)merge(...,tolerance=1),c(tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8))
  rm(list=ls(pattern="tmp"))
  save.image("~/ILRI/FeedSurfaces/FEASTdata/SpatialData/inputs/Feed_DrySeason/PhenologyModis/tmpWorkspace.RData")
  .rs.restartR() 
  memory.limit(size=64000)
  phenoMaturity2 <- projectRaster(phenoMaturity2, crs = crs(LUgrass))
  phenoMaturity2 <- resample(phenoMaturity2, LUgrass, method = "bilinear") #Resample from 500m to 300m using method for continuous data data
  phenoMaturity2 <- crop(phenoMaturity2, extent(LUgrass))
  phenoMaturity2 <- focal(phenoMaturity2, w=matrix(1,nrow=width, ncol=width), fun=mean, NAonly=TRUE, na.rm=TRUE) #!!!Very rough interpolation
  phenoMaturity2 <- focal(phenoMaturity2, w=matrix(1,nrow=width, ncol=width), fun=mean, NAonly=TRUE, na.rm=TRUE)
  phenoMaturity2 <- focal(phenoMaturity2, w=matrix(1,nrow=width, ncol=width), fun=mean, NAonly=TRUE, na.rm=TRUE)
  phenoMaturity2 <- mask(phenoMaturity2, LUgrass)
  writeRaster(phenoMaturity2, paste0(phenPath, "/outputTif/phenoMaturity2.tif"))
  rm(phenoMaturity2)
  gc()
  save.image("~/ILRI/FeedSurfaces/FEASTdata/SpatialData/inputs/Feed_DrySeason/PhenologyModis/tmpWorkspace.RData")
  .rs.restartR() 
  
  
  for(j in 1:jloop){
    assign(paste0("tmp", j), raster(filenamesTif[grep(pattern = "dormancy2", filenamesTif)][j]))
  }
  phenoDormancy2 <- Reduce(function(...)merge(...,tolerance=1),c(tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8))
  rm(list=ls(pattern="tmp"))
  save.image("~/ILRI/FeedSurfaces/FEASTdata/SpatialData/inputs/Feed_DrySeason/PhenologyModis/tmpWorkspace.RData")
  .rs.restartR() 
  memory.limit(size=64000)
  phenoDormancy2 <- projectRaster(phenoDormancy2, crs = crs(LUgrass))
  phenoDormancy2 <- crop(phenoDormancy2, extent(LUgrass))
  phenoDormancy2 <- resample(phenoDormancy2, LUgrass, method = "bilinear") #Resample from 500m to 300m using method for continuous data data
  phenoDormancy2 <- focal(phenoDormancy2, w=matrix(1,nrow=15, ncol=7), fun=mean, NAonly=TRUE, na.rm=TRUE) #!!!Very rough interpolation
  phenoDormancy2 <- focal(phenoDormancy2, w=matrix(1,nrow=width, ncol=width), fun=mean, NAonly=TRUE, na.rm=TRUE)
  phenoDormancy2 <- focal(phenoDormancy2, w=matrix(1,nrow=width, ncol=width), fun=mean, NAonly=TRUE, na.rm=TRUE)
  phenoDormancy2 <- mask(phenoDormancy2, LUgrass)
  writeRaster(phenoDormancy2, paste0(phenPath, "/outputTif/phenoDormancy2.tif"))
  rm(phenoDormancy2)
  gc()
  save.image("~/ILRI/FeedSurfaces/FEASTdata/SpatialData/inputs/Feed_DrySeason/PhenologyModis/tmpWorkspace.RData")
  .rs.restartR() 
  










