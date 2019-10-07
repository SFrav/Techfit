##Simulate dry season feed scarcity
#
#               ,,
#    ,,       ,\\//,          Detecting feed scarcity in the dry season
#  ,\\//,    ,\\\///,   ,,    Tested R versions: 3.4.1
# ,\\\///,   \\\\//// ,\\//,  Authors: Simon Fraval, John Mutua
# \\\\////    \\\/// ,\\\///, 
#  \\\///     ###### \\\\//// 
#  ######    ////\\\\ \\\///
#  ////\\\\  /////\\\\\######
# /////\\\\\//////\\\\////\\\\
#//////\\\\\\/,///\\\/////\\\\\
#//////_,_\\\\(_)    //////\\\\\\,
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
#ASCI credit: 'Pumpkins in the cornfield', anonymous, 1997
###########################################################
library(raster)
library(sf)
library(rgdal)
library(dplyr)
library(ggplot2)
library(ggspatial) #For mapping scalebar and north arrow


AOI3 <- st_read('SpatialData/inputs/AOI3.shp', layer = "AOI3")
AOI <- AOI3[AOI3$NAME_1 == "Nandi",]
AOI <- as_Spatial(AOI)

DMPorigin_jul20 <- raster('SpatialData/inputs/Feed_DrySeason/DMP/c_gls_DMP300-RT5_201506200000_GLOBE_PROBAV_V1.0.1.tif')
LUgrass <- raster('SpatialData/inputs/Feed_DrySeason/LandUse/LUgrass300.tif')
LUcrop <- raster('SpatialData/inputs/Feed_DrySeason/LandUse/LUcrops300.tif')
LUtrees <- raster('SpatialData/inputs/Feed_DrySeason/LandUse/LUtree300.tif')
LUforestType <- raster('SpatialData/inputs/Feed_DrySeason/LandUse/LUforestType300.tif')


DMPorigin_jul20 <- crop(DMPorigin_jul20, extent(AOI))
DMPorigin_jul20 <- mask(DMPorigin_jul20, AOI)
LUstack <- stack(LUgrass, LUcrop, LUtrees) #LUforestType
LUstack <- crop(LUstack, extent(AOI))
LUstack <- mask(LUstack, AOI)

feedFrac <- function(x,y){(x/100)*(y*3)} #*3 to make it total dry matter available per pixel - in kg
feedStack <- overlay(LUstack, DMPorigin_jul20, fun = feedFrac)

##Create annual DMP
growthDist1 <- c(seq(from = 0.2, to = 1, length.out = 90), #First growing season
                 seq(from = 1, to = 0.2, length.out = 10), #First harvest
                 seq(from = 0.2, to = 0.4, length.out = 60), #Short dry period
                 seq(from = 0.4, to = 1, length.out = 90), #Second growing period
                 seq(from = 1, to = 0.2, length.out = 10), #Second harvest
                 seq(from = 0.2, to = 0.1, length.out = 105)) #Long dry period

#growthDist2 <- c(seq(from = 0.2, to = 1, length.out = 90), #First growing season
#                 seq(from = 1, to = 0.2, length.out = 10), #First harvest
#                 seq(from = 0.2, to = 0.4, length.out = 60), #Short dry period
#                 seq(from = 0.4, to = 1, length.out = 90), #Second growing period
#                 seq(from = 1, to = 0.2, length.out = 10), #Second harvest
#                 seq(from = 0.2, to = 0.5, length.out = 5),
#                 seq(from = 0.5, to = 0.8, length.out = 100)) #Long dry period
#growthDist2 #Apply to latitude
#growthDist3

growingPeriod <- c(1:90, 161:250) #as.Date("1970-01-01") + 16770 #For MODIS data
nocropping <- c(91:160, 251:365)

DMPdailyGrassStack <- (overlay(LUstack$LUgrass300, DMPorigin_jul20, fun = feedFrac)+overlay(LUstack$LUtree300, DMPorigin_jul20, fun = feedFrac)*0.1)*growthDist1[1]
DMPdailyCropStack <- overlay(LUstack$LUcrops300, DMPorigin_jul20, fun = feedFrac)*growthDist1[1]
for(i in 2:length(growthDist1)){
  if(i %in% growingPeriod){
    #Assume 10% of forest DM can be extracted as feed 
  DMPdailyGrassStack <- stack(DMPdailyGrassStack, (overlay(LUstack$LUgrass300, DMPorigin_jul20, fun = feedFrac)+overlay(LUstack$LUtree300, DMPorigin_jul20, fun = feedFrac)*0.1)*growthDist1[i]) 
  }
  #Assume land not used for cropping is used for animal feed. Only 50% paletable + forest feed + grass
  else{DMPdailyGrassStack <- stack(DMPdailyGrassStack, (overlay(LUstack$LUgrass300, DMPorigin_jul20, fun = feedFrac)+overlay(LUstack$LUtree300, DMPorigin_jul20, fun = feedFrac)*0.1 + overlay(LUstack$LUcrops300, DMPorigin_jul20, fun = feedFrac)*0.5)*growthDist1[i])}
  #assign(paste0("DMPgrass", i), (overlay(LUstack$LUgrass300, DMPorigin_jul20, fun = feedFrac)+overlay(LUstack$LUtree300, DMPorigin_jul20, fun = feedFrac)*0.1)*growthDist1[i])
  #DMPdailyGrassStack <- stack(DMPdailyGrassStack, get(paste0("DMPgrass", i)))
  
  if(i %in% growingPeriod){
  DMPdailyCropStack <- stack(DMPdailyGrassStack, overlay(LUstack$LUcrops300, DMPorigin_jul20, fun = feedFrac)*growthDist1[i])
  }
}

residueFrac <- 0.4
utilisationFrac <- calc(LUcrop, fun = function(x){rnorm(x, 0.5, 0.05)})
utilisationFrac <- crop(utilisationFrac, extent(AOI))
utilisationFrac <- mask(utilisationFrac, AOI)

DMPdailyCropStack <- DMPdailyCropStack * residueFrac * utilisationFrac

rm(list = setdiff(ls(), c("DMPdailyCropStack", "DMPdailyGrassStack", "AOI", "growthDist1", "growingPeriod", "nocropping")))

DMPgrowing1mean <- calc(DMPdailyGrassStack[[growingPeriod[1:90]]], fun = mean)
DMPgrowing2mean <- calc(DMPdailyGrassStack[[growingPeriod[91:length(growingPeriod)]]], fun = mean)
DMPdry1mean <- calc(DMPdailyGrassStack[[91:160]], fun = mean)
DMPdry2mean <- calc(DMPdailyGrassStack[[251:365]], fun = mean)

DMPgrowing1sum <- calc(DMPdailyGrassStack[[growingPeriod[1:90]]], fun = sum)
DMPgrowing2sum <- calc(DMPdailyGrassStack[[growingPeriod[91:length(growingPeriod)]]], fun = sum)
DMPdry1sum <- calc(DMPdailyGrassStack[[91:160]], fun = sum)
DMPdry2sum <- calc(DMPdailyGrassStack[[251:365]], fun = sum)

DMPdiff1mean <- DMPdry1mean - DMPgrowing1mean
DMPdiff2mean <- DMPdry2mean - DMPgrowing2mean

DMPdiff1sum <- DMPdry1sum - DMPgrowing1sum
DMPdiff2sum <- DMPdry2sum - DMPgrowing2sum

DMPdryseason1perc <- 1 - abs(DMPdry1sum - DMPgrowing1sum)/DMPgrowing1sum
DMPdryseason2perc <- 1 - abs(DMPdry2sum - DMPgrowing2sum)/DMPgrowing1sum

DMPdryScarcityStack <- stack(DMPdiff1mean, DMPdiff2mean, DMPdiff1sum, DMPdiff2sum, DMPdryseason1perc, DMPdryseason2perc)

### Admin level
feedScarcityDat <- extract(DMPdryScarcityStack, fun = 'mean', na.rm=T, AOI, df = T)
feedScarcityDatSD <- extract(DMPdryScarcityStack, fun = 'sd', na.rm=T, AOI, df = T)
colnames(feedScarcityDatSD)[2:ncol(feedScarcityDatSD)] <- paste(colnames(feedScarcityDatSD)[2:ncol(feedScarcityDatSD)], "sd", sep = "_")

AOI@data <- cbind(AOI@data, feedScarcityDat)
AOI@data <- cbind(AOI@data, feedScarcityDatSD)

AOI <- st_as_sf(AOI)

#plot(AOI[23])

#plot(AOI[32])
AOI$layer.5 <- AOI$layer.5 *100
AOI$layer.6 <- AOI$layer.6 *100

AOI$DSscarcity <- ifelse(AOI$layer.5 < 60 | AOI$layer.6 < 60, "Yes", "No")

AOIplot <- AOI %>% select(layer.5, layer.6) %>% gather(VAR, shortagePerc, -geometry)