#############################################################          
# Techfit feasibility surface            .=     ,        =.
#                               _  _   /'/    )\,/,/(_   \ \
#                            `//-.|  (  ,\\)\//\)\/_  ) |
#                            //___\   `\\\/\\/\/\\///'  /
# Simon Fraval             ,-"~`-._ `"--'_   `"""`  _ \`'"~-,_
# V1: 23/9/2019           \       `-.  '_`.      .'_` \ ,-"~`/
#                          `.__.-'`/   (-\        /-) |-.__,'
#                            ||   |     \O)  /^\ (O/  |
#                           `\\  |         /   `\    /
#                             \\  \       /      `\ /
#                              `\\ `-.  /' .---.--.\
#                                `\\/`~(, '()      ('
#                                 /(O) \\   _,.-.,_)
#                                //  \\ `\'`      /
#                                / |  ||   `""""~"`
#                              /'  |__||
#                                    `o 
###########################################################

library(raster)
library(truncnorm)
library(sf)

library(rgdal)
library(dplyr)
library(ggplot2)
library(tmap) #interactive map

##Techfit data
techfit <- read.csv('TechfitOriginal.csv')
interventions <- techfit$Intervention
techfit <- techfit[-1,] #remove 'test' entry
row.names(techfit) <- 1:nrow(techfit)

###Spatial data
##Admin boundaries
Ethiopia <- st_read("C:/Users/Simon/Documents/GIS DataBase/Africa/AdminBoundaries/ETH_shp/gadm36_ETH_3.shp", layer ="gadm36_ETH_3")
Kenya <- st_read("C:/Users/Simon/Documents/GIS DataBase/Africa/AdminBoundaries/KEN_shp/gadm36_KEN_3.shp", layer ="gadm36_KEN_3")
Tanzania <- st_read("C:/Users/Simon/Documents/GIS DataBase/Africa/AdminBoundaries/TZA_shp/gadm36_TZA_3.shp", layer ="gadm36_TZA_3")
Uganda <- st_read("C:/Users/Simon/Documents/GIS DataBase/Africa/AdminBoundaries/UGA_shp/gadm36_UGA_3.shp", layer ="gadm36_UGA_3")
Rwanda <- st_read("C:/Users/Simon/Documents/GIS DataBase/Africa/AdminBoundaries/RWA_shp/gadm36_RWA_3.shp", layer ="gadm36_RWA_3")
Burundi <- st_read("C:/Users/Simon/Documents/GIS DataBase/Africa/AdminBoundaries/BDI_shp/gadm36_BDI_3.shp", layer ="gadm36_BDI_3")

EA_adm3 <- rbind(Ethiopia, Kenya, Tanzania, Uganda, Rwanda, Burundi)
rm(Ethiopia, Kenya, Tanzania, Uganda, Rwanda, Burundi)

##Add production system layer
systemMixedCrop <- st_read("SpatialData/inputs/Livestock_Production_Systems/lpsVector/lpsVector.shp", layer = "lpsVector") # raster("SpatialData/Land_Cover/EA_Pasture_Cropland_Area.tif") # use to subset EA_adm3

EA_adm3_cropLvst <- EA_adm3[systemMixedCrop, op = st_intersects]
EA_adm3_cropLvst <- sf:::as_Spatial(EA_adm3_cropLvst)


##Add constraint layers
constraintQuantity <- raster('SpatialData/inputs/Feed_Scarcity/FS_per_TLU.tif')

#constraintSeasonalityDryCool ; constraintSeasonalityGrowing ; constraintQuantityConstraint ; constraintQualityConstraint

##Add species layers
speciesCattleDairy <- raster("SpatialData/inputs/Livestock_Density/AF_Cattle1km_AD_2010_v2_1.tif")

#speciesCattleBreeding ; speciesCattleFattening ; sheepBreeding ; sheepFattening ; pigBreeding ; pigFattening ; 


##Add attributes layers
attributeLand <- raster("SpatialData/inputs/Land_Availability/crop_land_per_person.tif")
attributeWater <- raster("SpatialData/inputs/Water_Availability/Annual_Average_Precipitation/annual_average_precipitation_2008_2017.tif")

attributeCash_Input_Skill <- raster("C:/Users/Simon/Documents/GIS DataBase/TravelTime/2015_accessibility_to_cities_v1.0/2015_accessibility_to_cities_v1.0.tif")

#attributeCash ; attributeInput ; attributeSkill ; attributeMarket


##Crop and stack layers
speciesCattleDairy <- crop(speciesCattleDairy, extent(constraintQuantity))
speciesCattleDairy <- mask(speciesCattleDairy, constraintQuantity)
attributeCash_Input_Skill <- crop(attributeCash_Input_Skill, extent(constraintQuantity))
attributeCash_Input_Skill <- mask(attributeCash_Input_Skill, constraintQuantity)
feedConditions <- stack(constraintQuantity, speciesCattleDairy, attributeLand, attributeWater, attributeCash_Input_Skill)

feedConditionsDat <- extract(feedConditions, fun = 'mean', EA_adm3_cropLvst, df = T)
feedConditionsDatSD <- extract(feedConditions, fun = 'sd', EA_adm3_cropLvst, df = T)
colnames(feedConditionsDatSD)[2:ncol(feedConditionsDatSD)] <- paste(colnames(feedConditionsDatSD)[2:ncol(feedConditionsDatSD)], "sd", sep = "_")

EA_adm3_cropLvst@data <- cbind(EA_adm3_cropLvst@data, feedConditionsDat)
EA_adm3_cropLvst@data <- cbind(EA_adm3_cropLvst@data, feedConditionsDatSD)

EA_adm3_cropLvst2 <- st_as_sf(EA_adm3_cropLvst)
st_write(EA_adm3_cropLvst2, "SpatialData/outputs/EA_adm3_feedConditionsMeanSD2.shp")

EA_adm3_cropLvst2 <- st_read("SpatialData/outputs/EA_adm3_feedConditionsMeanSD2.shp", layer ="EA_adm3_feedConditionsMeanSD2")


###techfit_match
##Constraints
EA_adm3_cropLvst2$FS__TLU <- rnorm(length(EA_adm3_cropLvst2$FS__TLU), 1, 0.4)
EA_adm3_cropLvst2$feedQuantityBin <- ifelse(EA_adm3_cropLvst2$FS__TLU > 0.75, 1, 0) #greater than 75% of livestock needs then match to mitigation
EA_adm3_cropLvst2$feedQuantityBin[EA_adm3_cropLvst2$feedQuantityBin] <- 0

##Species bin
EA_adm3_cropLvst2$cattleBin <- ifelse(EA_adm3_cropLvst2$AF_C1_A > stats::quantile(EA_adm3_cropLvst2$AF_C1_A, probs = 0.05, na.rm=T), 1, 0) #greater than 5th percentile
EA_adm3_cropLvst2$cattleBin[is.na(EA_adm3_cropLvst2$cattleBin)] <- 0

##Attributes
EA_adm3_cropLvst2$suitabLand <- ifelse(EA_adm3_cropLvst2$crp_l__ < 0.2 & EA_adm3_cropLvst2$crp____ < 1 | is.na(EA_adm3_cropLvst2$crp_l__), 4, 
                                    ifelse(EA_adm3_cropLvst2$crp_l__ < 0.2 & EA_adm3_cropLvst2$crp____ > 1, 3,
                                    ifelse(EA_adm3_cropLvst2$crp_l__ > 0.2 & EA_adm3_cropLvst2$crp_l__ <= 1, 3,
                                    ifelse(EA_adm3_cropLvst2$crp_l__ > 1 & EA_adm3_cropLvst2$crp_l__ <= 3, 3,
                                    ifelse(EA_adm3_cropLvst2$crp_l__ >3 & EA_adm3_cropLvst2$crp_l__ <= 10, 2,
                                    ifelse(EA_adm3_cropLvst2$crp_l__ > 10, 1, NA)))))) #if less than set threshold and sd is less than 1
#!if SD high then show raster in plot
#! remove NA attribution to lowest category


EA_adm3_cropLvst2$suitabWater <- ifelse(EA_adm3_cropLvst2$a___200 < stats::quantile(EA_adm3_cropLvst2$a___200, probs = 0.05, na.rm=T) | is.na(EA_adm3_cropLvst2$a___200), 4, 
                                     ifelse(EA_adm3_cropLvst2$a___200 > stats::quantile(EA_adm3_cropLvst2$a___200, probs = 0.05, na.rm=T) & EA_adm3_cropLvst2$a___200 <= 900, 3,
                                            ifelse(EA_adm3_cropLvst2$a___200 > 900 & EA_adm3_cropLvst2$a___200 <=  1200, 2,
                                                   ifelse(EA_adm3_cropLvst2$a___200 > 1200, 1, NA)))) #greater than 5th percentile. Variability < 100mm, so no need to include


#@Inverse of others. Less is better
EA_adm3_cropLvst2$suitabInput <- ifelse(EA_adm3_cropLvst2$X2015__ < stats::quantile(EA_adm3_cropLvst2$X2015__, probs = 0.05, na.rm=T) | is.na(EA_adm3_cropLvst2$X2015__), 1, 
                                     ifelse(EA_adm3_cropLvst2$X2015__ >= stats::quantile(EA_adm3_cropLvst2$X2015__, probs = 0.05, na.rm=T) & EA_adm3_cropLvst2$X2015__ <= stats::quantile(EA_adm3_cropLvst2$X2015__, probs = 0.25, na.rm=T), 2,
                                            ifelse(EA_adm3_cropLvst2$X2015__ > stats::quantile(EA_adm3_cropLvst2$X2015__, probs = 0.25, na.rm=T) & EA_adm3_cropLvst2$X2015__ <=  stats::quantile(EA_adm3_cropLvst2$X2015__, probs = 0.5, na.rm=T), 3,
                                                   ifelse(EA_adm3_cropLvst2$X2015__ > stats::quantile(EA_adm3_cropLvst2$X2015__, probs = 0.5, na.rm=T), 4, NA)))) #greater than 5th percentile. Variability < 100mm, so no need to include


EA_adm3_cropLvst2$mitFeedQuantFeasiList <- NA
EA_adm3_cropLvst2$mitFeedQuantMax <- NA
for(i in 1:nrow(EA_adm3_cropLvst2)){
  if(EA_adm3_cropLvst2$feedQuantityBin[i] == 0 |  EA_adm3_cropLvst2$cattleBin[i] == 0) {EA_adm3_cropLvst2$mitFeedQuantFeasiList[i] <- NA} else {
      #tmp <- techfit[techfit$Cattle_buffalo.breeding > 0 | techfit$Cattle_buffalo.fattening | Dairy.cattle_buffalo > 0,) #all are above 0 for at least 1 cattle enterprise
      EA_adm3_cropLvst2$mitFeedQuantFeasiList[i] <- paste(row.names(techfit[techfit$Quantity.constraint > 0 & 
                                                    techfit$Attribute.1_Requirement.for.land >= EA_adm3_cropLvst2$suitabLand[i] & 
                                                    techfit$Attribute.2_Requirement.for.water >= EA_adm3_cropLvst2$suitabWater[i] &
                                                    techfit$Attribute.5 >= EA_adm3_cropLvst2$suitabInput[i],]), collapse="_")
      
      EA_adm3_cropLvst2$mitFeedQuantMax[i] <- max(techfit$Quantity.constraint[techfit$Quantity.constraint > 0 & 
                                                                       techfit$Attribute.1_Requirement.for.land >= EA_adm3_cropLvst2$suitabLand[i] & 
                                                                       techfit$Attribute.2_Requirement.for.water >= EA_adm3_cropLvst2$suitabWater[i] &
                                                                       techfit$Attribute.5 >= EA_adm3_cropLvst2$suitabInput[i]], na.rm=T)
  }
  
}

EA_adm3_cropLvst2$mitFeedQuantFeasiList[!(sapply(EA_adm3_cropLvst2$mitFeedQuantFeasiList, length))] <- NA
EA_adm3_cropLvst2$mitFeedQuantMax[is.infinite(EA_adm3_cropLvst2$mitFeedQuantMax)] <- NA


#By feed tech
for(i in 1:nrow(techfit)){
  EA_adm3_cropLvst2 <- cbind(EA_adm3_cropLvst2,
  ifelse(EA_adm3_cropLvst2$feedQuantityBin == 1 & EA_adm3_cropLvst2$cattleBin == 1 & EA_adm3_cropLvst2$suitabLand >= techfit$Attribute.1_Requirement.for.land[i] & EA_adm3_cropLvst2$suitabWater >= techfit$Attribute.2_Requirement.for.water[i] & EA_adm3_cropLvst2$suitabInput >= techfit$Attribute.5_Requirement.for.input.delivery[i],  EA_adm3_cropLvst2$FS__TLU, NA)
  )
  colnames(EA_adm3_cropLvst2)[length(EA_adm3_cropLvst2)-1] <- paste0("tech", i) #Geom always last, so -1
}


#techfitCattleMixedSys <- techfit[techfit$Cattle_buffalo.breeding > 0 | techfit$Cattle_buffalo.fattening >0 | techfit$Dairy.cattle_buffalo >0 & techfit$Intensive.mixed.croplivestock.system >0 ,]


##Intervention lists - interactive
EA_adm3_cropLvst2sub <- st_crop(EA_adm3_cropLvst2, c(xmin= 35, ymax = 1, xmax = 38, ymin = -1)) #@limit extent for quicker plotting
EA_adm3_cropLvst2sub <- select(EA_adm3_cropLvst2sub, -(GID_0:GID_3), -(VARNAME:ID))
EA_adm3_cropLvst2sub$Mitigation_potential <- as.factor(EA_adm3_cropLvst2sub$mitFeedQuantMax)
tmap_mode("view")
tm_shape(EA_adm3_cropLvst2sub) + tm_polygons("Mitigation_potential", id = "mitFeedQuantFeasiList", popup.vars = c("NAME_3","mitFeedQuantFeasiList", "suitabWater", "suitabLand", "suitabInput"))



###########
##Exploratory mapping
##Suitability maps
ggplot() + geom_sf(data = EA_adm3_cropLvst2sub, aes(fill = tech17)) + ggtitle(techfit$Intervention[17]) + guides(fill = F) # + labs(fill = "DM adequacy per TLU")

##Maximum  mitigation score
ggplot() + geom_sf(data = EA_adm3_cropLvst2sub, aes(fill = mitFeedQuantMax)) + ggtitle("Maximum mitigation potential")  + guides(fill = F)


constraintQuantity_df <- as(feedConditions$FS_per_TLU, "SpatialPixelsDataFrame")
constraintQuantity_df <- as.data.frame(constraintQuantity_df)
colnames(constraintQuantity_df) <- c("value", "x", "y")

DMP_TLU <- ggplot() + 
  geom_sf(data = EA_adm3_cropLvst2, aes(fill = FS_per_TLU), inherit.aes = FALSE) +
  geom_tile(data = constraintQuantity_df, aes(x = x, y = y, fill=value), alpha=0.8)

DMP_TLU <- ggplot() + 
  geom_sf(data = EA_adm3_cropLvst2, aes(fill = FS_per_TLU), inherit.aes = FALSE)

DMP_TLUsd <- ggplot() + 
  geom_sf(data = EA_adm3_cropLvst2, aes(fill = FS_per_TLUsd), inherit.aes = FALSE)


attributeLand_df <- as(attributeLand, "SpatialPixelsDataFrame")
attributeLand_df <- as.data.frame(attributeLand_df)
colnames(attributeLand_df) <- c("value", "x", "y")

DMP_TLU <- ggplot() + 
  geom_sf(data = EA_adm3_cropLvst2, aes(fill = FS_per_TLU), inherit.aes = FALSE) +
  geom_tile(data = attributeLand_df, aes(x = x, y = y, fill=value), alpha=0.8)

