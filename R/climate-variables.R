# READ ME ----
#This code computes the climate variables (mean annual temperatures, annual temperature 
#variation etc) that will be used to model seagrass response to heat stress in experiments. 
#The output is an excel csv file with climate variables calculated for each observation (see Methods)

# PREP ----
## LOAD PACKAGES ----
library(readxl)
library(tidyverse)
library(raster)
library(ncdf4)
library(sf)
library(sp)
library(IDPmisc)
#library(maptools)


## FUNCTIONS ----
source('R/functions.R')

## LOAD DATA ----
seagrass_studies <- readxl::read_excel('data/primary/seagrass-thermal-tolerance-studies.xlsx', 
                                       sheet = 2)
Edata <- seagrass_studies[, c("Authors", "Title of the publication", "Year of publication", 
                                        "Year of sample collection","Sample collection site", 
                                        "Collection bioregion", "Latitude of collection site",
                                        "Longitude of collection site", "Family", "Genus","Species", 
                                        "Type","Heat stress exposure duration (hours)", 
                                        "Heat stress temperature (ºC)","Accepted response variable", 
                                        "Survival")]
colnames(Edata) <- c("Study", "Title", "Publication", "Date", "Location", "Bioregion", 
                                "Latitude", "Longitude", "Family", "Genus", "Species",
                                "Type", "Time", "Temperature", "Variable", "Survival")
temp <- raster::stack("data/primary/sst.mon.mean.nc") %>% # load
  .[[((1960-1850)*12+1):nlayers(.)]] %>% # select period 1960 - present
  raster::rotate() #adjust longitude to -180,180
seagrasses <- st_read('data/primary/SEAGRASSES/SEAGRASSES.shp') #seagrass species distributions

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# because of the land sea mask in the temperature raster, some coastal cells where seagrass are present
# do not have temperature values
# fix this by assigning missing pixels as the mean of the surrounding 8 cells (i.e. queens move)

fill.na <- function(x, i=5) {
  if( is.na(x)[i] ) {
    return(mean(x, na.rm=TRUE))
  } else {
    return(x[i])
  }
} 

# https://gis.stackexchange.com/questions/181011/fill-the-gaps-using-nearest-neighbors

# check the function works i.e. fills in NA values but doesn't alter numeric values
rmat <- matrix(rnorm(64,0,1), nrow=8, ncol=8)
r <- raster(rmat)
r[1:2,1:2] <- NA
r2 <- focal(r, w = matrix(1,3,3), fun = fill.na,pad = TRUE, na.rm = FALSE)
plot(values(r),values(r2))
values(r)-values(r2)
mean(r[3,1:2])
values(r2)[which(is.na(values(r))==TRUE)]



## DATA PREP ----
# because of the land sea mask in the temperature raster, some coastal cells where seagrass are present
# do not have temperature values
# fix this by assigning missing pixels as the mean of the surrounding 8 cells using the function fill.na, available in functions.R
tempFilled <- temp
for(i in 1:nlayers(temp)) {
  tempFilled[[i]] <- focal(temp[[i]], w = matrix(1,3,3), fun = fill.na,pad = TRUE, na.rm = FALSE)
  print(i)
}
names(tempFilled) <- names(temp)
temp <- tempFilled

writeRaster(temp, "data/processed/tempFilled.nc")

#create new rasters to store the climate variables for each year
mtwa <- mat <-av <- raster(nrow = nrow(temp), ncol = ncol(temp), 
                           xmn = xmin(temp), xmx = xmax(temp), 
                           ymn = ymin(temp), ymx = ymax(temp),
                           crs = crs(temp))

year<-1960:2019
for (a in 1:length(year)) {
  focLayers <- grep(as.integer(year[a]), names(temp)) #find all layers from a specific year
  sub <- temp[[focLayers]] #subset the stack to only that year
  int_av <- max(sub) - min(sub) #calculate AV for that year
  av <- addLayer(av, int_av) #add AV for that year to a stack
  int_mat <- mean(sub) #calculate MAT for that year
  mat <- addLayer(mat, int_mat) #add MAT for that year to a stack
  int_mtwa <- max(sub) #calculate MTWA for that year
  mtwa <- addLayer(mtwa, int_mtwa) #add MTWA for that year to a stack
  print(a)
}
nlayers(av)
nlayers(mat)
nlayers(mtwa)
names(av)<-year
names(mat)<-year
names(mtwa)<-year

writeRaster(av, "data/processed/av.nc")
writeRaster(mat, "data/processed/mat.nc")
writeRaster(mtwa, "data/processed/mtwa.nc")

# EXTRACT POPULATION-LEVEL CLIMATE ----
temp_df <- cbind(Edata$Longitude, Edata$Latitude) %>% #put coordinates together
  sp::SpatialPoints(proj4string=crs(temp)) %>%  #create spatial points from coordinates
  raster::extract(temp, ., method = "bilinear") %>% #extract temperatures for each sampling location
  as.data.frame() #create df

#FIX NAs with values from closest non-NA cell
coords_withNAs <- cbind(Edata[which(is.na(temp_df[,1]), arr.ind=TRUE),]$Longitude, 
                        Edata[which(is.na(temp_df[,1]), arr.ind=TRUE),]$Latitude) #find coords with NAs
temp_df_fixed <- temp_df #new DF with corrected temperatures
#this loop finds the closest non-NAs cell for a given coordinate and repeats the same process
#for each raster layer (each represents a month and year), then it adds the new values
#to the DF
for (i in 1:nlayers(temp)) {
  r <- temp[[i]]
  sampled <- apply(X = coords_withNAs,MARGIN = 1,FUN = function(coords_withNAs) r@data@values[which.min(replace(distanceFromPoints(r, coords_withNAs), is.na(r), NA))])
  temp_df_fixed[which(is.na(temp_df[,1]), arr.ind=TRUE),i] <- sampled
}

#saving results
write.csv(temp_df_fixed, 'data/processed/temp_df_fixedv2.csv')

## CALCULATING POPULATION-LEVEL CLIMATE ----
# create new columns to store the climate variables
Edata $ mtwa_population <- NA
Edata $ mat_population <- NA
Edata $ av_population <- NA

for (i in 1:nrow(Edata)) {
  Edata$mtwa_population[i] <- max(NaRV.omit(temp_df_fixed[i,])) #mtwa
  #loop to calculate climate variables for each year until sample collection date
  #if there is no collection date, the paper publication year will be used instead
  year <- 1960:ifelse(is.na(Edata$Date), Edata$Publication, Edata$Date)[i]
  av_year <- rep(NA,length(year)) #to store temporary AV
  mat_year <- rep(NA,length(year)) #to store temporary MAT
  for (a in 1:length(year)) {
    z <- as.integer(year[a])
    av_year[a] <- max(NaRV.omit(unlist(temp_df_fixed[i,grep(z, names(temp_df_fixed))])))-
      min(NaRV.omit(unlist(temp_df_fixed[i,grep(z,names(temp_df_fixed))])))
    mat_year[a] <- mean(as.numeric(NaRV.omit(unlist(temp_df_fixed[i,grep(z,names(temp_df_fixed))]))))
  }
  Edata$av_population[i] <- mean(NaRV.omit(av_year)) #AV
  Edata$mat_population[i] <- mean(NaRV.omit(mat_year)) #MAT mean
  print(i)
}

Edata$difference_population <- Edata$Temperature-Edata$mat_population
Edata$magnitude_population <- Edata$difference_population/Edata$av_population

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

##CALCULATING SPECIES-LEVEL CLIMATE ---- only do this once for each combination of species and experiment year
EdataUn <- data.frame(Genus=Edata$Genus,Species=Edata$Species,Date=Edata$Date)
EdataUn$Date[which(is.na(Edata$Date))] <- Edata$Publication[which(is.na(Edata$Date))]
EdataUn <- unique(EdataUn)
EdataUn <- EdataUn[order(EdataUn$Species,EdataUn$Date),]
EdataUn$av_species <- NA
EdataUn$mat_species <- NA
EdataUn$mtwa_species <- NA

year<-1960:2019
focalSpecies<-unique(EdataUn$Species)
for(ii in 1:length(focalSpecies)){
  focalRow<-which(EdataUn$Species==focalSpecies[ii])
  
  spp <- as.character(paste(EdataUn$Genus, EdataUn$Species, sep = ' ')[focalRow[1]]) #extract species name
  poly <- seagrasses[seagrasses$binomial == spp,]
  #crop and mask the rasterbrick to the species distribution and downscale by 10 to avoid losing pixels at edge of range
  cropav <- raster::crop(av,poly)
  cropav<-disaggregate(cropav, fact=10)
  cropav<-raster::mask(cropav,poly)
  
  cropmat <- raster::crop(mat,poly)
  cropmat<-disaggregate(cropmat, fact=10)
  cropmat<-raster::mask(cropmat,poly)
  
  cropmtwa <- raster::crop(mtwa,poly)
  cropmtwa<-disaggregate(cropmtwa, fact=10)
  cropmtwa<-raster::mask(cropmtwa,poly)
  
  for(i in 1:length(focalRow)){
    
    #subset rasters to years prior to experiment
    
    focyear <- 1960:EdataUn$Date[focalRow[i]]
    num <- match(focyear,year) #find all layers from a specific year
    av_species<-subset(cropav,num)
    mat_species<-subset(cropmat,num)
    mtwa_species<-subset(cropmtwa,num)
    
    EdataUn$av_species[focalRow[i]] <- mean(NaRV.omit(values(av_species))) # mean of annual range in temperature across years and sites
    
    mat_species <- max(na.omit(mat_species)) # maximum of the mean annual temperatures across years for each site 
    mtwa_species <- max(na.omit(mtwa_species)) # maximum of the maximum monthly temperatures across years for each site 
    
    EdataUn$mat_species[focalRow[i]] <- mean(NaRV.omit(values(mat_species))) # mean across sites of the maximum of the mean annual temperatures across years
    EdataUn$mtwa_species[focalRow[i]] <- mean(NaRV.omit(values(mtwa_species)))
   
    print(focalRow[i])
    write.csv(EdataUn, 'data/processed/data_calculated_Unique.csv')
  }
}

# FULL DATASET ----

Edata$DateToUse<-Edata$Date
Edata$DateToUse[which(is.na(Edata$Date))]<-Edata$Publication[which(is.na(Edata$Date))]

Edata $ av_species <- NA
Edata $ mat_species <- NA
Edata $ mtwa_species <- NA
Edata $ difference_species <- NA
Edata $ magnitude_species <- NA

for (i in 1:nrow(Edata)){
  focRow<-which(EdataUn$Date==Edata$DateToUse[i] & EdataUn$Species==Edata$Species[i])
  Edata$av_species[i] <- EdataUn$av_species[focRow]
  Edata$mat_species[i]<-EdataUn$mat_species[focRow]
  Edata$mtwa_species[i]<-EdataUn$mtwa_species[focRow]
  
}  

Edata$difference_species<-Edata$Temperature-Edata$mat_species

Edata$magnitude_species <- Edata$difference_species/Edata$av_species

##SAVE ----
write.csv(Edata, 'data/processed/data_calculatedv2.csv')