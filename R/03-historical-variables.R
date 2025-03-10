##LOAD PACKAGES ----
library(raster)
library(tidyverse)
library(sf)
library(IDPmisc)

source('R/functions.R')

#LOAD DATA ----
dieoffs <- readxl::read_excel('data/primary/seagrass-heatwaves-studies.xlsx', sheet = 2)
Edata <- dieoffs[, c("Event ID", "Authors", "Title of the publication", "Year of publication", 
                     "Year of MHW","Location of MHW", "Latitude of MHW", "Longitude of MHW", 
                     "Family", "Genus", "Species", "Type",
                     "Duration of MHW period (hours)", "Temperature of MHW period (ºC)",
                     "Accepted response variable", "Survival")]
colnames(Edata) <- c("ID", "Study", "Title", "Publication", "Date", "Location",
                     "Latitude", "Longitude", "Family", "Genus", "Species",
                     "Type", "Time", "Temperature", "Variable", "Survival")

seagrasses <- sf::st_read('data/primary/SEAGRASSES/SEAGRASSES.shp') #seagrass species distributions
seagrasses $ type <- c('Temperate', 'Temperate', 'Temperate', 'Temperate', #4
                       'Temperate', 'Tropical', 'Tropical', 'Temperate', #8
                       'Temperate', 'Temperate', 'Temperate', 'Temperate', #12
                       'Temperate', 'Tropical', 'Tropical', 'Tropical', #16
                       'Tropical', 'Tropical', 'Tropical', 'Tropical', #20
                       'Tropical', 'Tropical', 'Temperate', 'Tropical', #24
                       'Temperate', 'Temperate', 'Tropical', 'Tropical', #28
                       'Temperate', 'Tropical', 'Temperate', 'Tropical', #32
                       'Tropical', 'Temperate', 'Temperate', 'Tropical', #36
                       'Tropical', 'Temperate', 'Temperate', 'Temperate', #40
                       'Temperate', 'Tropical', 'Temperate', 'Temperate', #44
                       'Temperate', 'Temperate', 'Temperate', 'Tropical', #48
                       'Temperate', 'Temperate', 'Temperate', 'Temperate', #52
                       'Tropical', 'Tropical', 'Tropical', 'Temperate', #56
                       'Temperate', 'Tropical', 'Temperate', 'Temperate', #60
                       'Temperate', 'Tropical', 'Temperate', 'Tropical', #64
                       'Tropical', 'Tropical', 'Temperate', 'Tropical', #68
                       'Temperate', 'Temperate', 'Temperate', 'Tropical', #72
                       'Temperate', 'Temperate')
seagrasses <- seagrasses[seagrasses$presence == 1,] #keep only presence = 1 (extant species on IUCN)
seagrasses <- seagrasses[seagrasses$binomial != "Ruppia polycarpa",]  #removing problematic extents

mat <- raster::stack("data/processed/mat.nc")
names(mat) <- 1960:2019
av <- raster::stack("data/processed/av.nc")
names(av) <- 1960:2019
mtwa <- raster::stack("data/processed/mtwa.nc")
names(mtwa) <- 1960:2019

temp_df_fixed <- read.csv('data/processed/temp_df_fixedv2.csv')

# Calculate climate variables ----
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

# CALCULATING SPECIES-LEVEL CLIMATE ---- 
#only do this once for each combination of species and experiment year
EdataUn <- data.frame(Genus=Edata$Genus,Species=Edata$Species,Date=Edata$Date)
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
    EdataUn$mtwa_species[focalRow[i]] <- max(NaRV.omit(values(mtwa_species)))
    
    print(focalRow[i])
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
write.csv(Edata, 'data/processed/heatwaves_calculated.csv')
