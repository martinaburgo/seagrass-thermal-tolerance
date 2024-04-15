## READ ME ----
#This file provides to code used for each plot in the manuscript and supplementary information 

###LOAD PACKAGES ----
library(ncdf4)
library(raster)
library(sp)
library(maptools)
library(scales)
library(lme4)
library(gridExtra)
library(maps)
library(ggExtra)
library(rgdal)
library(viridis)
library(tidyverse)
library(IDPmisc)

###LOAD DATA ----
data <- read.csv('clean_data/data_calculated.csv') %>%
  dplyr::select(!(X))
heatwaves <- read.csv('clean_data/heatwaves_calculated.csv') %>%
  dplyr::select(!(X))
world <- map_data("world") #basemap
projections_df <- read.csv('projections.csv') %>%
  dplyr::select(!(X))
seagrasses <- readOGR('raw_data/SEAGRASSES/SEAGRASSES.shp') #seagrass species distributions
seagrasses@data $ type <- c('Temperate', 'Temperate', 'Temperate', 'Temperate', #4
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
seagrasses <- seagrasses[seagrasses@data$presence == 1,] #keep only presence = 1 (extant species on IUCN)
seagrasses <- seagrasses[seagrasses$binomial != "Ruppia polycarpa",]  #removing problematic extents
species <- seagrasses@data[,c(2,28)] #DF with only spp name and type
env_current <- raster('env_current.tif')
magnitude <- median(unlist(aggregate(heatwaves$difference_population_mean,
                                     by = c(heatwaves["Location"], heatwaves['Date']), 
                                     FUN = mean)[3]))/env_current

seagrasses_effort <- raster('seagrasses_effort.tif')
seagrasses_presence <- raster('seagrasses_presence.tif')

#getting the studies coords
data.x <- data$Longitude
data.y <- data$Latitude
data.labs <- data.frame(data.x, data.y)

#getting observed heatwaves coords
heat.x <- heatwaves$Longitude
heat.y <- heatwaves$Latitude
heat.labs <- data.frame(heat.x,heat.y)

#models
model0 <- glmer(formula = "survival_adjusted ~ magnitude_population_mean + 
                magnitude_species_1_mean + LOGtime + (1 | Species) + (1 | Study)", 
                data = data, family = binomial)
model1 <- glmer(formula = "survival_adjusted ~ magnitude_population_mean + 
                magnitude_species_1_mean + Type + LOGtime + (1 | Species) + (1 | Study)", 
                data = data, family = binomial)
model2 <- glmer(formula = "survival_adjusted ~ difference_population_mean + 
                av_population + LOGtime + (1 | Species) + (1 | Study)", 
                data = data, family = binomial)
model3 <- glmer(formula = "survival_adjusted ~ difference_species_1_mean + 
                av_species + LOGtime + (1 | Species) + (1 | Study)", 
                data = data, family = binomial)
model4 <- glmer(formula = "survival_adjusted ~ magnitude_population_mean + 
                magnitude_species_1_mean + mtwa_species_2_max + LOGtime + (1 | Species) + (1 | Study)", 
                data = data, family = binomial)

temp <- raster::stack("raw_data/sst.mon.mean.nc") %>% #load
  .[[((1960-1850)*12+1):nlayers(.)]] %>% #select period 1960 - present
  raster::rotate() #adjust longitude to -180,180

###MAIN----
##FIG. 1 - conceptual figure ----
ji <- function(xy, origin=c(0,0), cellsize=c(4,4)) {
  t(apply(xy, 1, function(z) cellsize/2+origin+cellsize*(floor((z - origin)/cellsize))))
}

#FIG. 1c
#loading data
z_muelleri <- subset(data, data$Species == 'muelleri') #subset only for z_muelleri
temp <- raster::stack("raw_data/sst.mon.mean.nc") %>% #load
  .[[((1960-1850)*12+1):nlayers(.)]] %>% #select period 1960 - present
  raster::rotate() #adjust longitude to -180,180

#extracting data for single study, local scale
temp_df <- cbind(z_muelleri$Longitude[1], z_muelleri$Latitude[1]) %>% #put coordinates together
  sp::SpatialPoints(proj4string=crs(temp)) %>%  #create spatial points from coordinates
  raster::extract(temp, ., method = "bilinear") %>% #extract temperatures for each sampling location
  as.data.frame() #create df

par(cex = 1.6)
data_plot <- data.table::transpose(temp_df) %>%
  rename(temperature = V1) %>%
  round(2) %>%
  mutate(year = substr(names(temp_df), 2, 5) %>%
           as.numeric()) %>%
  .[1:660,] %>% #until collection date
  ts(data = .[,1], start = 1960, end = 2014, frequency = 12) %>%
  plot(ylab = "Temperature (ºC)", xlab = "Year")
par(cex = 1, cex.main = 1)

#FIG. 1b
#extract data
spp <- as.character(paste(z_muelleri$Genus, z_muelleri$Species, sep = ' ')[1])
poly <- seagrasses[seagrasses$binomial == spp,]  #removing problematic extents
#crop the 'poly' rasterbrick with temperatures to the species distribution
crop <- raster::crop(temp,poly) %>%
  raster::mask(.,poly) %>%
  .[[1:660]]
mat_species <- raster(nrow = nrow(crop), ncol = ncol(crop), 
                      xmn = xmin(crop), xmx = xmax(crop), 
                      ymn = ymin(crop), ymx = ymax(crop),
                      crs = crs(crop))
for (a in 1960:2014) {
  num <- grep(as.integer(a), names(crop)) #find all layers from a specific year
  sub <- crop[[num]] #subset the stack to only that year
  int_mat <- mean(na.omit(sub)) #calculate MAT for that year
  mat_species <- addLayer(mat_species, int_mat) #add MAT for that year to a stack
  }
MAT <- mean(na.omit(mat_species))

#map
df <- as.data.frame(MAT, xy = T) %>%
  na.omit()
aus_nz <- extent(poly)+5
aus_nz_poly <- subset(world, world$long > aus_nz@xmin &
                        world$long < aus_nz@xmax &
                        world$lat > aus_nz@ymin &
                        world$lat < aus_nz@ymax)
poly_df <-  as.data.frame(poly, xy = T) %>%
  na.omit()

rast_map <- ggplot(df, aes(x, y)) + 
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     text=element_text(size=24,  family="Helvetica"),
                     legend.direction = 'horizontal',
                     legend.position = c(0.4,0.11)) + 
  geom_tile(aes(fill=layer)) +
  scale_fill_viridis(option = "plasma", name = "MAT (ºC)") +
  xlab("Longitude") + ylab("Latitude")

rast_map + geom_polygon(data = aus_nz_poly, 
                        aes(x=long, y = lat, group = group), fill = 'grey85') +
  geom_polygon(data = aus_nz_poly, 
               aes(x=long, y = lat, group = group), col = 'black', fill = NA, size = 0.2) +
  geom_polygon(data = poly, aes(x=long, y = lat, group = group), col = 'black', fill = NA, 
               size = 0.3) +
  geom_point(data = z_muelleri, aes(x = Longitude, y = Latitude, size = 4),
             col= '#440154FF', alpha = 0.4)+ guides(size = F, alpha= F)


#FIG. 1a
#map
fig_1a_base <- ggplot() + geom_polygon(data = world, aes(x=long, y = lat, 
                                                     group = group), fill = 'grey85') +
  geom_polygon(data = world, aes(x=long, y = lat, 
                                 group = group), col = 'black', fill = NA, size = 0.2) + 
  coord_fixed(1.3) + theme_bw() + theme(panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(),
                                        text=element_text(size=24,  family="Helvetica"),
                                        legend.direction = 'horizontal',
                                        legend.position = c(0.5,0.11)) + 
  xlab("Longitude") + ylab("Latitude")

sample_loc_aggregate <- cbind(data[,6:7], data %>%
                                dplyr::select(Longitude, Latitude) %>%
                                rename(X = Longitude, Y = Latitude) %>%
                                cbind() %>%
                                ji() %>%
                                as.data.frame() %>%
                                mutate(Cell = paste(X, Y)))
fig_1a_base +  by(sample_loc_aggregate, sample_loc_aggregate$Cell,
                  function(d) c(d$X[1], d$Y[1], nrow(d))) %>%
  unlist() %>%
  matrix(nrow = 3) %>%
  as.data.frame() %>%
  data.table::transpose() %>%
  rename(X = V1, Y = V2, Cell = V3) %>%
  geom_point(mapping = aes(X, Y, size = Cell, color = Cell, alpha = 0.6)) +
  scale_size_continuous(range=c(4,12), breaks = c(5,15,35,50)) +
  scale_color_viridis(option = "viridis", breaks = c(5,15,35,50),
                      name = 'Number of \nobservations') + guides(size = F, alpha = F)

##FIG. 2 - survival vs magnitude + survival over temp difference and AV ----
#FIG. 2a
figure_2a <- expand.grid(LOGtime = seq(min(data$LOGtime), 
                          max(data$LOGtime), 
                          length.out = 50),
            magnitude_species_1_mean = seq(1, 
                                           max(data$magnitude_species_1_mean), 
                                           length.out = 50),
            magnitude_population_mean = mean(data$magnitude_population_mean)) %>%
  mutate(Survival = predict(model0, ., type = 'response', re.form = NA)) %>%
  ggplot(aes(magnitude_species_1_mean, LOGtime, fill= Survival)) + 
  geom_tile() + theme_bw() + theme(panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   text=element_text(size=24,  family="Helvetica"),
                                   legend.position = c(0.15, 0.22)) +
  scale_fill_viridis(option = "plasma", name = "Survival") + 
  ggtitle("Survival over time and magnitude - species level") +
  xlab("Magnitude of warming") + ylab("Duration (Hours, log scale)") 

#FIG. 2b
figure_2b <- data.frame(difference_population_mean = c(rep(seq(0,35),3)),
                        av_population = c(rep(5, length.out = 36),
                                          rep(10, length.out = 36),
                                          rep(15, length.out = 36)),
                        LOGtime = log(24*14)) %>% #2-week MHW
  mutate(Survival = predict(model2, ., type = 'response', re.form = NA)) %>%
  ggplot(aes(difference_population_mean, Survival, 
             color = as.factor(av_population))) + 
  geom_line(size = 1.2) + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=24,  family="Helvetica"),
        legend.position = c(0.15, 0.15)) +
  scale_color_viridis(option = "viridis", name = "AV (ºC)", discrete = T) + 
  ggtitle("Survival over temp difference - pop level") +
  xlab("Temperature difference (ºC)") + ylab("Survival")

#final
grid.arrange(figure_2a, figure_2b, nrow = 1)

data.frame(difference_population_mean = c(rep(seq(0,35),3)),
           av_population = c(rep(5, length.out = 36),
                             rep(10, length.out = 36),
                             rep(15, length.out = 36)),
           LOGtime = log(24*14)) %>% #2-week MHW
  mutate(Survival = predict(model2, ., type = 'response', re.form = NA)) %>%
  filter(difference_population_mean == 15) %>%
  dplyr::select(av_population, Survival)

##FIG. 3 - how predicted survival was calculated ----
IDs <- read.csv('heatwaves_predicted.csv') %>%
  dplyr::select(ID) %>%
  unique()
plots <- list()

for (a in 1:nrow(IDs)) {
  subset <- read.csv('heatwaves_predicted.csv') %>%
    filter(ID == IDs[a,])
  subset $ time <- exp(subset$LOGtime)
   time_series <- 0
   for (b in 1:length(subset$time)) {
     point <- subset$time[b]/24
     time_point <- point + time_series[b]
     time_series <- c(time_series, time_point)
     } #this loop calculates the 'breaks'/days since beginning of MHW when temperatures changed
  
  time <- time_series
  
  time <- c(0,cumsum(subset$time/24))
  
  survival <- c(1, subset$prediction)
  cummSurvival <- cumprod(survival)
  
  temp <- subset$Temperature
  temp <- c(subset$mat_species_1_mean[1], temp)
  
  tempScales<-scales::rescale(temp, to=c(0,1))
  tempScales<-tempScales[-1]
  
  tags <- paste(subset$Genus[1], subset$Species[1], '-', subset$Location[1], 
                '-', subset$Date[1],
                sep = ' ')
  
  par(mar=c(5,4,4,6) + 1)
  par(new=TRUE)
  plot(time,cummSurvival,type="n", ylim = c(0,1), 
       xlim = c(0,max(time_series)),
       xlab = 'Time since start of MHW (days)',
       ylab = 'Survival',
       main = tags,
       cex.lab = 2, cex.axis = 2, cex.main = 2)
  for(i in 1:length(time)){
    polygon(c(time[i],time[i+1],time[i+1],time[i]),
            c(0,0,tempScales[i],tempScales[i]),col="orange",border=FALSE)
  }  
  points(time,cummSurvival,type="l")
  points(time,cummSurvival,pch=21,bg="white", cex = 2.5)
  points(max(time),subset$Survival[1],pch=21,bg="black", cex = 2.5)
  
  axis(side=4, at = pretty(range(scales::rescale(temp, to=c(0,1)))),
       labels = c(round(seq(min(temp), max(temp), length.out = 6), 0)), cex.axis = 2)
  mtext("Temperature (ºC)", side=4, line=4, cex = 2)
  
  plots[[a]] <- recordPlot()
  graphics.off()
}

ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group), 
                        fill = 'grey85') + theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(size=24,  family="Helvetica"))  +
  xlab("Longitude") + ylab("Latitude") +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), 
               col = 'black', fill = NA, size = 0.2) +
  geom_point(data = heat.labs, aes(x = heat.x, y = heat.y), 
             color = "orange2", size = 4)

##FIG. 4 - observed survival vs predicted survival ----
heatwaves %>%
  group_by(ID) %>%
  mutate(overall_prediction = prod(prediction)) %>%
  mutate(duration = sum(exp(LOGtime))) %>%
  distinct(ID,.keep_all = T) %>%
  ggplot(mapping = aes(Survival, overall_prediction)) + geom_point(size = 3, 
                                                                   position = 'jitter') + 
  geom_abline() + xlim(0,1) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     text=element_text(size=24,  family="Helvetica")) +
  xlab("Observed survial") + ylab("Predicted survival")

##FIG. 5 - most vulnerable areas and current magnitude ----
#fig. 5a
magnitude %>%
  as.data.frame(xy = T) %>%
  na.omit() %>%
  ggplot(aes(x, y)) + theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(size=24,  family="Helvetica")) + 
  geom_tile(aes(fill=env_current), show.legend = T) +
  scale_fill_viridis(limits = c(0, 12), oob = squish, name = "Magnitude") +
  xlab("Longitude") + ylab("Latitude") + ggtitle("Magnitude - population") +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), 
               fill = 'grey85') +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), 
               col = 'black', fill = NA, size = 0.2)

#fig. 5b
raster('averages.tif') %>%
  as.data.frame(xy = T) %>%
  na.omit() %>%
  ggplot(aes(x, y)) + theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(size=24,  family="Helvetica")) + 
  geom_tile(aes(fill=averages), show.legend = T) + 
  scale_fill_gradientn(colours = c('#a50026','#d73027','#f46d43','#fdae61',
                                   '#fee08b','#ffffbf','#d9ef8b','#a6d96a',
                                   '#66bd63','#1a9850','#006837'), name = "Survival") +
  xlab("Longitude") + ylab("Latitude") + ggtitle("Current threat") +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), 
               fill = 'grey85') +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), 
               col = 'black', fill = NA, size = 0.2)

###EXTENDED DATA FIG ----
##Extended Data Fig. 1 - survival vs mtwa species ----
data.frame(mtwa_species_2_max = rep(seq(min(data$mtwa_species_2_max),
                                        max(data$mtwa_species_2_max),
                                        length.out = 40),3),
           magnitude_species_1_mean = c(rep(2, 40),
                                        rep(4, 40),
                                        rep(6, 40)),
           magnitude_population_mean = rep(mean(data$magnitude_population_mean),
                                           120),
           level = c(rep('Low', 40),
                     rep('Medium', 40),
                     rep('High', 40)),
           LOGtime = rep(log(24*14), 120)) %>% #2-week MHW
  mutate(Survival = predict(model4, ., type = 'response', re.form = NA)) %>%
  ggplot(aes(mtwa_species_2_max, Survival, col = level)) + 
  geom_line(size = 1.2) + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=24,  family="Helvetica"),
        legend.position = c(0.25, 0.86)) +
  scale_color_viridis(option = "plasma", name = "Magnitude (ºC)", discrete = T,
                      breaks = c('Low', 'Medium', 'High')) +
  xlab("MTWA (ºC)") + ylab("Survival")

##Extended Data Fig. 2 - survival over magnitude (population-level) ----
quantile(data$magnitude_population_mean)
data.frame(magnitude_species_1_mean = rep(mean(data$magnitude_species_1_mean), 150),
           magnitude_population_mean = round(c(rep(1, 50),
                                               rep(3, 50),
                                               rep(5, 50)),2),
           LOGtime = rep(seq(min(data$LOGtime), max(data$LOGtime), length.out = 50), 3)) %>% 
  mutate(Survival = predict(model0, ., type = 'response', re.form = NA)) %>%
  ggplot(aes(LOGtime, Survival, 
             color = as.factor(magnitude_population_mean))) + 
  geom_line(size = 1.2) + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=24,  family="Helvetica"),
        legend.position = c(0.28, 0.15)) +
  scale_color_viridis(option = "viridis", name = "Magnitude of warming (ºC)", discrete = T) + 
  ggtitle("Survival over mag - pop level") +
  xlab("Duration (Hours, log scale)") + ylab("Survival")

##Extended Data Fig. 3 - historical vs predicted with abiotic factors ----
heatwaves %>%
  group_by(ID) %>%
  mutate(overall_prediction = prod(prediction)) %>%
  mutate(duration = sum(exp(LOGtime))) %>%
  distinct(ID,.keep_all = T) %>%
  mutate(new_severity = ifelse(Severity == "N/A", NA,
                               as.character(Severity))) %>%
  ggplot() + geom_point(aes(Survival, overall_prediction, shape = Effect,
                            fill = new_severity, colour = new_severity),
             size = 3, position = 'jitter') + geom_abline() + xlim(0,1) + ylim(0,1) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=24,  family="Helvetica"),
        legend.position = c(0.85, 0.35)) +
  xlab("Observed survial") + ylab("Predicted survival") + 
  scale_shape_manual(values = c(24, 21, 25),
                     name = 'Effect on \n survival',
                     breaks = c('Positive','Neutral', 'Negative')) + 
  scale_fill_manual(values = c('High' = '#d73027',
                               'Medium' = '#fee08b',
                               'Low' = '#91cf60'),
                    name = 'Severity of \n abiotic factors',
                    na.translate = T, na.value = 'grey',
                    breaks = c('High', 'Medium', 'Low', NA),
                    aesthetics = "fill") +
  scale_colour_manual(values = c('High' = '#d73027',
                                 'Medium' = '#fee08b',
                                 'Low' = '#91cf60'),
                      name = 'Severity of \n abiotic factors',
                      na.translate = T, na.value = 'grey',
                      breaks = c('High', 'Medium', 'Low', NA),
                      aesthetics = "colour")
                     

###SUPPLEMENTARY INFORMATION ----
##Supplementary Fig. 1-2 - data distribution ----
#SI Figure 1
p <- ggplot(data, aes(Temperature, LOGtime, colour = Bioregion)) + 
  geom_point(size = 3) + theme_bw() + theme(panel.grid.major = element_blank(), 
                                            panel.grid.minor = element_blank(),
                                            text=element_text(size=16,  family="Helvetica"),
                                            legend.position = c(0.2, 0.2)) + 
  ylab("Duration (Hours, log scale)")
ggMarginal(p, type="histogram") 
  

#SI Figure 2
s <- data %>%
  select(Genus, Species) %>%
  mutate(name = paste0(substr(Genus, 1, 1), ". ", Species)) %>%
  select(name) %>%
  table() %>%
  as_data_frame()
names(s) <- c('spp', 'freq')

ggplot(s,aes(x=spp, y=freq)) + 
  geom_bar(stat = "identity") + theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(size=16,  family="Helvetica"),
        axis.text.x = element_text(angle = 90, colour = 'black')) +
  xlab('') + ylab('Number of observations')


##Supplementary Fig. 3 - sampling effort ----
seagrasses_effort_raster <- seagrasses_effort/seagrasses_presence*100

sampling_effort_map <- seagrasses_effort_raster %>% #calculate sampling effort
  as.data.frame(xy = T) %>%
  na.omit() %>%
  ggplot(aes(x, y)) + theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(size=16,  family="Helvetica")) +
  geom_tile(fill = 'grey') +
  xlab("Longitude") + ylab("Latitude") + 
  geom_tile(sampling_effort_df, mapping = aes(x, y, fill=layer)) + 
  scale_fill_gradientn(colours = c('#d73027','#f46d43','#fdae61','#fee090',
                                   '#e0f3f8','#abd9e9','#74add1','#4575b4'), 
                       name = "% of species") + 
  geom_polygon(data = world, aes(x=long, y = lat, group = group), 
               fill = 'grey85') +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), 
               col = 'black', fill = NA, size = 0.2) +
  geom_point(data = data.labs, aes(x = data.x, y = data.y), 
             color = "limegreen", size = 3)

seagrasses_presence_map <- seagrasses_presence %>%
  as.data.frame(xy = T) %>%
  na.omit() %>%
  ggplot(aes(x, y)) + theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(size=16,  family="Helvetica")) + 
  geom_tile(aes(fill=seagrasses_presence)) + 
  scale_fill_gradientn(colours = c('#f7fcb9','#d9f0a3','#addd8e',
                                   '#78c679','#41ab5d','#238443',
                                   '#006837','#004529'), 
                       name = "N. of species") +
  xlab("Longitude") + ylab("Latitude") + 
  geom_polygon(data = world, aes(x=long, y = lat, group = group), 
               fill = 'grey85') +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), 
               col = 'black', fill = NA, size = 0.2) +
  geom_point(data = data.labs, aes(x = data.x, y = data.y), 
             color = "limegreen", size = 3)

grid.arrange(seagrasses_presence_map,
             sampling_effort_map,
             nrow = 2)

##Supplementary Fig. 4 - current AV ----
env_current %>% 
  as.data.frame(xy = T) %>%
  na.omit() %>%
  ggplot(aes(x, y)) + geom_tile(aes(fill=env_current), show.legend = T) + theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(size=24,  family="Helvetica")) +
  scale_fill_viridis(option = "viridis", name = "AV (ºC)") +
  xlab("Longitude") + ylab("Latitude") + ggtitle("Current AV") +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), 
               fill = 'lightgrey') +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), 
               col = 'black', fill = NA, size = 0.2)

##Supplementary Fig. 5 - current magnitude risk ----
magnitudes_list <- list() #create list to store the rasters for each magnitude + spp
for (i in 1:nrow(species)) {
  spp <- as.character(species$binomial[i]) #get species name
  poly <- seagrasses[seagrasses@data$binomial == spp,] #subset spp polygon
  crop <- raster::crop(magnitude,poly) %>%
    raster::mask(.,poly)
  magnitudes_list[[i]] <- crop
}

#calculating average magnitude for each grid cell
magnitudes_stack <- stack()
for (i in 1:length(magnitudes_list)) {
  xx <- resample(magnitudes_list[[i]], env_current) #adjust the extent so that it's the same
  magnitudes_stack <- addLayer(magnitudes_stack, xx)
}

mean_magnitudes <- stackApply(magnitudes_stack,
                              indices = c(rep(1, nlayers(magnitudes_stack))),
                              fun = na.omit(mean))

mean_magnitudes %>% 
  as.data.frame(xy = T) %>%
  na.omit() %>%
  ggplot(aes(x, y)) + geom_tile(aes(fill=index_1), show.legend = T) + theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(size=24,  family="Helvetica")) +
  scale_fill_viridis(option = "viridis", name = "Average\nmagnitude") +
  xlab("Longitude") + ylab("Latitude") + ggtitle("Magnitude - species") +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), 
               fill = 'lightgrey') +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), 
               col = 'black', fill = NA, size = 0.2)

##Supplementary Fig. 6 - distribution of tropical and temperate seagrass species ----
raster_number <- raster(nrows = 180, ncols = 360, #create a raster as a reference to 'rasterize' the polygons
                        xmn = -180, xmx = 180, 
                        ymn = -90, ymx = 90)
crs(raster_number) <- '+proj=longlat +datum=WGS84 +no_defs'

trop_presence_raster <- seagrasses[seagrasses@data$type == 'Tropical',] %>%
  rasterize(raster_number, field = 1, 
            fun = sum, na.rm = T)
perc_trop <- trop_presence_raster/seagrasses_presence*100

prop_trop_spp_raster <- seagrasses_presence %>%
  reclassify(c(0, 19, 0)) %>%
  raster::merge(perc_trop, .)

prop_trop_spp_raster %>% 
  as.data.frame(xy = T) %>%
  na.omit() %>%
  ggplot(aes(x, y)) + geom_tile(aes(fill=layer), show.legend = T) + theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(size=24,  family="Helvetica")) +
  scale_fill_viridis(option = "plasma", name = "Tropical\nspecies (%)") +
  xlab("Longitude") + ylab("Latitude") +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), 
               fill = 'lightgrey') +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), 
               col = 'black', fill = NA, size = 0.2)

layerStats(stack(raster('averages.tif'), prop_trop_spp_raster,
                 mean_magnitudes, reclassify(magnitude, c(12, Inf, 12))), 
           'pearson', na.rm = T)
grid.arrange(projections_df %>%
  ggplot(aes(x=type, y=av_species)) +
  geom_boxplot() + theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(size=24,  family="Helvetica")) +
  xlab("") + ylab("AV (ºC)"), projections_df %>%
  ggplot(aes(x=type, y=mtwa_species)) +
  geom_boxplot() + theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(size=24,  family="Helvetica")) +
  xlab("") + ylab("MTWA (ºC)"), ncol = 2)
