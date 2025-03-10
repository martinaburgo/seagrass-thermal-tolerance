# READ ME ----
#This file provides to code used for each plot in the manuscript and supplementary information 

# Set up ----
## Packages ----
library(ncdf4)
library(raster)
library(sf)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(maptools)
library(scales)
library(lme4)
library(ggplot2)
library(gridExtra)
library(maps)
library(ggExtra)
library(rgdal)
library(emmeans)
library(tidyverse)
library(IDPmisc)

library(terra)         # For raster handling


## Data ----
data <- read.csv('data/processed/data_calculatedv2.csv') |> 
  dplyr::select(!(X)) |> 
  dplyr::filter(Temperature >= mat_population) |>
  dplyr::filter(Time >= 1) |> 
  dplyr::mutate(Bioregion = ifelse(Location == 'Peel-Harvey Estuary, WA', 'Temperate Southern Oceans',
                                   ifelse(Location == 'Swan River Estuary, Perth, WA', 'Temperate Southern Oceans',
                                          Bioregion))) #fix bioregions for two studies
heatwaves <- read.csv('data/processed/heatwaves_calculated.csv') |> 
  dplyr::select(!(X))
world <- map_data("world") #basemap
projections_df <- read.csv('projections.csv') |> 
  dplyr::select(!(X))
seagrasses <- st_read('data/primary/SEAGRASSES/SEAGRASSES.shp') #seagrass species distributions
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
species <- seagrasses |> dplyr::select(binomial) |> as.data.frame() #DF with only spp name and type
env_current <- raster('env_current.tif')
magnitude <- median(unlist(aggregate(heatwaves$difference_population_mean,
                                     by = c(heatwaves["Location"], heatwaves['Date']), 
                                     FUN = mean)[3]))/env_current

seagrasses_effort <- raster('seagrasses_effort.tif')
seagrasses_presence <- raster('seagrasses_presence.tif')

source('R/functions.R')

#getting the studies coords
data.x <- data$Longitude
data.y <- data$Latitude
data.labs <- data.frame(data.x, data.y)

#getting observed heatwaves coords
heat.x <- heatwaves$Longitude
heat.y <- heatwaves$Latitude
heat.labs <- data.frame(heat.x,heat.y)

#models
load('data/modelled/Model6Beta.RData')
load('data/modelled/Model7Beta.RData')
load('data/modelled/Model8Beta.RData')
load('data/modelled/Model9Beta.RData')

temp <- raster::stack("raw_data/sst.mon.mean.nc") %>% #load
  .[[((1960-1850)*12+1):nlayers(.)]] %>% #select period 1960 - present
  raster::rotate() #adjust longitude to -180,180

# Main figures ----
##FIG. 1 - conceptual figure ----
ji <- function(xy, origin=c(0,0), cellsize=c(4,4)) {
  t(apply(xy, 1, function(z) cellsize/2+origin+cellsize*(floor((z - origin)/cellsize))))
}

#FIG. 1c
#loading data
z_muelleri <- subset(data, data$Species == 'muelleri') #subset only for z_muelleri
temp <- raster::stack("data/primary/sst.mon.mean.nc") %>% #load
  .[[((1960-1850)*12+1):nlayers(.)]] %>% #select period 1960 - present
  raster::rotate() #adjust longitude to -180,180

#extracting data for single study, local scale
temp_df <- cbind(z_muelleri$Longitude[1], z_muelleri$Latitude[1]) %>% #put coordinates together
  sp::SpatialPoints(proj4string=crs(temp)) %>%  #create spatial points from coordinates
  raster::extract(temp, ., method = "bilinear") %>% #extract temperatures for each sampling location
  as.data.frame() #create df

ts_data <- data.table::transpose(temp_df) |> 
  dplyr::rename(temperature = V1) |> 
  round(2) |> 
  dplyr::mutate(year = substr(names(temp_df), 2, 5)  |> 
                  as.numeric()) %>% 
  .[1:660,] %>%  #until collection date
  ts(data = .[,1], start = 1960, end = 2014, frequency = 12) 
ts_data <- dplyr::tibble(Date = seq(from = as.Date("1960-01-01"), by = "month", length.out = length(ts_data)),
                         Temperature = as.numeric(ts_data))

 fig_1c <- ts_data |> 
  ggplot(aes(x = Date, y = Temperature)) +
  geom_line() +
  scale_x_date(breaks = seq(from = as.Date("1960-01-01"), 
                            to = as.Date("2014-01-01"), 
                            by = "10 years"),
               date_labels = "%Y") +  # Show only decades
  geom_hline(yintercept = mean(ts_data$Temperature), linetype = 'dashed') +
  geom_vline(xintercept = as.numeric(max(ts_data$Date)), size = 1.2) +
  # MTWA:
  geom_point(y = as.numeric(max(ts_data$Temperature)),
             x = ts_data$Date[which.max(ts_data$Temperature)] , size = 2) +
   # AV:
   geom_point(y = ts_data$Temperature[which(ts_data$Date == '2010-02-01')],
              x = ts_data$Date[which(ts_data$Date == '2010-02-01')], size = 2) +
   geom_point(y = ts_data$Temperature[which(ts_data$Date == '2010-07-01')],
              x = ts_data$Date[which(ts_data$Date == '2010-07-01')], size = 2) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      text=element_text(family="Helvetica"),
                      legend.direction = 'horizontal',
                      legend.position = c(0.4,0.11)) +
  labs(x = element_blank(), y = "Temperature (°C)", subtitle = 'c') +
  theme(plot.subtitle = element_text(hjust = 0.01, vjust = -5, margin = margin(2, 0, 0, 0), face = "bold"),
        plot.margin = margin(0, 120, 0, 0)) 

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
df <- as.data.frame(MAT, xy = T) |> 
  na.omit()
aus_nz <- extent(poly)+5
aus_nz_poly <- subset(world, world$long > aus_nz@xmin &
                        world$long < aus_nz@xmax &
                        world$lat > aus_nz@ymin &
                        world$lat < aus_nz@ymax)
poly_df <-  as.data.frame(poly, xy = TRUE)  |> 
  na.omit()

rast_map <- ggplot() + 
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     text = element_text(family = "Helvetica"),
                     legend.direction = 'horizontal',
                     legend.position = c(0.3,0.11),
                     legend.background = element_blank()) + 
  geom_tile(df, mapping = aes(x, y, fill=layer)) +
  scale_fill_viridis_c(option = "plasma", name = "MAT (ºC)")

fig_1b <- rast_map + geom_polygon(data = aus_nz_poly, 
                        aes(x=long, y = lat, group = group), fill = 'grey85') +
  geom_polygon(data = aus_nz_poly, 
               aes(x=long, y = lat, group = group), col = 'black', fill = NA, size = 0.2) +
  geom_sf(data =  sf::st_sf(poly), col = 'black', fill = NA, 
               size = 0.3) +
  geom_point(data = z_muelleri, aes(x = Longitude, y = Latitude, size = 3),
             col= 'black', alpha = 0.8)+ guides(size = F, alpha= F)  +
  ylab('Latitude') + xlab('Longitude') +
  labs(subtitle = 'b') +
  scale_x_continuous(labels = scales::label_number(accuracy = 1), n.breaks = 3) +  # Remove directional labels for Longitude
  scale_y_continuous(labels = scales::label_number(accuracy = 1), n.breaks = 3) +  # Remove directional labels for Latitude
  theme(plot.subtitle = element_text(hjust = 0.01, vjust = -5, margin = margin(2, 0, 0, 0), face = "bold")) 


#FIG. 1a
#map
fig_1a_base <- ggplot() + geom_polygon(data = world |> dplyr::filter(region != 'Antarctica'), aes(x=long, y = lat, 
                                                     group = group), fill = 'grey85') +
  geom_polygon(data = world |> dplyr::filter(region != 'Antarctica'), aes(x=long, y = lat, 
                                 group = group), col = 'black', fill = NA, size = 0.2) + 
  coord_fixed(1.3) + theme_bw() + theme(panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(),
                                        text=element_text(family="Helvetica"),
                                        legend.direction = 'horizontal',
                                        legend.position = c(0.15, 0.11),
                                        legend.background = element_blank())

sample_loc_aggregate <- cbind(data[,6:7], data |> 
                                dplyr::select(Longitude, Latitude) |> 
                                dplyr::rename(X = Longitude, Y = Latitude) |> 
                                cbind() |> 
                                ji() |> 
                                as.data.frame() |> 
                                dplyr::mutate(Cell = paste(X, Y)))
fig_1a <- fig_1a_base +  by(sample_loc_aggregate, sample_loc_aggregate$Cell,
                  function(d) c(d$X[1], d$Y[1], nrow(d))) |> 
  unlist() |> 
  matrix(nrow = 3) |> 
  as.data.frame() |> 
  data.table::transpose() |> 
  dplyr::rename(X = V1, Y = V2, Cell = V3) |> 
  geom_point(mapping = aes(X, Y, size = Cell, color = Cell, alpha = 0.7)) +
  scale_size_continuous(range = c(4, 12), breaks = c(5, 20, 35, 50)) +
  scale_color_viridis_c(option = "mako", breaks = c(5, 20, 35, 50), direction = -1,
                      name = 'Obs.') + guides(size = F, alpha = F)  +
  ylab('Latitude') + xlab('Longitude') +
  labs(subtitle = 'a') +
  theme(plot.subtitle = element_text(hjust = 0.01, vjust = -5, margin = margin(2, 0, 0, 0), face = "bold")) 

### Final Figure 1 ----
(fig_1a + fig_1b + 
   plot_layout(axis_titles = 'collect'))  / wrap_elements(fig_1c)
#save as 12x7 inch

##FIG. 2 - heatmap (temp and time) + survival over temp difference and AV (population) ----
#FIG. 2a
preds2a <- seagrass.brm7 |> 
  emmeans(~Time|difference_species, type = 'response', 
          at = list(Time = seq(min(seagrass$Time), max(seagrass$Time), length.out = 50),
                    difference_species = seq(0, 
                                             max(seagrass$difference_species), length.out = 50))) |> 
  dplyr::as_tibble()

figure_2a <- preds2a |> 
  ggplot(aes(x = Time, y = difference_species, fill = response)) + 
  geom_tile() + theme_bw() + theme(panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   text = element_text(size = 12,  family = "Helvetica"),
                                   legend.position = c(0.88, 0.78),
                                   plot.subtitle = element_text(hjust = 0.01, vjust = -5, 
                                                                margin = margin(2, 0, 0, 0), face = "bold")) +
  scale_fill_viridis_c(option = "plasma", name = "Survival", direction = -1) +
  scale_x_continuous(name = 'Duration (days)',
    breaks = c(7*24, 14*24, 30*24, 60*24, 90*24), # Set the tick positions
    labels = c(7, 14, 30, 60, 90) # Set the corresponding labels
  )  +
  labs(subtitle = 'a') + 
  #ggtitle("Survival over time and intensity - species level") +
  ylab("∆MAT (ºC)")

#FIG. 2b
preds2b <- seagrass.brm6 |> 
  emmeans(~av_population|difference_population|Time, type = 'response', 
          at = list(av_population = c(5, 15),
                    difference_population = seq(min(seagrass$difference_population), 
                                             max(seagrass$difference_population), length.out = 50),
                    Time = 14*24)) |> 
  dplyr::as_tibble()

figure_2b <- preds2b |> 
  ggplot(aes(difference_population, response, 
             color = as.factor(av_population),
             fill = as.factor(av_population))) + 
  geom_ribbon(mapping = aes(ymin = lower.HPD, ymax = upper.HPD), alpha = 0.2, colour = NA) +
  geom_line(size = 1) + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=12,  family="Helvetica"),
        legend.position = c(0.9, 0.87),
        plot.subtitle = element_text(hjust = 0.01, vjust = -5, margin = margin(2, 0, 0, 0), face = "bold")) +
  scale_color_viridis_d(option = "D", name = "AV (ºC)") + 
  scale_fill_viridis_d(option = "D", name = "AV (ºC)")  +
  labs(subtitle = 'b') + 
  #ggtitle("Survival over temp difference - pop level") +
  xlab("∆MAT (ºC)") + ylab("Survival")

#Fig 2c - MTWA effect
preds2c <- seagrass.brm6 |> 
  emmeans(~mtwa_population|av_population|difference_population|Time, type = 'response', 
          at = list(av_population = c(10),
                    difference_population = seq(min(seagrass$difference_population), 
                                                max(seagrass$difference_population), length.out = 50),
                    Time = c(14*24),
                    mtwa_population = c(20, 30))) |> 
  dplyr::as_tibble()

figure_2c <- preds2c |> 
  ggplot(aes(difference_population, response, 
             color = as.factor(mtwa_population),
             fill = as.factor(mtwa_population))) + 
  geom_ribbon(mapping = aes(ymin = lower.HPD, ymax = upper.HPD), alpha = 0.2, colour = NA) +
  geom_line(size = 1) + theme_bw() + 
  #facet_grid(~av_population) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=12,  family="Helvetica"),
        legend.position = c(0.87, 0.87),
        plot.subtitle = element_text(hjust = 0.01, vjust = -5, margin = margin(2, 0, 0, 0), face = "bold")) +
  scale_color_viridis_d(option = "D", name = "MTWA (ºC)") + 
  scale_fill_viridis_d(option = "D", name = "MTWA (ºC)")  +
  labs(subtitle = 'c') +
  #ggtitle("Survival over temp difference - pop level") +
  xlab("∆MAT (ºC)") + ylab("Survival")


#final
figure_2a +figure_2b + figure_2c + plot_layout(axes = 'collect')

##FIG. 3 - how predicted survival was calculated ----
IDs <- read.csv('data/processed//heatwaves_predicted.csv') |> 
  dplyr::select(ID) |> 
  unique()
plots <- list()

for (a in 1:nrow(IDs)) {
  subset <- read.csv('data/processed//heatwaves_predicted.csv') |> 
    filter(ID == IDs[a,])
  subset $ time <- exp(subset$Time)
   time_series <- 0
   for (b in 1:length(subset$Time)) {
     point <- subset$Time[b]/24
     time_point <- point + time_series[b]
     time_series <- c(time_series, time_point)
     } #this loop calculates the 'breaks'/days since beginning of MHW when temperatures changed
  
  time <- time_series
  
  time <- c(0,cumsum(subset$Time/24))
  
  survival <- c(1, subset$Predicted)
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
# Load world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Custom color palette
custom_colors <- c('#a50026','#d73027','#f46d43','#fdae61',
                   '#fee08b','#ffffbf','#d9ef8b','#a6d96a',
                   '#66bd63','#1a9850','#006837')

# Function to process raster stack (cap values at 1)
cap_raster_values <- function(raster_stack) {
  raster_stack[raster_stack > 1] <- 1
  return(raster_stack)
}

# Apply the function to both raster stacks
pred_1w_stack <- cap_raster_values(pred_1w_stack)
pred_2w_stack <- cap_raster_values(pred_2w_stack)

# Function to convert raster to data frame for ggplot
raster_to_df <- function(r, name) {
  r_df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  colnames(r_df) <- c("x", "y", "value")
  r_df$layer <- name
  return(r_df)
}

# Extract layers from raster stacks and convert to data frames
layers <- c("diff5", "diff10", "diff15")

df_list <- list()
for (i in seq_along(layers)) {
  df_list[[layers[i]]] <- list(
    pred_1w = raster_to_df(pred_1w_stack[[layers[i]]], layers[i]),
    pred_2w = raster_to_df(pred_2w_stack[[layers[i]]], layers[i])
  )
}

# Function to create a ggplot for each raster layer
plot_raster <- function(df, title) {
  ggplot() +
    geom_tile(data = df, aes(x = x, y = y, fill = value)) +
    geom_sf(data = world, fill = "grey85", color = "grey20", size = 0.1) +
    scale_fill_gradientn(colors = custom_colors, limits = c(0, 1), na.value = "transparent") +
    coord_sf() +
    labs(title = title, x = "Longitude", y = "Latitude") +
    theme_bw() +
    theme(
      legend.position = "none",  # Remove legend
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      axis.ticks = element_line(), # Ensure ticks are visible
      #panel.grid.major = element_line(color = "gray90", size = 0.2), # Light gridlines
      plot.margin = margin(5, 5, 5, 5), # Reduce extra space
      plot.annotate = element_text(size = 16, hjust = 0, vjust = 1, face = "bold"), # Position annotation text
      plot.title.position = "plot", # Title at the top of the plot
      plot.subtitle.position = "plot" # Subtitle at the top of the plot
    )
}

# Generate plots for each raster layer with formatted titles
plots <- list()

for (i in seq_along(layers)) {
  layer_name <- gsub("diff", "", layers[i]) # Remove 'diff' prefix
  formatted_title <- paste0("∆MAT = ", layer_name, "°C")
  
  plots[[paste0("p1_", layers[i])]] <- plot_raster(df_list[[layers[i]]]$pred_1w, formatted_title)
  plots[[paste0("p2_", layers[i])]] <- plot_raster(df_list[[layers[i]]]$pred_2w, formatted_title)
}

# Arrange plots in 2 columns with column titles
final_plot <- (plots$p1_diff5 + plots$p2_diff5) / 
  (plots$p1_diff10 + plots$p2_diff10) / 
  (plots$p1_diff15 + plots$p2_diff15) +
  plot_annotation(
    subtitle = c("1-week MHW", "2-week MHW"), # Column subtitles
    tag_levels = "a",  # Adds 'a' and 'b' subtitles for columns
    theme = theme(plot.subtitle = element_text(size = 14, hjust = 0.5, face = "bold"))
  ) 

# Show plot
print(final_plot)

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

# Extended data figures ----
##Extended Data Fig. 1 - historical vs predicted with abiotic factors ----
dieoffs_final  |> 
  mutate(new_severity = ifelse(severity == 0, NA,
                               as.character(abiotic_factor_s_severity))) |> 
  ggplot() + geom_abline() + 
  geom_point(aes(Survival, Predicted, shape = effect_on_survival,
                            fill = new_severity, colour = new_severity),
             size = 2) + xlim(0,1) + ylim(0,1) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=12,  family="Helvetica"),
        legend.position = c(0.85, 0.35)) +
  xlab("Observed survival") + ylab("Predicted survival") + 
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
                     

# Supplementary Information ----
##Supplementary Fig. 1 - data distribution ----
#SI Figure 1
p <- ggplot(data, aes(Temperature, log(Time), colour = Bioregion)) + 
  geom_point(size = 3) + theme_bw() + theme(panel.grid.major = element_blank(), 
                                            panel.grid.minor = element_blank(),
                                            text = element_text(size=12,  
                                                                family="Helvetica"),
                                            legend.position = c(0.2, 0.3),
                                            plot.subtitle = element_text(hjust = 0.01, vjust = -5, 
                                                                         margin = margin(2, 0, 0, 0), face = "bold")) + 
  ylab("Duration (Hours, log scale)") + xlab('Temperature (ºC)')
p <- ggMarginal(p, type = "histogram")
  

#SI Figure 2
s <- data  |> 
  dplyr::select(Genus, Species) |> 
  dplyr::mutate(name = paste0(substr(Genus, 1, 1), ". ", Species)) |> 
  dplyr::select(name)  |> 
  table()  |> 
  as.data.frame()
names(s) <- c('spp', 'freq')

s2 <- ggplot(s, aes(x = spp, y = freq)) + 
  geom_bar(stat = "identity") + theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(size = 12,  family="Helvetica"),
        axis.text.x = element_text(angle = 90, colour = 'black'),
        plot.subtitle = element_text(hjust = 0.01, vjust = -5, 
                                     margin = margin(2, 0, 0, 0), face = "bold")) +
  labs(subtitle = 'b') +
  xlab('') + ylab('Number of observations')

wrap_elements(p) + s2 + plot_layout(ncol = 1)

ggsave(file = paste0(FIGS_PATH, "/SF1.png"), 
       width = 600, 
       height = 600/1.6, 
       units = "mm", 
       dpi = 300)


##Supplementary Fig. 2 - sampling effort ----
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
