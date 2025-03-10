##LOAD PACKAGES ----
library(raster)
library(brms)
library(tidyverse)
library(emmeans)
library(tidybayes)
library(sf)

source('R/functions.R')


#LOAD DATA ----
#model
load(file = 'data/modelled/Model7Beta.RData')

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

# climate data
mat <- raster::stack("data/processed/mat.nc") |>
  mean()
av <- raster::stack("data/processed/av.nc") |>
  mean()
mtwa <- raster::stack("data/processed/mtwa.nc") |>
  max()

# DATA CALC ----
## Create df to store values
env_seagrasses_df <- data.frame(binomial = seagrasses$binomial,
                                av_species = NA,
                                mat_species = NA,
                                mtwa_species = NA,
                                Type = NA)

# extract climate variables
for (i in 1:nrow(env_seagrasses_df)) {
  env_seagrasses_df[i, 'mat_species'] <- crop(mat, 
                                              seagrasses[seagrasses$binomial == env_seagrasses_df[i, 'binomial'],]) |>
    values() |> 
    na.omit() |> 
    mean()
  
  env_seagrasses_df[i, 'av_species'] <- crop(av,
                                             seagrasses[seagrasses$binomial == env_seagrasses_df[i, 'binomial'],]) |>
    values() |> 
    na.omit() |> 
    mean()
  
  env_seagrasses_df[i, 'mtwa_species'] <- crop(mtwa,
                                               seagrasses[seagrasses$binomial == env_seagrasses_df[i, 'binomial'],]) |>
    values() |> 
    na.omit() |> 
    max()
  
  env_seagrasses_df[i, 'Type'] <- seagrasses[seagrasses$binomial == env_seagrasses_df[i, 'binomial'], ]$type
  
  print(i)
}

# PREDICT ----
## ONE WEEK at 5, 10, and 15 degree above MAT ----
pred_1w <- env_seagrasses_df |>
  mutate(diff5 = NA,
         diff10 = NA,
         diff15 = NA)

for (i in 1:nrow(pred_1w)) {
  pred_1w[i, 'diff5'] <- seagrass.brm7 |> 
    emmeans(~ Time|difference_species|av_species|mtwa_species|Type, 
            at = with(pred_1w[i,],
                      list(Time = 24*7,
                           difference_species = mat_species + 5,
                           av_species = av_species,
                           mtwa_species = mtwa_species,
                           Type = Type))) |>
    gather_emmeans_draws() |>
    mutate(.value = exp(.value)) |>
    summarise(median_hdci(.value)) |>
    as.data.frame() |>
    dplyr::select(y)
  
  pred_1w[i, 'diff10'] <- seagrass.brm7 |> 
    emmeans(~ Time|difference_species|av_species|mtwa_species|Type, 
            at = with(pred_1w[i,],
                      list(Time = 24*7,
                           difference_species = mat_species + 10,
                           av_species = av_species,
                           mtwa_species = mtwa_species,
                           Type = Type))) |>
    gather_emmeans_draws() |>
    mutate(.value = exp(.value)) |>
    summarise(median_hdci(.value)) |>
    as.data.frame() |>
    dplyr::select(y)
  
  pred_1w[i, 'diff15'] <- seagrass.brm7 |> 
    emmeans(~ Time|difference_species|av_species|mtwa_species|Type, 
            at = with(pred_1w[i,],
                      list(Time = 24*7,
                           difference_species = mat_species + 15,
                           av_species = av_species,
                           mtwa_species = mtwa_species,
                           Type = Type))) |>
    gather_emmeans_draws() |>
    mutate(.value = exp(.value)) |>
    summarise(median_hdci(.value)) |>
    as.data.frame() |>
    dplyr::select(y)
  
  print(i)
}

write.csv(pred_1w, 'data/processed/pred_1w.csv')

## TWO WEEKS at 5, 10, and 15 degree above MAT ----
pred_2w <- env_seagrasses_df |>
  mutate(diff5 = NA,
         diff10 = NA,
         diff15 = NA)

for (i in 1:nrow(pred_2w)) {
  pred_2w[i, 'diff5'] <- seagrass.brm7 |> 
    emmeans(~ Time|difference_species|av_species|mtwa_species|Type, 
            at = with(pred_2w[i,],
                      list(Time = 24*14,
                           difference_species = mat_species + 5,
                           av_species = av_species,
                           mtwa_species = mtwa_species,
                           Type = Type))) |>
    gather_emmeans_draws() |>
    mutate(.value = exp(.value)) |>
    summarise(median_hdci(.value)) |>
    as.data.frame() |>
    dplyr::select(y)
  
  pred_2w[i, 'diff10'] <- seagrass.brm7 |> 
    emmeans(~ Time|difference_species|av_species|mtwa_species|Type, 
            at = with(pred_2w[i,],
                      list(Time = 24*14,
                           difference_species = mat_species + 10,
                           av_species = av_species,
                           mtwa_species = mtwa_species,
                           Type = Type))) |>
    gather_emmeans_draws() |>
    mutate(.value = exp(.value)) |>
    summarise(median_hdci(.value)) |>
    as.data.frame() |>
    dplyr::select(y)
  
  pred_2w[i, 'diff15'] <- seagrass.brm7 |> 
    emmeans(~ Time|difference_species|av_species|mtwa_species|Type, 
            at = with(pred_2w[i,],
                      list(Time = 24*14,
                           difference_species = mat_species + 15,
                           av_species = av_species,
                           mtwa_species = mtwa_species,
                           Type = Type))) |>
    gather_emmeans_draws() |>
    mutate(.value = exp(.value)) |>
    summarise(median_hdci(.value)) |>
    as.data.frame() |>
    dplyr::select(y)
  
  print(i)
}

write.csv(pred_2w, 'data/processed/pred_2w.csv')

# STACKS for PLOTs ----
## ONE WEEK
pred_1w_diff5 <- mat

for (i in 1:nrow(pred_1w)) {
  int_raster <- mat
  values(int_raster) <- pred_1w[i, 'diff5']
  pred_1w_diff5 <- mask(int_raster, seagrasses[seagrasses$binomial == pred_1w[i, 'binomial'],]) |>
    addLayer(pred_1w_diff5)
  print(i)
}
pred_1w_diff5 <- pred_1w_diff5[[-nlayers(pred_1w_diff5)]]

pred_1w_diff10 <- mat

for (i in 1:nrow(pred_1w)) {
  int_raster <- mat
  values(int_raster) <- pred_1w[i, 'diff10']
  pred_1w_diff10 <- mask(int_raster, seagrasses[seagrasses$binomial == pred_1w[i, 'binomial'],]) |>
    addLayer(pred_1w_diff10)
  print(i)
}
pred_1w_diff10 <- pred_1w_diff10[[-nlayers(pred_1w_diff10)]]


pred_1w_diff15 <- mat

for (i in 1:nrow(pred_1w)) {
  int_raster <- mat
  values(int_raster) <- pred_1w[i, 'diff15']
  pred_1w_diff15 <- mask(int_raster, seagrasses[seagrasses$binomial == pred_1w[i, 'binomial'],]) |>
    addLayer(pred_1w_diff15)
  print(i)
}
pred_1w_diff15 <- pred_1w_diff15[[-nlayers(pred_1w_diff15)]]

# Final stack
pred_1w_stack <- addLayer(stackApply(pred_1w_diff5,
                    indices = c(rep(1, nlayers(pred_1w_diff5))),
                    fun = na.omit(mean)),
         stackApply(pred_1w_diff10,
                    indices = c(rep(1, nlayers(pred_1w_diff10))),
                    fun = na.omit(mean)),
         stackApply(pred_1w_diff15,
                    indices = c(rep(1, nlayers(pred_1w_diff15))),
                    fun = na.omit(mean)))
names(pred_1w_stack) <- c('diff5', 'diff10', 'diff15')

writeRaster(pred_1w_stack, filename = 'data/processed/pred_1w_stack.grd', overwrite = TRUE)


## TWO WEEKS
pred_2w_diff5 <- mat

for (i in 1:nrow(pred_2w)) {
  int_raster <- mat
  values(int_raster) <- pred_2w[i, 'diff5']
  pred_2w_diff5 <- mask(int_raster, seagrasses[seagrasses$binomial == pred_2w[i, 'binomial'],]) |>
    addLayer(pred_2w_diff5)
  print(i)
}
pred_2w_diff5 <- pred_2w_diff5[[-nlayers(pred_2w_diff5)]]

pred_2w_diff10 <- mat

for (i in 1:nrow(pred_2w)) {
  int_raster <- mat
  values(int_raster) <- pred_2w[i, 'diff10']
  pred_2w_diff10 <- mask(int_raster, seagrasses[seagrasses$binomial == pred_2w[i, 'binomial'],]) |>
    addLayer(pred_2w_diff10)
  print(i)
}
pred_2w_diff10 <- pred_2w_diff10[[-nlayers(pred_2w_diff10)]]


pred_2w_diff15 <- mat

for (i in 1:nrow(pred_2w)) {
  int_raster <- mat
  values(int_raster) <- pred_2w[i, 'diff15']
  pred_2w_diff15 <- mask(int_raster, seagrasses[seagrasses$binomial == pred_2w[i, 'binomial'],]) |>
    addLayer(pred_2w_diff15)
  print(i)
}
pred_2w_diff15 <- pred_2w_diff15[[-nlayers(pred_2w_diff15)]]

# Final stack
pred_2w_stack <- addLayer(stackApply(pred_2w_diff5,
                                     indices = c(rep(1, nlayers(pred_2w_diff5))),
                                     fun = na.omit(mean)),
                          stackApply(pred_2w_diff10,
                                     indices = c(rep(1, nlayers(pred_2w_diff10))),
                                     fun = na.omit(mean)),
                          stackApply(pred_2w_diff15,
                                     indices = c(rep(1, nlayers(pred_2w_diff15))),
                                     fun = na.omit(mean)))
names(pred_2w_stack) <- c('diff5', 'diff10', 'diff15')

writeRaster(pred_2w_stack, filename = 'data/processed/pred_2w_stack.grd', overwrite = TRUE)

# Estimates ----
raster_data <- pred_1w_stack$diff5

## Tropics vs Temperate ----
# Define extents
tropics_extent <- extent(-180, 180, -23.5, 23.5)  # Between -23.5 and 23.5
outside_tropics_south_extent <- extent(-180, 180, -90, -23.5)  # South of -23.5
outside_tropics_north_extent <- extent(-180, 180, 23.5, 90)  # North of 23.5

# Convert extents to SpatialPolygons and rasterize
mask_tropics <- rasterize(as(tropics_extent, "SpatialPolygons"), raster_data)
mask_outside_south <- rasterize(as(outside_tropics_south_extent, "SpatialPolygons"), raster_data)
mask_outside_north <- rasterize(as(outside_tropics_north_extent, "SpatialPolygons"), raster_data)

# Mask the raster to extract values
raster_tropics <- mask(raster_data, mask_tropics)
raster_outside_south <- mask(raster_data, mask_outside_south)
raster_outside_north <- mask(raster_data, mask_outside_north)

# Combine the two "outside tropics" regions
raster_outside_tropics <- mosaic(raster_outside_south, raster_outside_north, fun = "mean")

# Compute mean and standard deviation
mean_tropics <- cellStats(raster_tropics, stat = 'mean')
sd_tropics <- cellStats(raster_tropics, stat = 'sd')
n_tropics <- sum(!is.na(values(raster_tropics)))
sem_tropics <- sd_tropics / sqrt(n_tropics)  # Standard Error of the Mean

mean_outside_tropics <- cellStats(raster_outside_tropics, stat = 'mean')
sd_outside_tropics <- cellStats(raster_outside_tropics, stat = 'sd')
n_outside_tropics <- sum(!is.na(values(raster_outside_tropics)))
sem_outside_tropics <- sd_outside_tropics / sqrt(n_outside_tropics)  # Standard Error of the Mean

# Print results
cat("Tropics (-23.5 to 23.5) - Mean:", mean_tropics, " SEM:", sem_tropics, "\n")
cat("Outside Tropics (< -23.5 or > 23.5) - Mean:", mean_outside_tropics, " SEM:", sem_outside_tropics, "\n")

## Northern vs Southern ----
# Define extents for Northern and Southern Hemisphere
north_extent <- extent(-180, 180, 0, 90)
south_extent <- extent(-180, 180, -90, 0)

# Convert extents to SpatialPolygons and rasterize them
mask_north <- rasterize(as(north_extent, "SpatialPolygons"), raster_data)
mask_south <- rasterize(as(south_extent, "SpatialPolygons"), raster_data)

# Mask the raster to extract values for each hemisphere
raster_north <- mask(raster_data, mask_north)
raster_south <- mask(raster_data, mask_south)

# Compute mean and standard deviation for both hemispheres
mean_north <- cellStats(raster_north, stat = 'mean')
sd_north <- cellStats(raster_north, stat = 'sd')
n_north <- sum(!is.na(values(raster_north)))  # Count valid (non-NA) cells
sem_north <- sd_north / sqrt(n_north)  # Standard Error of the Mean

mean_south <- cellStats(raster_south, stat = 'mean')
sd_south <- cellStats(raster_south, stat = 'sd')
n_south <- sum(!is.na(values(raster_south)))  # Count valid (non-NA) cells
sem_south <- sd_south / sqrt(n_south)  # Standard Error of the Mean

# Print results
cat("Northern Hemisphere (lat > 0) - Mean:", mean_north, " SEM:", sem_north, "\n")
cat("Southern Hemisphere (lat < 0) - Mean:", mean_south, " SEM:", sem_south, "\n")

