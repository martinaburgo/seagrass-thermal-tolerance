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

writeRaster(pred_1w_stack, filename = 'data/processed/pred_1w_stack.grd')


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

writeRaster(pred_2w_stack, filename = 'data/processed/pred_2w_stack.grd')
