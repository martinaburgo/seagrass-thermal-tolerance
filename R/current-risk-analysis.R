##LOAD PACKAGES ----
library(ncdf4)
library(raster)
library(brms)
library(tidyverse)
library(emmeans)
library(tidybayes)

source('R/functions.R')


#LOAD DATA ----
load(file = 'data/modelled/Model6Beta.RData')
load(file = 'data/modelled/Model7Beta.RData')
load(file = 'data/modelled/Model8Beta.RData')
seagrass <- read.csv(file = 'data/processed/seagrass_subset.csv')[, -1]

mat_population <- raster::stack('data/processed/mat.nc') |>
  calc(mean)
av_population <- raster::stack('data/processed/av.nc') |>
  calc(mean)
mtwa_population <- raster::stack('data/processed/mtwa.nc') |>
  calc(max)

# DATA CALC ----
## Diff 5 ----

difference_population <- mat_population + 5
enviro_5 <- raster::stack(difference_population, av_population, mtwa_population)
names(enviro_5) <- c('difference_population', 'av_population', 'mtwa_population')

enviro_5_sub <- enviro_5 |>
  as.data.frame() |>
  mutate(Survival = NA,
         ID = row_number()) |>
  dplyr::filter(!is.na(difference_population) | !is.na(av_population) | !is.na(mtwa_population))


for (i in 1:nrow(enviro_5_sub)) {
  diff <- enviro_5_sub$difference_population[i]
  av <- enviro_5_sub$av_population[i]
  mtwa <- enviro_5_sub$mtwa_population[i]
  
  enviro_5_sub$Survival[i] <- seagrass.brm6 |> 
    emmeans(~difference_population|av_population|mtwa_population, 
            at = list(difference_population = diff,
                      av_population = av,
                      mtwa_population = mtwa)) |>
    gather_emmeans_draws() |>
    mutate(.value = plogis(.value)) |>
    summarise(median_hdci(.value)) |>
    as.data.frame() |>
    dplyr::select(y)
}


write_csv(enviro_5_sub, file = 'data/processed/enviro_5.csv')


## Diff 10 ----

difference_population <- mat_population + 10
enviro_10 <- raster::stack(difference_population, av_population, mtwa_population)
names(enviro_10) <- c('difference_population', 'av_population', 'mtwa_population')

enviro_10_sub <- enviro_10 |>
  as.data.frame() |>
  mutate(Survival = NA,
         ID = row_number()) |>
  dplyr::filter(!is.na(difference_population) | !is.na(av_population) | !is.na(mtwa_population))


for (i in 1:nrow(enviro_10_sub)) {
  diff <- enviro_10_sub$difference_population[i]
  av <- enviro_10_sub$av_population[i]
  mtwa <- enviro_10_sub$mtwa_population[i]
  
  enviro_10_sub$Survival[i] <- seagrass.brm6 |> 
    emmeans(~difference_population|av_population|mtwa_population, 
            at = list(difference_population = diff,
                      av_population = av,
                      mtwa_population = mtwa)) |>
    gather_emmeans_draws() |>
    mutate(.value = plogis(.value)) |>
    summarise(median_hdci(.value)) |>
    as.data.frame() |>
    dplyr::select(y)
}

write_csv(enviro_10_sub, file = 'data/processed/enviro_10.csv')

## Diff 15 ----

difference_population <- mat_population + 15
enviro_15 <- raster::stack(difference_population, av_population, mtwa_population)
names(enviro_15) <- c('difference_population', 'av_population', 'mtwa_population')

enviro_15_sub <- enviro_15 |>
  as.data.frame() |>
  mutate(Survival = NA,
         ID = row_number()) |>
  dplyr::filter(!is.na(difference_population) | !is.na(av_population) | !is.na(mtwa_population))


for (i in 1:nrow(enviro_15_sub)) {
  diff <- enviro_15_sub$difference_population[i]
  av <- enviro_15_sub$av_population[i]
  mtwa <- enviro_15_sub$mtwa_population[i]
  
  enviro_15_sub$Survival[i] <- seagrass.brm6 |> 
    emmeans(~difference_population|av_population|mtwa_population, 
            at = list(difference_population = diff,
                      av_population = av,
                      mtwa_population = mtwa)) |>
    gather_emmeans_draws() |>
    mutate(.value = plogis(.value)) |>
    summarise(median_hdci(.value)) |>
    as.data.frame() |>
    dplyr::select(y)
}

write_csv(enviro_15_sub, file = 'data/processed/enviro_15.csv')


## Create survival raster stack ----
survival <- raster::raster(nrow = nrow(enviro_5[[1]]),
                           ncol = ncol(enviro_5[[1]]),
                           ext = extent(enviro_5[[1]]),
                           crs = crs(enviro_5[[1]]))
demo <- enviro_5_sub
row.names(demo) <- demo$ID

survival@data $ surv_5 <- merge(enviro_5 |> 
        as.data.frame(), demo |>
        dplyr::select(!ID), 
      by = 'row.names', all = TRUE) |>
  dplyr::select(Survival)

