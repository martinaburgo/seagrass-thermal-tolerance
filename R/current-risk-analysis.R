##LOAD PACKAGES ----
library(ncdf4)
library(raster)
library(brms)
library(tidyverse)

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

difference_population <- mat_population + 5
difference_10 <- mat_population + 10
difference_15 <- mat_population + 5

predict(brm8.form, stack(difference_population, av_population, mtwa_population), allow.new.levels = F)
