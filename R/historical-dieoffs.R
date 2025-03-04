##LOAD PACKAGES ----
library(brms)
library(tidyverse)
library(emmeans)
library(tidybayes)

source('R/functions.R')


#LOAD DATA ----
load(file = 'data/modelled/Model9Beta.RData') #model

dieoffs <- read.csv('data/processed/heatwaves_calculated.csv') |>
  dplyr::select(!(X))

# CALCULATE partial survival
dieoffs <- dieoffs |>
  mutate(partial_survival = NA)

for (i in 1:nrow(dieoffs)) {
  dieoffs[i, 'partial_survival'] <- seagrass.brm9 |> 
    emmeans(~ Time|difference_species|av_population|mtwa_population|Type, 
            at = with(dieoffs[i,],
                      list(Time = Time,
                           difference_species = difference_species,
                           av_population = av_population,
                           mtwa_population = mtwa_population,
                           Type = Type))) |>
    gather_emmeans_draws() |>
    mutate(.value = exp(.value)) |>
    summarise(median_hdci(.value)) |>
    as.data.frame() |>
    dplyr::select(ymin)
  
  print(i)
}

dieoffs_final <- dieoffs |>
  dplyr::distinct(ID, .keep_all = TRUE) |>
  full_join(dieoffs |>
              mutate(partial_survival = ifelse(partial_survival > 1, 1, partial_survival)) |>
              group_by(ID) |>
              summarise(Predicted = prod(partial_survival),
                        Duration = sum(Time)))
  
cor(dieoffs_final$Survival, dieoffs_final$Predicted)

dieoffs_final |>
  ggplot(mapping = aes(Survival, Predicted)) + geom_point(size = 3, 
                                                                   position = 'jitter') + 
  geom_abline() + xlim(0, 1) + ylim(0, 1) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     text=element_text(size=12,  family="Helvetica")) +
  xlab("Observed survial") + ylab("Predicted survival")
