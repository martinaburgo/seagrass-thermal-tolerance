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
                           Temperature = Temperature,
                           difference_species = difference_species,
                           difference_population = difference_population,
                           av_population = av_population,
                           av_species = av_species,
                           mtwa_species = mtwa_species,
                           mtwa_population = mtwa_population,
                           Type = Type)), type = 'response') |>
    gather_emmeans_draws() |>
    mutate(.value = exp(.value)) |>
    summarise(median_hdci(.value)) |>
    as.data.frame() |>
    dplyr::select(y)
  
  print(i)
}

dieoffs_final <- dieoffs |>
  dplyr::distinct(ID, .keep_all = TRUE) |>
  full_join(dieoffs |>
              group_by(ID) |>
              summarise(Predicted = prod(partial_survival),
                        Duration = sum(Time))) |>
  mutate(Predicted = ifelse(Predicted > 1, 1, Predicted))
  
cor(dieoffs_final$Survival, dieoffs_final$Predicted)

dieoffs_final <- dieoffs_final |>
  add_column(readxl::read_xlsx('data/primary/seagrass-heatwaves-studies.xlsx', sheet = 2) |> 
               janitor::clean_names() |> 
               dplyr::distinct(event_id, .keep_all = TRUE) |> 
               dplyr::select(additional_abiotic_factor_s, abiotic_factor_s_severity, effect_on_survival)) |> 
  mutate(severity = ifelse(abiotic_factor_s_severity == 'Low', 1, 
                           ifelse(abiotic_factor_s_severity == 'Medium', 2,
                                  ifelse(abiotic_factor_s_severity == 'High', 3, 0)))) |> 
  mutate(severity = as.numeric(severity))

dieoffs_final |> 
  filter(!(effect_on_survival == "Negative" & severity == 3)) |> 
  droplevels() |> 
  ggplot(mapping = aes(Survival, Predicted)) + 
  geom_point(size = 2) + 
  geom_abline() + xlim(0, 1) + ylim(0, 1) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     text=element_text(size=12,  family="Helvetica")) +
  xlab("Observed survival") + ylab("Predicted survival")

write.csv(dieoffs_final, 'data/processed/heatwaves_predicted.csv')

trial <- dieoffs_final |> 
  filter(!(effect_on_survival == "Negative" & severity == 3)) |> 
  droplevels()
cor(trial$Survival, trial$Predicted)

trial |> 
  ggplot(mapping = aes(Survival, Predicted)) + 
  geom_point(size = 2) + 
  geom_abline() + xlim(0, 1) + ylim(0, 1) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     text=element_text(size=12,  family="Helvetica")) +
  xlab("Observed survival") + ylab("Predicted survival")
