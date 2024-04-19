# Preparations ----

## Load the necessary libraries ----
library(tidyverse)   #for data wrangling
library(stringr)
library(car)
library(brms)
library(patchwork)
library(corrplot)
library(loo)
library(tidybayes)
library(DHARMa)     #for residual diagnostics
library(rstan)      #for interfacing with STAN
library(emmeans)

source('R/functions.R')

## Read in the data ----
seagrass <- read_csv('data/processed/data_calculatedv2.csv', trim_ws = TRUE)[, -1]
seagrass |> str()

# Declare categorical variables and create a column with the scientific name of each seagrass:
seagrass <- seagrass |>
  mutate(Study = factor(Study),
         Location = factor(Location),
         Bioregion = factor(Bioregion),
         Family = factor(Family),
         Genus = factor(Genus),
         Species = factor(Species),
         Type = factor(Type),
         Variable = factor(Variable),
         Name = factor(paste(stringr::str_extract(Genus, "^.{1}"), '.', Species, sep = '')))

# Exploratory data analysis ----
seagrass |>
  ggplot(aes(y = Survival, x = Temperature)) +
  geom_point() + theme_classic()

# Some studies also looked at the lower thermal limit. Since we're interested in the upper thermal limit, any study with experimental temperature 
#lower than MAT will be excluded.Some studies exposed seagrasses to hight heat stress for very short periods of time (< 1 hour) and these studies 
#were excluded as they do not represent realistic in-situ conditions.
seagrass <- seagrass |> 
  filter(Temperature >= mat_population) |>
  filter(Time >= 1)
write.csv(seagrass, file = 'data/processed/seagrass_subset.csv') #saving the subset data used for analysis

seagrass |> ggplot(aes(x = log10(Time), y = Temperature, color = Survival)) +
  geom_point(size = 3) + scale_x_continuous('Time (log10)') + 
  scale_color_viridis_c(option = 'C') + theme_classic()

seagrass |>
  ggplot(aes(y = Survival, x = difference_species, color = av_species)) +
  geom_point(size = 3) + theme_classic()

#And plotting the potential predictors to check if they're correlated:
scatterplotMatrix(~Survival+log(Time)+Temperature+mat_population+av_population+mtwa_population+difference_population + magnitude_population + av_species + 
                    mat_species + mtwa_species + difference_species + magnitude_species, 
                  data = seagrass,diagonal = list(method = 'boxplot'))

seagrass |> dplyr::select(Time, Temperature, mat_population, av_population, mtwa_population, difference_population, magnitude_population, av_species, 
                          mat_species, mtwa_species, difference_species, magnitude_species) |> cor() |> corrplot()

#As expected, most climate variables are correlated. To avoid problems associated with collinearity, magnitude will be used in the models to 
#understand the effect of environmental conditions on seagrass thermal tolerance.

#There's two ways we could analyse data:

#1.  Binomial distribution: binning survival in either 1 (survival \> 0.5) or 0 (survival \<= 0)

#  Using a beta distribution (zero and one inflated)


# BINOMIAL DISTRIBUTION ----
# First, we'll adjust the survival to fit a binomial distribution:
seagrass <- seagrass |>
  mutate(surv_adj = ifelse(Survival > 0.5, 1, 0))

## Model 1: effect of temperature and duration of heat stress ----

### Fit model ----
#The model will be run using 'Study' as varying effect. This is to account for the variation methodology between study.

brm1.form <- bf(surv_adj|trials(1) ~ scale(Temperature)*scale(log(Time)) + (1|Study), family = binomial(link = 'logit')) 

brm1.form |> get_prior(data = seagrass) # to check which priors we need

priors1 <- prior(normal(0, 2.5), class = 'Intercept') +
  prior(normal(0, 1), class = 'b') + 
  prior(student_t(3, 0, 2.5), class = 'sd')


seagrass.brm1 <- brm(brm1.form, prior = priors1, data = seagrass, 
                     sample_prior = 'yes', 
                     iter = 5000, 
                     warmup = 1000, 
                     chains = 3, cores = 3, 
                     thin = 5, 
                     refresh = 0, 
                     backend = 'rstan') 

seagrass.brm1 |> SUYR_prior_and_posterior()

### Diagnostics ----
(seagrass.brm1$fit |> stan_trace()) + (seagrass.brm1$fit |> stan_ac()) + (seagrass.brm1$fit |> stan_rhat()) + (seagrass.brm1$fit |> stan_ess())

seagrass.brm1 |> pp_check(type = 'dens_overlay', ndraws = 100)

seagrass.resids <- make_brms_dharma_res(seagrass.brm1, integerResponse = FALSE)
testUniformity(seagrass.resids)
plotResiduals(seagrass.resids, form = factor(rep(1, nrow(seagrass))))
plotResiduals(seagrass.resids, quantreg = FALSE)
testDispersion(seagrass.resids)

save(seagrass.brm1, seagrass, priors1, brm1.form, file = '../data/modelled/Model1Binomial.RData')

### Model investigation ----
seagrass.brm1 |> 
  conditional_effects(spaghetti = TRUE, ndraws = 100) |>
  plot(points = TRUE, ask = FALSE)
#to help interpret the model

seagrass.brm1 |>
  as.tibble() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1)) |>
  as.tibble() |>
  dplyr::slice(1:5)


## Model 2: effect of magnitude and time on survival (pop) ----

### Fit model ----
brm2.form <- bf(surv_adj|trials(1) ~ scale(magnitude_population)*scale(log(Time)) + mtwa_population + (1|Study), family = binomial(link = 'logit')) 

brm2.form |> get_prior(data = seagrass) # to check which priors we need

priors2 <- prior(normal(0, 2.5), class = 'Intercept') +
  prior(normal(0, 10), class = 'b') + 
  prior(student_t(3, 0, 2.5), class = 'sd')


seagrass.brm2 <- brm(brm2.form, prior = priors2, data = seagrass, 
                     sample_prior = 'yes', 
                     iter = 5000, 
                     warmup = 1000, 
                     chains = 3, cores = 3, 
                     thin = 5, 
                     refresh = 0, 
                     backend = 'rstan') 

seagrass.brm2 |> SUYR_prior_and_posterior()

### Diagnostics ----
(seagrass.brm2$fit |> stan_trace()) + (seagrass.brm2$fit |> stan_ac()) + (seagrass.brm2$fit |> stan_rhat()) + (seagrass.brm2$fit |> stan_ess())

seagrass.brm2 |> pp_check(type = 'dens_overlay', ndraws = 100)

seagrass.resids <- make_brms_dharma_res(seagrass.brm2, integerResponse = FALSE) 
testUniformity(seagrass.resids)
plotResiduals(seagrass.resids, form = factor(rep(1, nrow(seagrass))))
plotResiduals(seagrass.resids, quantreg = TRUE)
testDispersion(seagrass.resids)

save(seagrass.brm2, seagrass, priors2, brm2.form, file = '../data/modelled/Model2Binomial.RData')

### Model investigation ----
seagrass.brm2 |> 
  conditional_effects(spaghetti = TRUE, ndraws = 100) |>
  plot(points = TRUE, ask = FALSE)
#to help interpret the model

seagrass.brm2 |>
  as.tibble() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1)) |>
  as.tibble() |>
  dplyr::slice(1:6)

## Model 3: effect of magnitude and time on survival (spp) ----

### Fit model ----

brm3.form <- bf(surv_adj|trials(1) ~ scale(magnitude_species)*scale(log(Time)) + mtwa_species + (1|Study), family = binomial(link = 'logit')) 

brm3.form |> get_prior(data = seagrass) # to check which priors we need

priors3 <- prior(normal(0, 2.5), class = 'Intercept') +
  prior(normal(0, 10), class = 'b') + 
  prior(student_t(3, 0, 2.5), class = 'sd')


seagrass.brm3 <- brm(brm3.form, prior = priors3, data = seagrass, 
                     sample_prior = 'yes', 
                     iter = 5000, 
                     warmup = 1000, 
                     chains = 3, cores = 3, 
                     thin = 5, 
                     refresh = 0, 
                     backend = 'rstan') 

seagrass.brm3 |> SUYR_prior_and_posterior()


### Diagnostics ----
(seagrass.brm3$fit |> stan_trace()) + (seagrass.brm3$fit |> stan_ac()) + (seagrass.brm3$fit |> stan_rhat()) + (seagrass.brm3$fit |> stan_ess())

seagrass.brm3 |> pp_check(type = 'dens_overlay', ndraws = 100)

seagrass.resids <- make_brms_dharma_res(seagrass.brm3, integerResponse = FALSE) 
testUniformity(seagrass.resids)
plotResiduals(seagrass.resids, form = factor(rep(1, nrow(seagrass))))
plotResiduals(seagrass.resids, quantreg = TRUE)
testDispersion(seagrass.resids)

save(seagrass.brm3, seagrass, priors3, brm3.form, file = '../data/modelled/Model3Binomial.RData')

### Model investigation ----
seagrass.brm3 |> 
  conditional_effects(spaghetti = TRUE, ndraws = 100) |>
  plot(points = TRUE, ask = FALSE)

seagrass.brm3 |>
  as.tibble() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1)) |>
  as.tibble() |>
  dplyr::slice(1:6)

## Compare ----
loo::loo_compare(rstan::loo(seagrass.brm1),
                 rstan::loo(seagrass.brm2),
                 rstan::loo(seagrass.brm3))


## Model 4: effect of local climate on survival ----

### Fit model ----
brm4.form <- bf(surv_adj|trials(1) ~ scale(log(Time)) + scale(difference_population) + scale(av_population) + scale(mtwa_population) + (1|Study), family = binomial(link = 'logit')) 

brm4.form |> get_prior(data = seagrass) # to check which priors we need

priors4 <- prior(normal(0, 2.5), class = 'Intercept') +
  prior(normal(0, 10), class = 'b') + 
  prior(student_t(3, 0, 2.5), class = 'sd')

seagrass.brm4 <- brm(brm4.form, prior = priors4, data = seagrass, 
                     sample_prior = 'yes', 
                     iter = 5000, 
                     warmup = 1000, 
                     chains = 3, cores = 3, 
                     thin = 5, 
                     refresh = 0, 
                     backend = 'rstan') 

seagrass.brm4 |> SUYR_prior_and_posterior()


### Diagnostics ----
(seagrass.brm4$fit |> stan_trace()) + (seagrass.brm4$fit |> stan_ac()) + (seagrass.brm4$fit |> stan_rhat()) + (seagrass.brm4$fit |> stan_ess())

seagrass.brm4 |> pp_check(type = 'dens_overlay', ndraws = 100)

seagrass.resids <- make_brms_dharma_res(seagrass.brm4, integerResponse = FALSE) 
testUniformity(seagrass.resids)
plotResiduals(seagrass.resids, form = factor(rep(1, nrow(seagrass))))
plotResiduals(seagrass.resids, quantreg = TRUE)
testDispersion(seagrass.resids)


save(seagrass.brm4, seagrass, priors4, brm4.form, file = '../data/modelled/Model4Binomial.RData')

### Model investigation ----
seagrass.brm4 |> 
  conditional_effects(spaghetti = TRUE, ndraws = 100) |>
  plot(points = TRUE, ask = FALSE)

seagrass.brm4 |>
  as.tibble() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1)) |>
  as.tibble() |>
  dplyr::slice(1:7)

## Model 5: effect of species climate on survival ----

### Fit model ----
brm5.form <- bf(surv_adj|trials(1) ~ scale(log(Time)) + scale(difference_species) + scale(av_species) + scale(mtwa_species) + (1|Study), family = binomial(link = 'logit')) 

brm5.form |> get_prior(data = seagrass) # to check which priors we need

priors5 <- prior(normal(0, 2.5), class = 'Intercept') +
  prior(normal(0, 10), class = 'b') + 
  prior(student_t(3, 0, 2.5), class = 'sd')

seagrass.brm5 <- brm(brm5.form, prior = priors5, data = seagrass, 
                     sample_prior = 'yes', 
                     iter = 5000, 
                     warmup = 1000, 
                     chains = 3, cores = 3, 
                     thin = 5, 
                     refresh = 0, 
                     backend = 'rstan') 

seagrass.brm5 |> SUYR_prior_and_posterior()

### Diagnostics ----
(seagrass.brm5$fit |> stan_trace()) + (seagrass.brm5$fit |> stan_ac()) + (seagrass.brm5$fit |> stan_rhat()) + (seagrass.brm5$fit |> stan_ess())

seagrass.brm5 |> pp_check(type = 'dens_overlay', ndraws = 100)

seagrass.resids <- make_brms_dharma_res(seagrass.brm5, integerResponse = FALSE) 
testUniformity(seagrass.resids)
plotResiduals(seagrass.resids, form = factor(rep(1, nrow(seagrass))))
plotResiduals(seagrass.resids, quantreg = FALSE)
testDispersion(seagrass.resids)

save(seagrass.brm5, seagrass, priors5, brm5.form, file = '../data/modelled/Model5Binomial.RData')

### Model investigation ----
seagrass.brm5 |> 
  conditional_effects(spaghetti = TRUE, ndraws = 100) |>
  plot(points = TRUE, ask = FALSE)

seagrass.brm5 |>
  as.tibble() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median,
                  HDInterval::hdi,
                  rhat, length, ess_bulk, ess_tail,
                  Pl = ~mean(.x < 1),
                  Pg = ~mean(.x > 1)) |>
  as.tibble() |>
  dplyr::slice(1:6)

loo::loo_compare(rstan::loo(seagrass.brm1),
                 rstan::loo(seagrass.brm4),
                 rstan::loo(seagrass.brm5))


# BETA DISTRIBUTION ----
## Model 6: effect of local climate on survival ----

### Fit the model ----
brm6.form <- bf(Survival ~ scale(log(Time)) * scale(difference_population) + scale(av_population) + scale(mtwa_population) + (1|Study),
                family = zero_one_inflated_beta())

get_prior(brm6.form, data = seagrass)

seagrass |> 
  summarise(logit(median(Survival)),
            logit(mad(Survival)))

priors6 <- prior(normal(-2, 2), class = 'Intercept') +
  prior(normal(0, 5), class = 'b') + 
  prior(student_t(3, 0,5), class = 'sd')  +
  prior(gamma(0.01, 0.01), class = 'phi') +
  prior(beta(1, 1), class = 'zoi') +
  prior(beta(1, 1), class = 'coi')

seagrass.brm6p <- brm(brm6.form, prior = priors6, data = seagrass, 
                      sample_prior = 'only', 
                      iter = 5000, 
                      warmup = 1000, 
                      chains = 3, cores = 3, 
                      thin = 5, 
                      control = list(adapt_delta = 0.99, max_treedepth = 20),
                      refresh = 100, 
                      backend = 'rstan') 

seagrass.brm6p |> conditional_effects() |> plot(points = TRUE, ask = FALSE)

seagrass.brm6 <- update(seagrass.brm6p, sample_prior = 'yes') # add data

seagrass.brm6 |> SUYR_prior_and_posterior()

### Diagnostics ----
seagrass.brm6$fit |> stan_trace()
seagrass.brm6$fit |> stan_ac()
seagrass.brm6$fit |> stan_rhat()
seagrass.brm6$fit |> stan_ess()

seagrass.brm6 |> pp_check(type = 'dens_overlay', ndraws = 100)

seagrass.resids <- make_brms_dharma_res(seagrass.brm6, integerResponse = FALSE)
wrap_elements(~testUniformity(seagrass.resids)) +
  wrap_elements(~plotResiduals(seagrass.resids, form = factor(rep(1, nrow(seagrass))))) + 
  wrap_elements(~plotResiduals(seagrass.resids, quantreg = TRUE)) + 
  wrap_elements(~testDispersion(seagrass.resids)) 

save(seagrass.brm6, seagrass, priors6, brm6.form, file = '../data/modelled/Model6Beta.RData')

### Model investigation ----
seagrass.brm6 |> 
  conditional_effects(spaghetti = TRUE, ndraws = 100) |>
  plot(points = TRUE, ask = FALSE)

seagrass.brm6 |>
  as_draws_df() |>
  mutate(across(everything(), exp)) |> #to go from log scale to odds ratio
  summarise_draws(median, HDInterval::hdi, rhat, length, ess_bulk, ess_tail,
                  Pl = ~ mean(.x < 1),
                  Pg = ~ mean(.x > 1)) |>
  as_tibble() |>
  dplyr::slice(1:7)



## Model 7: effect of species-level climate on survival ----

### Fit the model ----
brm7.form <- bf(Survival ~ scale(log(Time)) * scale(difference_species) + scale(av_species) + scale(mtwa_species) + (1|Study),
                family = zero_one_inflated_beta())

get_prior(brm7.form, data = seagrass)

seagrass |> 
  summarise(logit(median(Survival)),
            logit(mad(Survival)))

priors7 <- prior(normal(-2, 2), class = 'Intercept') +
  prior(normal(0, 5), class = 'b') + 
  prior(student_t(3, 0,5), class = 'sd')  +
  prior(gamma(0.01, 0.01), class = 'phi') +
  prior(beta(1, 1), class = 'zoi') +
  prior(beta(1, 1), class = 'coi')

seagrass.brm7p <- brm(brm7.form, prior = priors7, data = seagrass, 
                      sample_prior = 'only', 
                      iter = 5000, 
                      warmup = 1000, 
                      chains = 3, cores = 3, 
                      thin = 5, 
                      control = list(adapt_delta = 0.99, max_treedepth = 50),
                      refresh = 100, 
                      backend = 'rstan') 

seagrass.brm7p |> conditional_effects() |> plot(points = TRUE, ask = FALSE)

seagrass.brm7 <- update(seagrass.brm7p, sample_prior = 'yes')

seagrass.brm7 |> SUYR_prior_and_posterior()

### Diagnostics ----
seagrass.brm7$fit |> stan_trace()
seagrass.brm7$fit |> stan_ac()
seagrass.brm7$fit |> stan_rhat()
seagrass.brm7$fit |> stan_ess()

seagrass.brm7 |> pp_check(type = 'dens_overlay', ndraws = 100) 

seagrass.resids <- make_brms_dharma_res(seagrass.brm7, integerResponse = FALSE)
wrap_elements(~testUniformity(seagrass.resids)) +
  wrap_elements(~plotResiduals(seagrass.resids, form = factor(rep(1, nrow(seagrass))))) + 
  wrap_elements(~plotResiduals(seagrass.resids, quantreg = TRUE)) + 
  wrap_elements(~testDispersion(seagrass.resids)) 

save(seagrass.brm7, seagrass, priors7, brm7.form, file = '../data/modelled/Model7Beta.RData')

### Model investigation ----
seagrass.brm7 |> 
  conditional_effects(spaghetti = TRUE, ndraws = 100) |>
  plot(points = TRUE, ask = FALSE)

seagrass.brm7 |>
  as_draws_df() |>
  mutate(across(everything(), exp)) |> #to go from log scale to odds ratio
  summarise_draws(median, HDInterval::hdi, rhat, length, ess_bulk, ess_tail,
                  Pl = ~ mean(.x < 1),
                  Pg = ~ mean(.x > 1)) |>
  as_tibble() |>
  dplyr::slice(1:7)


## Model 8: effect of temp and type ----

### Fit the model ----
brm8.form <- bf(Survival ~ scale(log(Time)) * scale(Temperature) + Type + (1|Study),
                family = zero_one_inflated_beta())

get_prior(brm8.form, data = seagrass)

seagrass |> 
  group_by(Type) |>
  summarise(logit(median(Survival)),
            logit(mad(Survival)))

priors8 <- prior(normal(-4, 4), class = 'Intercept') +
  prior(normal(0, 5), class = 'b') + 
  prior(student_t(3, 0,5), class = 'sd')  +
  prior(gamma(0.01, 0.01), class = 'phi') +
  prior(beta(1, 1), class = 'zoi') +
  prior(beta(1, 1), class = 'coi')

seagrass.brm8p <- brm(brm8.form, prior = priors8, data = seagrass, 
                      sample_prior = 'only', 
                      iter = 5000, 
                      warmup = 1000, 
                      chains = 3, cores = 3, 
                      thin = 5, 
                      control = list(adapt_delta = 0.99, max_treedepth = 20),
                      refresh = 100, 
                      backend = 'rstan') 


seagrass.brm8p |> conditional_effects() |> plot(points = TRUE, ask = FALSE)

seagrass.brm8 <- update(seagrass.brm8p, sample_prior = 'yes')

seagrass.brm8 |> SUYR_prior_and_posterior()


### Diagnostics ----
seagrass.brm8$fit |> stan_trace()
seagrass.brm8$fit |> stan_ac()
seagrass.brm8$fit |> stan_rhat()
seagrass.brm8$fit |> stan_ess()

seagrass.brm8 |> pp_check(type = 'dens_overlay', ndraws = 100)

seagrass.resids <- make_brms_dharma_res(seagrass.brm8, integerResponse = FALSE)
wrap_elements(~testUniformity(seagrass.resids)) +
  wrap_elements(~plotResiduals(seagrass.resids, form = factor(rep(1, nrow(seagrass))))) + 
  wrap_elements(~plotResiduals(seagrass.resids, quantreg = TRUE)) + 
  wrap_elements(~testDispersion(seagrass.resids)) 

save(seagrass.brm8, seagrass, priors8, brm8.form, file = '../data/modelled/Model8Beta.RData')

### Model investigation ----
seagrass.brm8 |> 
  conditional_effects(spaghetti = TRUE, ndraws = 100) |>
  plot(points = TRUE, ask = FALSE)

seagrass.brm8 |>
  as_draws_df() |>
  mutate(across(everything(), exp)) |> #to go from log scale to odds ratio
  summarise_draws(median, HDInterval::hdi, rhat, length, ess_bulk, ess_tail,
                  Pl = ~ mean(.x < 1),
                  Pg = ~ mean(.x > 1)) |>
  as_tibble() |>
  dplyr::slice(1:5)


## Compare models ----
loo::loo_compare(rstan::loo(seagrass.brm6),
                 rstan::loo(seagrass.brm7),
                 rstan::loo(seagrass.brm8))
