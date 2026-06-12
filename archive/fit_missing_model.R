library("tidyverse")
library("brms")

load("dat_test.rds")
dat_test <- dat_test %>%
      mutate(
        Class = as.factor(Class),
        log_sampling = log(Sampling_effort)
        )

dat_no_grassland <- dat_test %>%
        filter(Biome %in% c("Mangroves", "Mediterranean Forests, Woodlands & Scrub","Temperate Broadleaf & Mixed Forests","Tropical & Subtropical Dry Broadleaf Forests","Tropical & Subtropical Moist Broadleaf Forests"))

brm_no_grass <-
  brm(data = dat_no_grassland, family = negbinomial,
      Measurement ~ 1 + offset(log_sampling) + scaled_yr,
      prior = c(prior(normal(0, 10), class = Intercept),
                prior(normal(0, 1), class = b),
                prior(gamma(0.01, 0.01), class = shape) # the brms default
                ),
      cores = 4,
      chains = 4,
      iter = 1000)

save(brm_no_grass, file = "brm_no_grass.rda")
