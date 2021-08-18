## Goal: analyze treatment effects on binary infection data 

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse) # version used: 1.2.1
library(brms) # version used: 2.7.0
library(tidybayes) # version used: 1.0.4
library(bayesplot) # version used: 1.6.0

# import data
dat <- read_csv("./output/concentration_analysis_all_data.csv")
# samples with uncertainty about presence/absence have been removed (remove == 1 from qPCR_raw_data_processing)
# samples above the standard curve are excluded
# the only difference between this and the concentration dataset is that it includes values that are effectively 0


#### edit data ####

dat2 <- dat %>%
  filter(inoc %in% c("PAV", "coinfection", "RPV")) %>%
  mutate(co = ifelse(inoc == "coinfection", 1, 0),
         present = ifelse(quant_zero == 1, 0, 1))

# check for same sample in multiple qPCR groups
dups <-dat2 %>%
  group_by(target, sample) %>%
  mutate(dup = duplicated(sample)) %>%
  filter(dup == T) %>%
  select(sample, target, dup) %>%
  ungroup()

dups %>%
  left_join(dat2) %>%
  select(sample, target, quant_adj) %>%
  data.frame() 
# 22 duplicates, 44 samples

# remove duplicates and qPCR-specific info
dat3 <- dat2 %>%
  select(target:mass_ext_mg, quant_adj, co, present) %>%
  distinct(sample, target, keep_all = F)

# data by virus
datp <- dat2 %>%
  filter(target == "PAV" & inoc != "RPV") 

datr <- dat2 %>%
  filter(target == "RPV" & inoc != "PAV")

# save data
write_csv(datp, "./output/infection_analysis_pav_data.csv")
write_csv(datr, "./output/infection_analysis_rpv_data.csv")


#### visualize data ####

# PAV binomial
datp %>%
  ggplot(aes(x = dpi, y = present)) +
  stat_summary(geom = "point", fun.y = "mean") +
  stat_summary(geom = "line", fun.y = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1) +
  facet_grid(inoc~nutrient)

datp %>%
  group_by(dpi, inoc, nutrient) %>%
  summarise(tot = sum(present)) %>%
  ggplot(aes(x = as.factor(dpi), y = tot)) +
  geom_bar(stat = "identity") +
  facet_grid(inoc~nutrient)
  
# RPV binomial
datr %>%
  ggplot(aes(x = dpi, y = present)) +
  stat_summary(geom = "point", fun.y = "mean") +
  stat_summary(geom = "line", fun.y = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1) +
  facet_grid(inoc~nutrient)

datr %>%
  group_by(dpi, inoc, nutrient) %>%
  summarise(tot = sum(present)) %>%
  ggplot(aes(x = as.factor(dpi), y = tot)) +
  geom_bar(stat = "identity") +
  facet_grid(inoc~nutrient)


#### statistical models ####

## PAV model 
# uninformative priors
# you can't currently implement autoregressive models, time included as a random effect
m.bu.p <- brm(data = datp, family = bernoulli,
              present ~ co * high_N * high_P + (1|time),
              prior <- c(prior(normal(0, 10), class = Intercept),
                         prior(normal(0, 10), class = b),
                         prior(cauchy(0, 1), class = sd)),
              iter = 6000, warmup = 1000, chains = 3, cores = 2)

# inspect model
summary(m.bu.p) 
plot(m.bu.p) 
pp_check(m.bu.p, nsamples = 100)

# save model
save(m.bu.p, file = "./output/infection_analysis_uninformative_pav.rda")


## RPV model
# uninformative priors
m.bu.r <- update(m.bu.p, newdata = datr)

# inspect model
summary(m.bu.r) 
plot(m.bu.r) 
pp_check(m.bu.r, nsamples = 100)

# save model
save(m.bu.r, file = "./output/infection_analysis_uninformative_rpv.rda")


## PAV model 
# informative priors
# Lacroix et al. 2014
m.bi.p <- update(m.bu.p,
                 prior = c(prior(normal(0.82, 0.525), class = Intercept),
                           prior(normal(0.001, 0.002), class = b, coef = high_N),
                           prior(normal(0.018, 0.017), class = b, coef = high_P),
                           prior(normal(-0.874, 0.69), class = b, coef = co),
                           prior(normal(0, 1), class = b, coef = high_N:high_P),
                           prior(normal(0.004, 0.003), class = b, coef = co:high_N),
                           prior(normal(-0.001, 0.022), class = b, coef = co:high_P),
                           prior(normal(0, 1), class = b, coef = co:high_N:high_P),
                           prior(cauchy(0, 1), class = sd)),
              iter = 6000, warmup = 1000, chains = 3, cores = 2)

# inspect model
summary(m.bi.p) # pretty strong changes in parameter estimates because of the priors
plot(m.bi.p)
pp_check(m.bi.p, nsamples = 100)

# save model
save(m.bi.p, file = "./output/infection_analysis_informative_pav.rda")


## RPV model 
# informative priors
# Lacroix et al. 2014
m.bi.r <- update(m.bu.r,
                 prior = c(prior(normal(1.162, 0.555), class = Intercept),
                           prior(normal(-0.002, 0.002), class = b, coef = high_N),
                           prior(normal(-0.051, 0.016), class = b, coef = high_P),
                           prior(normal(-1.785, 0.729), class = b, coef = co),
                           prior(normal(0, 1), class = b, coef = high_N:high_P),
                           prior(normal(0.005, 0.003), class = b, coef = co:high_N),
                           prior(normal(0.028, 0.023), class = b, coef = co:high_P),
                           prior(normal(0, 1), class = b, coef = co:high_N:high_P),
                           prior(cauchy(0, 1), class = sd)),
                 iter = 6000, warmup = 1000, chains = 3, cores = 2)

# inspect model
summary(m.bi.r) # pretty strong changes in values because of priors
plot(m.bi.r) 
pp_check(m.bi.r, nsamples = 100)

# save model
save(m.bi.r, file = "./output/infection_analysis_informative_rpv.rda")
