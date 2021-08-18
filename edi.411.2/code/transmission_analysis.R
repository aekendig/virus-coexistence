## Goal: combine transmisison data and concentration data and analyze treatment effects

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse) # version used: 1.2.1
library(cowplot) # version used: 0.9.4
library(brms)  # version used: 2.7.0
library(tidybayes)  # version used: 1.0.4
library(bayesplot) # versium used: 1.6.0

# import data
qdatp <- read_csv("./output/concentration_analysis_pav_data.csv")
qdatr <- read_csv("./output/concentration_analysis_rpv_data.csv")
tdat <- read_csv("./data/transmission_data.csv") # came from MeanqPCR_VCEBaseDataset_Trans_ManUp_123016.csv (unnecessary columns and rows removed)
sdat <- read_csv("./data/sample_exp_molc_data.csv")

# import models
# load("./output/lacroix_output/lacroix_transmission_pav.rda")
# load("./output/lacroix_output/lacroix_transmission_rpv.rda")

# functions
plot_intervals <- function(data) {
  ggplot(data, aes(y = parameter, yend = parameter)) + 
    geom_vline(xintercept = 1) +
    geom_segment(aes(x = ll, xend = hh), size = 1) + 
    geom_segment(aes(x = l, xend = h), size = 2) +
    geom_point(aes(x = m), size = 3, color = "red") +
    theme_bw() +
    facet_wrap(~model)
}
# modified from Tristan Mahr (https://www.tjmahr.com/ggplot2-how-to-do-nothing/)


#### edit concentration data ####

# add aphid mass to dataset and calculate new concentration metrics
sdat2 <- sdat %>%
  filter(material == "shoot") %>%
  select(time, inoc, nutrient, round, replicate, aphid_mass_g) %>%
  mutate(aphid_mass_mg = aphid_mass_g * 1000,
         aphid_mass_mg = replace_na(aphid_mass_mg, 0)) %>%
  select(-aphid_mass_g)

qdatp2 <- qdatp %>%
  left_join(sdat2) %>%
  mutate(quant_t = conc * aphid_mass_mg)

qdatr2 <- qdatr %>%
  left_join(sdat2) %>%
  mutate(quant_t = conc * aphid_mass_mg)

# replicates were combined in round 1 for transmission trial - average their concentration values and combine the amount of tissue provided

# subset data for round 1 plants and average over replicates
qdatpR1 <- qdatp2 %>%
  filter(round == 1) %>%
  group_by(dpi, time, inoc, nutrient, round) %>%
  mutate(total_aphid_mass_mg = sum(aphid_mass_mg),
         conc2 = conc * (aphid_mass_mg / total_aphid_mass_mg)) %>%
  summarise(aphid_mass_mg = unique(total_aphid_mass_mg),
            conc = sum(conc2),
            quant_t = sum(quant_t)) %>%
  ungroup() %>%
  mutate(replicate = 1)

qdatrR1 <- qdatr2 %>%
  filter(round == 1) %>%
  group_by(dpi, time, inoc, nutrient, round) %>%
  mutate(total_aphid_mass_mg = sum(aphid_mass_mg),
         conc2 = conc * (aphid_mass_mg / total_aphid_mass_mg)) %>%
  summarise(aphid_mass_mg = unique(total_aphid_mass_mg),
            conc = sum(conc2),
            quant_t = sum(quant_t)) %>%
  ungroup() %>%
  mutate(replicate = 1)

# recombine with rest of data
qdatp3 <- qdatp2 %>%
  filter(round != 1) %>%
  select(dpi, time, inoc, nutrient, round, replicate, aphid_mass_mg, conc, quant_t) %>%
  full_join(qdatpR1)

qdatr3 <- qdatr2 %>%
  filter(round != 1) %>%
  select(dpi, time, inoc, nutrient, round, replicate, aphid_mass_mg, conc, quant_t) %>%
  full_join(qdatrR1)


#### edit transmission data ####

# examine data
unique(tdat$material)
tdat %>% filter(is.na(material)) %>% data.frame() # remove these
unique(tdat$PAV_t)
unique(tdat$RPV_t)
unique(tdat$inoc)
unique(tdat$nutrient)

# remove root samples and unknown samples, create ID column
tdat2 <- tdat %>%
  filter(material == "shoot") %>%
  mutate(ID = paste(round, time, inoc, nutrient, replicate, nutrient_t, sep = "."))

# check duplicates
tdat2 %>%
  group_by(ID) %>%
  filter(n() > 1) %>%
  data.frame()
# 1.2.coinfection.N+P.2.N+P: different extraction dates
# 1.8.RPV.N.2.low: take out the one with PAV_t = 0.5 (checked gel)
# 2.4.PAV.P.3.N: take out the one with RPV_t = 0 (checked gel)
# 2.7.RPV.N.3.N: take out the one with PAV_t = 0 (checked gel)

tdat3 <- tdat2 %>%
  filter(!(ID == "1.2.coinfection.N+P.2.N+P" & ext_date_t == "3/23/16")) %>%
  filter(!(ID == "1.8.RPV.N.2.low" & PAV_t == 0.5)) %>%
  filter(!(ID == "2.4.PAV.P.3.N" & RPV_t == 0)) %>%
  filter(!(ID == "2.7.RPV.N.3.N" & PAV_t == 0))

# sample size (Table S2)
tdat3 %>%
  group_by(inoc, nutrient, nutrient_t) %>%
  summarise(PAV = sum(!is.na(PAV_t)), RPV = sum(!is.na(RPV_t))) %>%
  data.frame()


#### combine concentration and transmission data ####

# only keep overlapping data
# make N and P columns for transmission
# coinfection column for source plant
# transmission columns
# 0 concentration and quantity if missing
# scale concentration and quantity
datp <- inner_join(qdatp3, tdat3) %>%
  mutate(high_N = ifelse(nutrient %in% c("N", "N+P"), 1, 0),
         high_P = ifelse(nutrient %in% c("P", "N+P"), 1, 0),
         high_N_t = ifelse(nutrient_t %in% c("N", "N+P"), 1, 0),
         high_P_t = ifelse(nutrient_t %in% c("P", "N+P"), 1, 0),
         co = ifelse(inoc == "coinfection", 1, 0),
         t_up = ceiling(PAV_t),
         t_dn = floor(PAV_t))

datr <- inner_join(qdatr3, tdat3) %>%
  mutate(high_N = ifelse(nutrient %in% c("N", "N+P"), 1, 0),
         high_P = ifelse(nutrient %in% c("P", "N+P"), 1, 0),
         high_N_t = ifelse(nutrient_t %in% c("N", "N+P"), 1, 0),
         high_P_t = ifelse(nutrient_t %in% c("P", "N+P"), 1, 0),
         co = ifelse(inoc == "coinfection", 1, 0),
         t_up = ceiling(RPV_t),
         t_dn = floor(RPV_t))

# check for NA's
datp %>% 
  select(t_up, t_dn, conc, quant_t, high_N, high_P, high_N_t, high_P_t, co) %>%
  is.na() %>%
  colSums()

datr %>% 
  select(t_up, t_dn, conc, quant_t, high_N, high_P, high_N_t, high_P_t, co) %>%
  is.na() %>%
  colSums()


#### sample sizes and accidental infections ####

# accidental inoculations from all samples (Table S2)
datp %>%
  filter(inoc == "PAV") %>%
  group_by(nutrient, nutrient_t) %>%
  summarise(RPV_up = sum(ceiling(RPV_t) == 1), RPV_dn = sum(floor(RPV_t) == 1)) %>%
  data.frame()

datr %>%
  filter(inoc == "RPV") %>%
  group_by(nutrient, nutrient_t) %>%
  summarise(PAV_up = sum(ceiling(PAV_t) == 1), PAV_dn = sum(floor(PAV_t) == 1)) %>%
  data.frame()

# remove accidental infections and calculate standardized concentration and quantity
datp2 <- datp %>%
  mutate(RPV_up = ceiling(RPV_t)) %>%
  filter(co == 1 | RPV_up == 0) %>%
  select(-RPV_up) %>%
  mutate(conc_s = (conc - mean(conc)) / sd(conc),
         quant_s = (quant_t - mean(quant_t)) / sd(quant_t))

datr2 <- datr %>%
  mutate(PAV_up = ceiling(PAV_t)) %>%
  filter(co == 1 | PAV_up == 0) %>%
  select(-PAV_up) %>%
  mutate(conc_s = (conc - mean(conc)) / sd(conc),
         quant_s = (quant_t - mean(quant_t)) / sd(quant_t))

# new sample sizes (Table S2)
datp2 %>%
  group_by(co, nutrient, nutrient_t) %>%
  summarize(n = n()) %>%
  data.frame()

datr2 %>%
  group_by(co, nutrient, nutrient_t) %>%
  summarize(n = n()) %>%
  data.frame()

# new source plant concentrations
datp2 %>%
  select(co, nutrient, dpi, replicate, conc) %>%
  unique() %>%
  group_by(co, nutrient) %>%
  summarize(conc_mean = mean(conc), n = n()) %>%
  data.frame()

qdatp %>%
  group_by(co, nutrient) %>%
  summarize(conc_mean = mean(conc), n = n()) %>%
  data.frame()

datr2 %>%
  select(co, nutrient, dpi, replicate, conc) %>%
  unique() %>%
  group_by(co, nutrient) %>%
  summarize(conc_mean = mean(conc), n = n()) %>%
  data.frame()

qdatr %>%
  group_by(co, nutrient) %>%
  summarize(conc_mean = mean(conc), n = n()) %>%
  data.frame()

# traits of unintended infection plants
datp %>%
  anti_join(datp2 %>%
              select(-c(conc_s, quant_s))) %>%
  group_by(nutrient) %>%
  summarise(n = n()) # fewer N's

datp %>%
  anti_join(datp2 %>%
              select(-c(conc_s, quant_s))) %>%
  group_by(nutrient_t) %>%
  summarise(n = n())

datp %>%
  anti_join(datp2 %>%
              select(-c(conc_s, quant_s))) %>%
  group_by(dpi) %>%
  summarise(n = n()) # middle dpi

datr %>%
  anti_join(datr2 %>%
              select(-c(conc_s, quant_s))) %>%
  group_by(nutrient) %>%
  summarise(n = n()) # fewer low and N+P

datr %>%
  anti_join(datr2 %>%
              select(-c(conc_s, quant_s))) %>%
  group_by(nutrient_t) %>%
  summarise(n = n())

datr %>%
  anti_join(datr2 %>%
              select(-c(conc_s, quant_s))) %>%
  group_by(dpi) %>%
  summarise(n = n()) # middle dpi


#### visualize sources of variation ####

datp2 %>%
  ggplot(aes(x = round, y = t_up)) +
  stat_summary(geom = "point", size = 2, fun.y = "mean") +
  stat_summary(geom = "errorbar", width = 0.1, fun.data = "mean_cl_boot") +
  facet_wrap(~inoc) # 3 lower for single

datr2 %>%
  ggplot(aes(x = round, y = t_up)) +
  stat_summary(geom = "point", size = 2, fun.y = "mean") +
  stat_summary(geom = "errorbar", width = 0.1, fun.data = "mean_cl_boot") +
  facet_wrap(~inoc) # 3 higher for coinfection

datp2 %>%
  ggplot(aes(x = as.factor(time), y = t_up)) +
  stat_summary(geom = "point", size = 2, fun.y = "mean") +
  stat_summary(geom = "errorbar", width = 0.1, fun.data = "mean_cl_boot") +
  facet_wrap(~inoc) # increases with time

datp2 %>%
  ggplot(aes(x = as.factor(time), y = conc_s)) +
  stat_summary(geom = "point", size = 2, fun.y = "mean") +
  stat_summary(geom = "errorbar", width = 0.1, fun.data = "mean_cl_boot") +
  facet_wrap(~inoc) # so does conc

cor.test(datp2$conc_s, datp2$time) # 0.25

datp2 %>%
  ggplot(aes(x = as.factor(time), y = quant_s)) +
  stat_summary(geom = "point", size = 2, fun.y = "mean") +
  stat_summary(geom = "errorbar", width = 0.1, fun.data = "mean_cl_boot") +
  facet_wrap(~inoc) # and quantity

cor.test(datp2$quant_s, datp2$time) #0.27

datr2 %>%
  ggplot(aes(x = as.factor(time), y = t_up)) +
  stat_summary(geom = "point", size = 2, fun.y = "mean") +
  stat_summary(geom = "errorbar", width = 0.1, fun.data = "mean_cl_boot") +
  facet_wrap(~inoc) # constant with time

datr2 %>%
  ggplot(aes(x = as.factor(time), y = conc_s)) +
  stat_summary(geom = "point", size = 2, fun.y = "mean") +
  stat_summary(geom = "errorbar", width = 0.1, fun.data = "mean_cl_boot") +
  facet_wrap(~inoc) # conc increases

cor.test(datr2$conc_s, datr2$time) # 0.19

datr2 %>%
  ggplot(aes(x = as.factor(time), y = quant_s)) +
  stat_summary(geom = "point", size = 2, fun.y = "mean") +
  stat_summary(geom = "errorbar", width = 0.1, fun.data = "mean_cl_boot") +
  facet_wrap(~inoc) # so does quantity

cor.test(datr2$conc_s, datr2$time) # 0.19

datp2 %>%
  ggplot(aes(x = nutrient, y = t_up)) +
  stat_summary(geom = "point", size = 2, fun.y = "mean") +
  stat_summary(geom = "errorbar", width = 0.1, fun.data = "mean_cl_boot") +
  facet_grid(nutrient_t ~ inoc) # look similar

datr2 %>%
  ggplot(aes(x = nutrient, y = t_up)) +
  stat_summary(geom = "point", size = 2, fun.y = "mean") +
  stat_summary(geom = "errorbar", width = 0.1, fun.data = "mean_cl_boot") +
  facet_grid(nutrient_t ~ inoc) # look similar

datp2 %>%
  filter(inoc == "PAV") %>%
  ggplot(aes(x = conc_s, y = t_up)) +
  geom_point() +
  facet_grid(nutrient_t ~ nutrient) 
# sample sizes small when doing all of the combinations

datp2 %>%
  ggplot(aes(x = conc_s, y = t_up)) +
  geom_point() +
  facet_grid(inoc ~ nutrient) 

datr2 %>%
  ggplot(aes(x = conc_s, y = t_up)) +
  geom_point() +
  facet_grid(inoc ~ nutrient) 


#### models with uniformative priors and simplest format ####

# can't do autoregressive models with bernoulli
# 
# mpuc <- brm(data = datp2, family = bernoulli,
#            t_up ~ conc_s + co * (high_N * high_P + high_N_t * high_P_t) + (1|round) + (1|time),
#            prior = c(prior(normal(0, 10), class = Intercept),
#                      prior(normal(0, 10), class = b),
#                      prior(cauchy(0, 1), class = sd)),
#            iter = 6000, warmup = 1000, chains = 3, cores = 2,
#            control = list(adapt_delta = 0.9999))
# summary(mpuc)
# save(mpuc, file = "./output/transmission_pav_up_concentration.rda")

# mpuq <- update(mpuc, formula. = t_up ~ quant_s + co * (high_N * high_P + high_N_t * high_P_t) + (1|round) + (1|time), newdata = datp2)
# summary(mpuq)
# save(mpuq, file = "./output/transmission_pav_up_quantity.rda")

# mruc <- update(mpuc, formula. = t_up ~ conc_s + co * (high_N * high_P + high_N_t * high_P_t) + (1|round) +(1|time), newdata = datr2)
# summary(mruc)
# save(mruc, file = "./output/transmission_rpv_up_concentration.rda")
  
# mruq <- update(mruc, formula. = t_up ~ quant_s + co * (high_N * high_P + high_N_t * high_P_t) + (1|round) + (1|time), newdata = datr2, control = list(adapt_delta = 0.99999))
# summary(mruq)
# save(mruq, file = "./output/transmission_rpv_up_quantity.rda")
  
# mpdc <- update(mpuc, formula. = t_dn ~ conc_s + co * (high_N * high_P + high_N_t * high_P_t) + (1|round) + (1|time), newdata = datp2, control = list(adapt_delta = 0.99999999))
# summary(mpdc)
# save(mpdc, file = "./output/transmission_pav_down_concentration.rda")
  
# mpdq <- update(mpdc, formula. = t_dn ~ quant_s + co * (high_N * high_P + high_N_t * high_P_t) + (1|round) + (1|time), newdata = datp2, control = list(adapt_delta = 0.9999999))
# summary(mpdq)
# save(mpdq, file = "./output/transmission_pav_down_quantity.rda")
  
# mrdc <- update(mruc, formula. = t_dn ~ conc_s + co * (high_N * high_P + high_N_t * high_P_t) + (1|round) + (1|time), newdata = datr2, control = list(adapt_delta = 0.9999999))
# summary(mrdc)
# save(mrdc, file = "./output/transmission_rpv_down_concentration.rda")
  
# mrdq <- update(mruc, formula. = t_dn ~ quant_s + co * (high_N * high_P + high_N_t * high_P_t) + (1|round) + (1|time), newdata = datr2)
# summary(mrdq)
# save(mrdq, file = "./output/transmission_rpv_down_quantity.rda")

  
#### check models ####

# convergence of chains
# plot(mpuc)
# plot(mpuq)
# plot(mpdc)
# plot(mpdq)

# plot(mruc)
# plot(mruq)
# plot(mrdc)
# plot(mrdq)

# posterior predictive check
# pp_check(mpuc, nsamples = 50)
# pp_check(mpuq, nsamples = 50)
# pp_check(mpdc, nsamples = 50)
# pp_check(mpdq, nsamples = 50)

# pp_check(mruc, nsamples = 50)
# pp_check(mruq, nsamples = 50)
# pp_check(mrdc, nsamples = 50)
# pp_check(mrdq, nsamples = 50)

# compare concentration and quantity
# loo_mpu <- list(loo(mpuc), loo(mpuq))
# loo_mpu
# loo::loo_compare(loo_mpu) # pretty much equal
# 
# loo_mpd <- list(loo(mpdc), loo(mpdq))
# loo_mpd
# loo::loo_compare(loo_mpd) # equal
# 
# loo_mru <- list(loo(mruc), loo(mruq))
# loo_mru
# loo::loo_compare(loo_mru) # equal
# 
# loo_mrd <- list(loo(mrdc), loo(mrdq))
# loo_mrd
# loo::loo_compare(loo_mrd) # equal

# compare estimates for up and down
# params = c("b_Intercept", "b_conc_s", "b_co", "b_high_N", "b_high_P", "b_high_N_t", "b_high_P_t", "b_high_N:high_P", "b_high_N_t:high_P_t", "b_co:high_N", "b_co:high_P", "b_co:high_N_t", "b_co:high_P_t", "b_co:high_N:high_P", "b_co:high_N_t:high_P_t")
# 
# apuc <- as.array(mpuc) %>%
#   mcmc_intervals_data(pars = params) %>%
#   mutate(model = "PAV up conc uninformative")
# 
# apdc <- as.array(mpdc) %>%
#   mcmc_intervals_data(pars = params) %>%
#   mutate(model = "PAV down conc uninformative")
# 
# plot_intervals(full_join(apuc, apdc)) 
# # rounding up leads to larger effect sizes
# 
# aruc <- as.array(mruc) %>%
#   mcmc_intervals_data(pars = params) %>%
#   mutate(model = "RPV up conc uninformative")
# 
# ardc <- as.array(mrdc) %>%
#   mcmc_intervals_data(pars = params) %>%
#   mutate(model = "RPV down conc uninformative")
# 
# plot_intervals(full_join(aruc, ardc))
# very similar


#### models with informative priors ####

# compare estimates for prior model and model without priors
# summary(mp)
# aplc <- (as.array(mp)) %>%
#   mcmc_intervals_data(pars = c("b_Intercept", "b_conc_s", "b_co", "b_high_N", "b_high_P", "b_co:high_N", "b_high_N:high_P")) %>%
#   mutate(model = "PAV Lacroix")
# plot_intervals(full_join(apuc, aplc))
# some estimates in the opposite direction, but error bars generally larger

# summary(mr)
# arlc <- (as.array(mr)) %>%
#   mcmc_intervals_data(pars = c("b_Intercept", "b_conc_s", "b_co", "b_high_N", "b_high_P", "b_co:high_N", "b_co:high_P", "b_high_N:high_P", "b_co:high_N:high_P")) %>%
#   mutate(model = "RPV Lacroix")
# plot_intervals(full_join(aruc, arlc))
# estimates closer to zero

# PAV model
# summary(mp)
# 
# mpuci <- update(mpuc,
#                 prior = c(prior(normal(1.53, 0.60), class = Intercept),
#                           prior(normal(-0.26, 0.28), class = b, coef = conc_s),
#                           prior(normal(6.05, 4.17), class = b, coef = co),
#                           prior(normal(0.69, 0.94), class = b, coef = high_N),
#                           prior(normal(0.63, 1.42), class = b, coef = high_P),
#                           prior(normal(-6.51, 4.23), class = b, coef = co:high_N),
#                           prior(normal(0, 10), class = b, coef = co:high_P),
#                           prior(normal(-0.20, 1.63), class = b, coef = high_N:high_P),
#                           prior(normal(0, 10), class = b, coef = co:high_N:high_P),
#                           prior(normal(0, 10), class = b, coef = high_N_t),
#                           prior(normal(0, 10), class = b, coef = high_P_t),
#                           prior(normal(0, 10), class = b, coef = co:high_N_t),
#                           prior(normal(0, 10), class = b, coef = co:high_P_t),
#                           prior(normal(0, 10), class = b, coef = high_N_t:high_P_t),
#                           prior(normal(0, 10), class = b, coef = co:high_N_t:high_P_t),
#                           prior(cauchy(0, 1), class = sd)))
# 
# summary(mpuci)
# plot(mpuci)
# 
# save(mpuci, file = "./output/transmission_pav_up_concentration_informative.rda")
# 
# apuci <- as.array(mpuci) %>%
#   mcmc_intervals_data(pars = params) %>%
#   mutate(model = "PAV up conc informative")
# 
# plot_intervals(full_join(apuc, apuci)) # slightly smaller errorbars, similar mean values
# 
# # RPV model
# summary(mr)
# 
# mruci <- update(mruc,
#                 prior = c(prior(normal(0.54, 0.30), class = Intercept),
#                           prior(normal(0.19, 0.13), class = b, coef = conc_s),
#                           prior(normal(-0.19, 0.44), class = b, coef = co),
#                           prior(normal(0.62, 0.42), class = b, coef = high_N),
#                           prior(normal(0.13, 0.71), class = b, coef = high_P),
#                           prior(normal(0.33, 0.60), class = b, coef = co:high_N),
#                           prior(normal(-0.10, 1.12), class = b, coef = co:high_P),
#                           prior(normal(-0.01, 0.86), class = b, coef = high_N:high_P),
#                           prior(normal(-0.60, 1.28), class = b, coef = co:high_N:high_P),
#                           prior(normal(0, 10), class = b, coef = high_N_t),
#                           prior(normal(0, 10), class = b, coef = high_P_t),
#                           prior(normal(0, 10), class = b, coef = co:high_N_t),
#                           prior(normal(0, 10), class = b, coef = co:high_P_t),
#                           prior(normal(0, 10), class = b, coef = high_N_t:high_P_t),
#                           prior(normal(0, 10), class = b, coef = co:high_N_t:high_P_t),
#                           prior(cauchy(0, 1), class = sd)))
# 
# summary(mruci)
# plot(mruci)
# 
# save(mruci, file = "./output/transmission_rpv_up_concentration_informative.rda")
# 
# aruci <- as.array(mruci) %>%
#   mcmc_intervals_data(pars = params) %>%
#   mutate(model = "RPV up conc informative")
# 
# plot_intervals(full_join(aruc, aruci)) # smaller error bars, similar mean values


#### models with concentration interactions and informative priors ####

# PAV model
# summary(mp)
mpcii <- brm(data = datp2, family = bernoulli,
             formula = t_up ~ conc_s * (high_N * high_P + high_N_t * high_P_t) + co * (high_N * high_P + high_N_t * high_P_t) + (1|round) + (1|time),
             prior = c(prior(normal(0, 10), class = Intercept),
                       prior(normal(-0.22, 0.25), class = b, coef = conc_s),
                       prior(normal(6.05, 4.13), class = b, coef = co),
                       prior(normal(0.64, 0.78), class = b, coef = high_N),
                       prior(normal(-1.18, 0.71), class = b, coef = high_P),
                       prior(normal(0, 10), class = b, coef = high_N_t),
                       prior(normal(0, 10), class = b, coef = high_P_t),
                       prior(normal(2.07, 0.93), class = b, coef = high_N:high_P),
                       prior(normal(0, 10), class = b, coef = high_N_t:high_P_t),
                       prior(normal(0, 10), class = b, coef = conc_s:high_N),
                       prior(normal(0, 10), class = b, coef = conc_s:high_P),
                       prior(normal(0, 10), class = b, coef = conc_s:high_N_t),
                       prior(normal(0, 10), class = b, coef = conc_s:high_P_t),
                       prior(normal(0, 10), class = b, coef = conc_s:high_N:high_P),
                       prior(normal(0, 10), class = b, coef = conc_s:high_N_t:high_P_t),
                       prior(normal(-7.67, 4.16), class = b, coef = high_N:co),
                       prior(normal(0, 10), class = b, coef = high_P:co),
                       prior(normal(0, 10), class = b, coef = high_N_t:co),
                       prior(normal(0, 10), class = b, coef = high_P_t:co),
                       prior(normal(0, 10), class = b, coef = high_N:high_P:co),
                       prior(normal(0, 10), class = b, coef = high_N_t:high_P_t:co),
                       prior(cauchy(0, 1), class = sd)),
             iter = 6000, warmup = 1000, chains = 3, cores = 2,
             control = list(adapt_delta = 0.9999))
save(mpcii, file = "./output/transmission_pav_up_concentration_interaction_informative.rda")

summary(mpcii)
plot(mpcii)
pp_check(mpcii, nsamples = 100)

# RPV model
# summary(mr)
mrcii <- brm(data = datr2, family = bernoulli,
             formula = t_up ~ conc_s * (high_N * high_P + high_N_t * high_P_t) + co * (high_N * high_P + high_N_t * high_P_t) + (1|round) + (1|time),
             prior = c(prior(normal(0, 10), class = Intercept),
                       prior(normal(0.20, 0.13), class = b, coef = conc_s),
                       prior(normal(-0.17, 0.45), class = b, coef = co),
                       prior(normal(0.64, 0.42), class = b, coef = high_N),
                       prior(normal(0.12, 0.70), class = b, coef = high_P),
                       prior(normal(0, 10), class = b, coef = high_N_t),
                       prior(normal(0, 10), class = b, coef = high_P_t),
                       prior(normal(0, 0.86), class = b, coef = high_N:high_P),
                       prior(normal(0, 10), class = b, coef = high_N_t:high_P_t),
                       prior(normal(0, 10), class = b, coef = conc_s:high_N),
                       prior(normal(0, 10), class = b, coef = conc_s:high_P),
                       prior(normal(0, 10), class = b, coef = conc_s:high_N:high_P),
                       prior(normal(0, 10), class = b, coef = conc_s:high_N_t),
                       prior(normal(0, 10), class = b, coef = conc_s:high_P_t),
                       prior(normal(0, 10), class = b, coef = conc_s:high_N_t:high_P_t),
                       prior(normal(0.31, 0.60), class = b, coef = high_N:co),
                       prior(normal(-0.69, 0.97), class = b, coef = high_P:co),
                       prior(normal(-0.01, 1.14), class = b, coef = high_N:high_P:co),
                       prior(normal(0, 10), class = b, coef = high_N_t:co),
                       prior(normal(0, 10), class = b, coef = high_P_t:co),
                       prior(normal(0, 10), class = b, coef = high_N_t:high_P_t:co),
                       prior(cauchy(0, 1), class = sd)),
             iter = 6000, warmup = 1000, chains = 3, cores = 2,
             control = list(adapt_delta = 0.9999, max_treedepth = 15))
save(mrcii, file = "./output/transmission_rpv_up_concentration_interaction_informative.rda")

summary(mrcii)
plot(mrcii)
pp_check(mrcii, nsamples = 100)


#### models with concentration interactions and uinformative priors ####

# PAV model
mpciu <- brm(data = datp2, family = bernoulli,
             formula = t_up ~ conc_s * (high_N * high_P + high_N_t * high_P_t) + co * (high_N * high_P + high_N_t * high_P_t) + (1|round) + (1|time),
             prior = c(prior(normal(0, 10), class = Intercept),
                       prior(normal(0, 10), class = b),
                       prior(cauchy(0, 1), class = sd)),
             iter = 6000, warmup = 1000, chains = 3, cores = 2,
             control = list(adapt_delta = 0.9999))
save(mpciu, file = "./output/transmission_pav_up_concentration_interaction_uninformative.rda")

summary(mpciu)
plot(mpciu)
pp_check(mpciu, nsamples = 100)

# RPV model
mrciu <- update(mpciu, newdata = datr2)
save(mrciu, file = "./output/transmission_rpv_up_concentration_interaction_uninformative.rda")

summary(mrciu)
plot(mrciu)
pp_check(mrciu, nsamples = 100)

#### save data for plotting ####

# save file
write_csv(datp2, "./output/transmission_analysis_pav_data.csv")
write_csv(datr2, "./output/transmission_analysis_rpv_data.csv")

