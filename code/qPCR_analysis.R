## Goal: analyze qPCR data


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(nlme)

# import data
dat <- read_csv("intermediate-data/qPCR_expt_data_cleaned.csv")


#### edit data ####

# make wide by target
# for each invasion trial, determine whether resident established
# for each single inoculation, determine whether there was contamination
# identify samples missing one of the viruses
# make column for roles
dat2 <- dat %>%
  select(-c(quant_adj, quant.ul)) %>% # remove other quantity calculations
  pivot_wider(names_from = target,
              values_from = quant.mg,
              names_glue = "{target}_quant.mg") %>% # make wide by target
  mutate(resident_est = case_when(invasion == "I" & first_inoculation == "PAV" & PAV_quant.mg > 0  ~ 1,
                                  invasion == "I" & first_inoculation == "RPV" & RPV_quant.mg > 0  ~ 1,
                                  TRUE ~ 0),
         single_cont = case_when(invasion == "S" & first_inoculation == "PAV" & RPV_quant.mg > 0  ~ 1,
                                 invasion == "S" & first_inoculation == "RPV" & PAV_quant.mg > 0  ~ 1,
                                 TRUE ~ 0),
         missing_quant = case_when(is.na(RPV_quant.mg) | is.na(PAV_quant.mg) ~ 1,
                                   TRUE ~ 0),
         PAV_role = case_when(invasion == "S" & first_inoculation == "PAV" ~ "single",
                              invasion == "I" & first_inoculation == "PAV" ~ "resident",
                              invasion == "I" & first_inoculation == "RPV" ~ "invader",
                              TRUE ~ NA_character_),
         RPV_role = case_when(invasion == "S" & first_inoculation == "RPV" ~ "single",
                              invasion == "I" & first_inoculation == "RPV" ~ "resident",
                              invasion == "I" & first_inoculation == "PAV" ~ "invader",
                              TRUE ~ NA_character_),
         dpiSI = case_when(invasion == "S" ~ dpiR, # DPI that can be used for single and invaders
                           invasion == "I" ~ dpiI)) # incorrect for residents

# samples with each issue
dat2 %>%
  filter(missing_quant == 1)
# 56 missing one of the viruses

dat2 %>%
  filter(invasion == "I" & resident_est == 0) 
# 32 failed invasions

dat2 %>%
  filter(invasion == "S" & single_cont == 1)
# 92 contaminated single infections out of 197

dat2 %>%
  filter(invasion == "S") %>%
  ggplot(aes(x = log10(RPV_quant.mg), color = as.factor(single_cont))) +
  geom_density()
# clear separation between contamination and inoculation

dat2 %>%
  filter(invasion == "S") %>%
  ggplot(aes(x = log10(PAV_quant.mg), color = as.factor(single_cont))) +
  geom_density()
# no clear separation

# remove failed invasions
# contaminated single inoculations
# and missing quantities
dat3 <- dat2 %>%
  filter(!(invasion == "I" & resident_est == 0) & 
           !(invasion == "S" & single_cont == 1) &
           missing_quant != 1)

# sample sizes
dat3 %>%
  group_by(nutrient, first_inoculation, second_inoculation, time, invasion) %>%
  count() %>%
  ggplot(aes(x = n)) +
  geom_histogram()
# may need to remove some treatments for analyses

# PAV quantities visualization
dat3 %>%
  filter(!is.na(PAV_role)) %>%
  ggplot(aes(x = time, y = PAV_quant.mg, color = PAV_role)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  facet_wrap(~ nutrient)
# no single PAV data
# PAV growth may be higher with N

# RPV quantities visualization
dat3 %>%
  filter(!is.na(RPV_role)) %>%
  ggplot(aes(x = time, y = RPV_quant.mg, color = RPV_role)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  facet_wrap(~ nutrient)
# single and resident are very close
# RPV growth may be higher with N

# PAV single - why are they missing?
(PAV_single_summary <- dat2 %>%
  filter(PAV_role == "single") %>%
  summarise(missing = sum(missing_quant),
            contaminated = sum(single_cont),
            total = n(),
            maxRPV = max(RPV_quant.mg, na.rm = T)) %>%
  mutate(log_maxRPV = log10(maxRPV)))
# 92 out of 100 had some RPV
# max RPV is 171937, consistent with density figure

# same for RPV?
dat2 %>%
  filter(RPV_role == "single") %>%
  summarise(missing = sum(missing_quant),
            contaminated = sum(single_cont),
            total = n())
# no contamination

# replace contaminated single inoculations (use dat2)
# remove missing quantities
# make RPV quant NA below threshold
dat4 <- dat2 %>%
  filter(missing_quant != 1) %>%
  mutate(RPV_quant.mg = case_when(RPV_quant.mg < PAV_single_summary$maxRPV ~ NA_real_,
                                  TRUE ~ RPV_quant.mg),
         resident_est = case_when(invasion == "I" & first_inoculation == "PAV" & PAV_quant.mg > 0  ~ 1,
                                  invasion == "I" & first_inoculation == "RPV" & RPV_quant.mg > 0  ~ 1,
                                  TRUE ~ 0))

# failed invasions removed with missing quantity
dat2 %>%
  filter(invasion == "I" & resident_est == 0 & missing_quant == 1)
# 27

# new failed invasion count
dat4 %>%
  filter(invasion == "I" & resident_est == 0) 
# 6 (only added one additional)

# remove failed invasions
dat5 <- dat4 %>%
  filter(!(invasion == "I" & resident_est == 0))
# recovered almost 100 samples

# sample sizes
dat5 %>%
  group_by(nutrient, first_inoculation, second_inoculation, time, invasion) %>%
  count() %>%
  ggplot(aes(x = n)) +
  geom_histogram()
# still have some small sample sizes

# PAV quantities visualization
dat5 %>%
  filter(!is.na(PAV_role)) %>%
  ggplot(aes(x = time, y = PAV_quant.mg, color = PAV_role)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  facet_wrap(~ nutrient)

# RPV quantities visualization
dat5 %>%
  filter(!is.na(RPV_role)) %>%
  ggplot(aes(x = time, y = RPV_quant.mg, color = RPV_role)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  facet_wrap(~ nutrient)
# invasions still visible

# divide datasets by virus
PAVdat <- dat5 %>%
  filter(!is.na(PAV_role) & !is.na(PAV_quant.mg)) %>%
  mutate(quant.mg = PAV_quant.mg)

RPVdat <- dat5 %>%
  filter(!is.na(RPV_role) & !is.na(RPV_quant.mg)) %>%
  mutate(quant.mg = RPV_quant.mg)


#### estimate growth parameters ####

# visualize growth of invaders and single as one process
PAVdat %>%
  filter(PAV_role %in% c("single", "invader")) %>%
  ggplot(aes(dpiSI, PAV_quant.mg, color = PAV_role)) +
  stat_summary(geom = "line", fun = "mean") +
  geom_point() +
  facet_wrap(~ nutrient)
# invaders don't reach K

RPVdat %>%
  filter(RPV_role %in% c("single", "invader")) %>%
  ggplot(aes(dpiSI, RPV_quant.mg, color = RPV_role)) +
  stat_summary(geom = "line", fun = "mean") +
  geom_point() +
  facet_wrap(~ nutrient)
# invaders don't reach K

# select resident
PAVRdat <- PAVdat %>%
  filter(PAV_role == "resident") %>%
  select(set, nutrient, highN, highP, time, dpiI, replicate, quant.mg)

RPVRdat <- RPVdat %>%
  filter(RPV_role == "resident") %>%
  select(set, nutrient, highN, highP, time, dpiI, replicate, quant.mg)

# estimate carrying capacity
PAVmodK <- lm(quant.mg ~ highN * highP, data = PAVRdat)
summary(PAVmodK) # no sig difference among treatments
PAVK <- PAVRdat %>%
  select(nutrient, highN, highP) %>%
  unique() %>%
  mutate(K = predict(PAVmodK, newdata = .))

RPVmodK <- lm(quant.mg ~ highN * highP, data = RPVRdat)
summary(RPVmodK)  # no sig difference among treatments
RPVK <- RPVRdat %>%
  select(nutrient, highN, highP) %>%
  unique() %>%
  mutate(K = predict(RPVmodK, newdata = .))

# select invader
PAVIdat <- PAVdat %>%
  filter(PAV_role == "invader") %>%
  select(set, nutrient, highN, highP, time, dpiI, replicate, quant.mg) %>%
  left_join(PAVK)

RPVIdat <- RPVdat %>%
  filter(RPV_role == "invader") %>%
  select(set, nutrient, highN, highP, time, dpiI, replicate, quant.mg) %>%
  left_join(RPVK)

### parameters vary by treatment (4) ###

# simulated data/starting values
PAVIdat %>%
  group_by(nutrient) %>%
  summarise(minVal = min(quant.mg))

PAVIsim <- PAVIdat %>%
  select(dpiI) %>%
  unique() %>%
  expand_grid(PAVK) %>%
  mutate(N0 = 10,
         r = 0.1,
         quant.mg = K*N0/(N0 + (K-N0) * exp(-r * dpiI)))
# I think N0 = 1 and r = 0.4 are better starting values
# but I get a singular gradient with these

ggplot(PAVIdat, aes(dpiI, quant.mg)) +
  geom_point() +
  geom_line(data = PAVIsim) +
  facet_wrap(~nutrient)
  
RPVIsim <- RPVIdat %>%
  select(dpiI) %>%
  unique() %>%
  expand_grid(RPVK) %>%
  mutate(N0 = 100,
         r = 0.6,
         quant.mg = K*N0/(N0 + (K-N0) * exp(-r * dpiI)))

ggplot(RPVIdat, aes(dpiI, quant.mg)) +
  geom_point() +
  geom_line(data = RPVIsim) +
  facet_wrap(~nutrient)

# formula
form1 <- formula(quant.mg ~ K*N0/(N0 + (K-N0) * exp(-r * dpiI)) | nutrient)

# model
(PAVmod1 <- nlsList(form1, data = PAVIdat, 
                    start = list(N0 = 10, r = 0.1)))

(RPVmod1 <- nlsList(form1, data = RPVIdat, 
                    start = list(N0 = 100, r = 0.6)))

#### start here ####
# fit the below models
# AIC or something to pick best one

# parameters vary by N level (2)
# parameters vary by P level (2)
# parameters don't vary