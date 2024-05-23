#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse) # version 2.0.0
library(patchwork) # version 1.1.2
library(glmmTMB) # version 1.1.8
library(DHARMa) # version 0.4.6

# import data
dat <- read_csv("intermediate-data/qPCR_expt_data_cleaned.csv")

# load figure settings
source("code/model_settings.R")

# function for exporting model summary
mod_sum <- function(mod, filename){
  
  # fixed effect coefficients
  mod_coef <- data.frame(coef(summary(mod))$cond) %>%
    mutate(Variable = row.names(coef(summary(mod))$cond))
  colnames(mod_coef) <- c("Estimate", "Std.error", "Z", "P", "Variable") 
  mod_coef <- mod_coef %>%
    relocate(Variable)
  
  # random effects
  mod_ran <- data.frame(lme4::formatVC(summary(mod)$varcor$cond))
  
  # sample size
  mod_n <- tibble(N = as.numeric(summary(mod)$nobs))
  
  # output
  write.table(mod_coef, paste0("output/", filename,".csv"), col.names = T, row.names = F, sep=",")
  write.table(mod_ran, paste0("output/", filename,".csv"), col.names = T, row.names = F, sep=",", append=TRUE)
  write.table(mod_n, paste0("output/", filename,".csv"), col.names = T, row.names = F, sep=",", append=TRUE)

}


#### edit data ####

# make wide by target
# for each invasion trial, determine whether resident established
# for each single inoculation, determine whether there was contamination
# identify samples missing one of the viruses (not in dataset because of contamination, above standard curve, etc.)
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
         PAV_role = case_when(invasion == "S" & first_inoculation == "PAV" ~ "only",
                              invasion == "I" & first_inoculation == "PAV" ~ "resident",
                              invasion == "I" & first_inoculation == "RPV" ~ "invader",
                              TRUE ~ NA_character_),
         RPV_role = case_when(invasion == "S" & first_inoculation == "RPV" ~ "only",
                              invasion == "I" & first_inoculation == "RPV" ~ "resident",
                              invasion == "I" & first_inoculation == "PAV" ~ "invader",
                              TRUE ~ NA_character_),
         dpiUI = case_when(invasion == "S" ~ dpiR, # DPI that can be used for single and invaders
                           invasion == "I" ~ dpiI)) # incorrect for residents

# samples with each issue
dat2 %>%
  filter(missing_quant == 1)
# 56 missing one of the viruses

dat2 %>%
  filter(invasion == "I" & resident_est == 0) %>%
  count(first_inoculation)
# 32 failed establishments

dat2 %>%
  filter(invasion == "S") %>%
  count(first_inoculation)
# 197 single infections, 100 PAV

dat2 %>%
  filter(single_cont == 1) %>%
  count(first_inoculation)
# 92 contaminated single infections, all in PAV

dat2 %>%
  filter(invasion == "S") %>%
  ggplot(aes(x = log10(RPV_quant.mg), color = as.factor(single_cont))) +
  geom_density()
# clear separation between contamination and intentional inoculation

dat2 %>%
  filter(invasion == "S") %>%
  ggplot(aes(x = log10(PAV_quant.mg), color = as.factor(single_cont))) +
  geom_density()
# no clear separation

# save single inoculation contamination
dat2 %>%
  filter(single_cont == 1) %>%
  count(first_inoculation, dpiR, nutrient) %>%
  pivot_wider(names_from = dpiR, values_from = n) %>%
  write_csv("output/single_inoculation_contamination.csv")

# remove failed resident establishment
# contaminated single inoculations
# and missing quantities
dat3 <- dat2 %>%
  filter(!(invasion == "I" & resident_est == 0) & 
           !(invasion == "S" & single_cont == 1) &
           missing_quant != 1)

# PAV quantities visualization
dat3 %>%
  filter(!is.na(PAV_role)) %>%
  ggplot(aes(x = time, y = PAV_quant.mg, color = PAV_role)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  facet_wrap(~ nutrient)
# not enough remaining single PAV across treatments

# PAV single infections
(PAV_single_summary <- dat2 %>%
  filter(PAV_role == "only") %>%
  summarise(missing = sum(missing_quant),
            contaminated = sum(single_cont),
            total = n(),
            maxRPV = max(RPV_quant.mg, na.rm = T)) %>%
  mutate(log_maxRPV = log10(maxRPV)))
# 92 out of 100 had some RPV
# max RPV is 171937, consistent with density figure

# RPV with NA titer
(init_RPV_NA <- filter(dat2, is.na(RPV_quant.mg)) %>% 
  count(RPV_role, PAV_role))

# RPV with NA titer after lowering the threshold
filter(dat2, is.na(RPV_quant.mg) | RPV_quant.mg < PAV_single_summary$maxRPV) %>% 
  count(RPV_role, PAV_role) %>%
  left_join(init_RPV_NA %>%
              rename(n_init = n)) %>%
  mutate(n_new = n - n_init)
# 20 new added

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
# 6

# remove failed resident establishment
# rename nutrient levels
dat5 <- dat4 %>%
  filter(!(invasion == "I" & resident_est == 0)) %>%
  mutate(Nutrient = fct_recode(nutrient, "low" = "L", "+N+P" = "NP", "+N" = "N", "+P" = "P") %>%
           fct_relevel("low", "+N", "+P"))

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

# block-level replication
dat5 %>%
  distinct(set, tube_label) %>%
  count(set)

# divide datasets by virus
PAVdat <- dat5 %>%
  filter(!is.na(PAV_role) & !is.na(PAV_quant.mg)) %>%
  mutate(quant.mg = PAV_quant.mg,
         log_quant.mg = log(quant.mg + 1),
         role = PAV_role,
         quant.plant = quant.mg * shoot_mass.g,
         log_quant.plant = log(quant.plant + 1))

RPVdat <- dat5 %>%
  filter(!is.na(RPV_role) & !is.na(RPV_quant.mg)) %>%
  mutate(quant.mg = RPV_quant.mg,
         log_quant.mg = log(quant.mg + 1),
         role = RPV_role,
         quant.plant = quant.mg * shoot_mass.g,
         log_quant.plant = log(quant.plant + 1))

# set for invasions
bioInvDat <- dat5 %>%
  filter(invasion == "I") %>%
  inner_join(full_join(PAVdat, RPVdat) %>%
               filter(invasion == "I") %>%
               distinct(sample))

# set for residents of invasion & single infections
bioResDat <- dat5 %>%
  inner_join(full_join(PAVdat, RPVdat) %>%
               distinct(sample)) %>%
  mutate(inoculation = paste(invasion, first_inoculation, sep = "_"))

ggplot(bioInvDat, aes(x = dpp, y = shoot_mass.g, color = nutrient)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_se") +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ first_inoculation)

ggplot(bioResDat, aes(x = dpp, y = shoot_mass.g, color = nutrient)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_se") +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ inoculation)


#### biomass models ####

# nutrients, established virus
bioInvMod <- glmmTMB(shoot_mass.g ~ highN * highP * first_inoculation + (1|set) + (1|time), 
                     data = bioInvDat)
summary(bioInvMod)
bioInvResid <- simulateResiduals(bioInvMod, n = 1000)
plot(bioInvResid)

bioResMod <- glmmTMB(shoot_mass.g ~ highN * highP * inoculation + (1|set) + (1|time), 
                     data = bioResDat)
summary(bioResMod)
bioResResid <- simulateResiduals(bioResMod, n = 1000)
plot(bioResResid)

# Table 4
mod_sum(bioInvMod, "invasion_biomass_model")

# Table S11
mod_sum(bioResMod, "resident_biomass_model")

# values for text
# N effect
noNBio <- fixef(bioInvMod)$cond[1]
NBio <- fixef(bioInvMod)$cond[1] + fixef(bioInvMod)$cond[2]
100 * (NBio - noNBio) / noNBio


#### invasion models ####

# subset data
PAVIdat <- PAVdat %>%
  filter(role == "invader")

RPVIdat <- RPVdat %>%
  filter(role == "invader")

# model with titer concentration
PAVImod <- glmmTMB(log_quant.mg ~ dpiI + highN:dpiI + highP:dpiI + highN:highP:dpiI + (1|set), 
                    data = PAVIdat)
summary(PAVImod)
PAVIResid <- simulateResiduals(PAVImod, n = 1000)
plot(PAVIResid)
testZeroInflation(PAVIResid)

PAVImodb <- glmmTMB(quant.mg ~ dpiI + highN:dpiI + highP:dpiI + highN:highP:dpiI + (1|set), 
                    data = PAVIdat, family = "tweedie")
summary(PAVImodb)
PAVIResidb <- simulateResiduals(PAVImodb, n = 1000)
plot(PAVIResidb)

RPVImod <- glmmTMB(log_quant.mg ~ dpiI + highN:dpiI + highP:dpiI + highN:highP:dpiI + (1|set), 
                   data = RPVIdat)
summary(RPVImod)
RPVIResid <- simulateResiduals(RPVImod, n = 1000)
plot(RPVIResid)

RPVImodb <- glmmTMB(quant.mg ~ dpiI + highN:dpiI + highP:dpiI + highN:highP:dpiI + (1|set), 
                    data = RPVIdat, family = "lognormal")
summary(RPVImodb)
RPVIResidb <- simulateResiduals(RPVImodb, n = 1000)
plot(RPVIResidb)
plotResiduals(RPVIResidb, RPVIdat$dpiI)
plotResiduals(RPVIResidb, as.factor(RPVIdat$nutrient))
plotResiduals(RPVIResidb, as.factor(paste(RPVIdat$dpiI, RPVIdat$nutrient)))
plotResiduals(RPVIResidb, as.factor(RPVIdat$set))

# use log-link models
PAVImod <- PAVImodb
RPVImod <- RPVImodb

# model with plant-scale titer
PAVImod2 <- glmmTMB(quant.plant ~ dpiI + highN:dpiI + highP:dpiI + highN:highP:dpiI + (1|set), 
                    data = PAVIdat, family = "tweedie")
summary(PAVImod2)
PAVIResid2 <- simulateResiduals(PAVImod2, n = 1000)
plot(PAVIResid2)

RPVImod2 <- glmmTMB(quant.plant ~ dpiI + highN:dpiI + highP:dpiI + highN:highP:dpiI + (1|set), 
                    data = RPVIdat, family = "lognormal")
summary(RPVImod2)
RPVIResid2 <- simulateResiduals(RPVImod2, n = 1000)
plot(RPVIResid2)

# visualize model 1
PAVIsim1 <- PAVIdat %>%
  select(nutrient, highN, highP, dpiI) %>%
  unique() %>%
  mutate(set = NA) %>%
  mutate(quant.mg = predict(PAVImod, newdata = ., re.form = NA, type = "response"),
         log_quant.mg = predict(PAVImod, newdata = ., re.form = NA))

ggplot(PAVIdat, aes(dpiI, quant.mg)) +
  geom_point() +
  geom_line(data = PAVIsim1) +
  facet_wrap(~nutrient)

ggplot(PAVIdat, aes(dpiI, log_quant.mg)) +
  geom_point() +
  geom_line(data = PAVIsim1) +
  facet_wrap(~nutrient)

RPVIsim1 <- RPVIdat %>%
  select(nutrient, highN, highP, dpiI) %>%
  unique() %>%
  mutate(set = NA) %>%
  mutate(quant.mg = predict(RPVImod, newdata = ., re.form = NA, type = "response"),
         log_quant.mg = predict(RPVImod, newdata = ., re.form = NA))

ggplot(RPVIdat, aes(dpiI, quant.mg)) +
  geom_point() +
  geom_line(data = RPVIsim1) +
  facet_wrap(~nutrient)

ggplot(RPVIdat, aes(dpiI, log_quant.mg)) +
  geom_point() +
  geom_line(data = RPVIsim1) +
  facet_wrap(~nutrient)

# visualize model 2
PAVIsim2 <- PAVIdat %>%
  select(nutrient, highN, highP, dpiI) %>%
  unique() %>%
  mutate(set = NA) %>%
  mutate(quant.plant = predict(PAVImod2, newdata = ., re.form = NA, type = "response"),
         log_quant.plant = predict(PAVImod2, newdata = ., re.form = NA))

ggplot(PAVIdat, aes(dpiI, quant.plant)) +
  geom_point() +
  geom_line(data = PAVIsim2) +
  facet_wrap(~nutrient)

ggplot(PAVIdat, aes(dpiI, log_quant.plant)) +
  geom_point() +
  geom_line(data = PAVIsim2) +
  facet_wrap(~nutrient)

RPVIsim2 <- RPVIdat %>%
  select(nutrient, highN, highP, dpiI) %>%
  unique() %>%
  mutate(set = NA) %>%
  mutate(quant.plant = predict(RPVImod2, newdata = ., re.form = NA, type = "response"),
         log_quant.plant = predict(RPVImod2, newdata = ., re.form = NA))

ggplot(RPVIdat, aes(dpiI, quant.plant)) +
  geom_point() +
  geom_line(data = RPVIsim2) +
  facet_wrap(~nutrient)

ggplot(RPVIdat, aes(dpiI, log_quant.plant)) +
  geom_point() +
  geom_line(data = RPVIsim2) +
  facet_wrap(~nutrient)

# invasion growth rates
fixef(PAVImod)$cond[2]
fixef(RPVImod)$cond[2]

# initial titer values
exp(fixef(PAVImod)$cond[1])
exp(fixef(RPVImod)$cond[1])

# Table 2
mod_sum(PAVImod, "pav_invasion_model")

# Table 3
mod_sum(RPVImod, "rpv_invasion_model")

# Table S7
mod_sum(PAVImod2, "pav_invasion_per_plant_model")

# Table S8
mod_sum(RPVImod2, "rpv_invasion_per_plant_model")


#### compare resident to single ####

# select data
filter(PAVdat, role %in% c("only", "resident") & quant.mg < 100) %>% select(quant.mg) # only one quant = 0
filter(PAVdat, role %in% c("only", "resident") & quant.plant < 100) %>% select(quant.plant) 
# can't fit log-normal with quant = 0 and choice of + # affects model fit

PAVURdat <- PAVdat %>%
  filter(role %in% c("only", "resident") & quant.mg > 0) %>%
  mutate(dpiUR = dpiR - min(dpiR),
         role = fct_relevel(role, "resident"))

RPVURdat <- RPVdat %>%
  filter(role %in% c("only", "resident")) %>%
  mutate(dpiUR = dpiR - min(dpiR),
         role = fct_relevel(role, "resident"))

# visualize
ggplot(PAVURdat, aes(dpiUR, log_quant.mg, color = nutrient)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ role)
# role affects intercept and linearity

ggplot(PAVURdat, aes(dpiUR, quant.mg, color = nutrient)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ role)

ggplot(PAVURdat, aes(dpiUR, log_quant.plant, color = nutrient)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ role)

ggplot(RPVURdat, aes(dpiUR, log_quant.mg, color = nutrient)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ role)
# nutrient affects peak

ggplot(RPVURdat, aes(dpiUR, log_quant.plant, color = nutrient)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ role)

# polynomial model
PAVURmod1 <- glmmTMB(quant.mg ~ highN * highP * role  * (dpiUR + I(dpiUR^2)) + (1|set), 
                     data = PAVURdat, family = "lognormal")
# warning okay if no convergence warning
summary(PAVURmod1)
PAVURResid1 <- simulateResiduals(PAVURmod1, n = 1000)
plot(PAVURResid1)

RPVURmod1 <- glmmTMB(quant.mg ~ highN * highP * role  * (dpiUR + I(dpiUR^2)) + (1|set), 
                     data = RPVURdat, family = "lognormal")
# warning okay if no convergence warning
summary(RPVURmod1)
RPVURResid1 <- simulateResiduals(RPVURmod1, n = 1000)
plot(RPVURResid1)

# linear model
PAVURmod2 <- glmmTMB(quant.mg ~ highN * highP * role  * dpiUR + (1|set), 
                     data = PAVURdat, family = "lognormal")
summary(PAVURmod2)
PAVURResid2 <- simulateResiduals(PAVURmod2, n = 1000)
plot(PAVURResid2)

RPVURmod2 <- glmmTMB(quant.mg ~ highN * highP * role  * dpiUR + (1|set), 
                     data = RPVURdat, family = "lognormal")
summary(RPVURmod2)
RPVURResid2 <- simulateResiduals(RPVURmod2, n = 1000)
plot(RPVURResid2)

# compare models
AIC(PAVURmod1, PAVURmod2) # linear
AIC(PAVURmod1) - AIC(PAVURmod2)
AIC(RPVURmod1, RPVURmod2) # polynomial
AIC(RPVURmod1) - AIC(RPVURmod2)

# check model fit
PAVURdat %>%
  select(role, Nutrient, highN, highP) %>%
  expand_grid(tibble(dpiR = min(PAVURdat$dpiR):max(PAVURdat$dpiR),
                     dpiUR = min(PAVURdat$dpiUR):max(PAVURdat$dpiUR))) %>%
  unique() %>%
  mutate(set = NA) %>%
  mutate(quant.mg = predict(PAVURmod2, newdata = ., re.form = NA, type = "response"),
         quant.se = predict(PAVURmod2, newdata = ., se.fit = T, re.form = NA, type = "response")$se.fit) %>%
  ggplot(aes(dpiR, quant.mg, color = Nutrient, fill = Nutrient, linetype = Nutrient)) +
  geom_ribbon(aes(ymin = quant.mg - quant.se, ymax = quant.mg + quant.se), 
              color = NA, alpha= 0.1) +
  geom_line() +
  stat_summary(data = PAVURdat, geom = "errorbar", fun.data = "mean_se", width = 0, alpha = 0.8, 
               position = position_dodge(dodge_size), linetype = "solid") +
  stat_summary(data = PAVURdat, geom = "point", size = 2, fun = "mean", position = position_dodge(dodge_size),
               shape = 21, color = "black") +
  facet_wrap(~ role, scales = "free")

# save model
PAVURmod <- PAVURmod2
RPVURmod <- RPVURmod1

# Table S5
mod_sum(PAVURmod, "pav_established_model")

# Table S6
mod_sum(RPVURmod, "rpv_established_model")

# polynomial model with total titer
PAVURmod3 <- glmmTMB(quant.plant ~ highN * highP * role  * (dpiUR + I(dpiUR^2)) + (1|set), 
                     data = PAVURdat, family = "lognormal")
summary(PAVURmod3)
PAVURResid3 <- simulateResiduals(PAVURmod3, n = 1000)
plot(PAVURResid3)

RPVURmod3 <- glmmTMB(quant.plant ~ highN * highP * role  * (dpiUR + I(dpiUR^2)) + (1|set), 
                     data = RPVURdat, family = "lognormal")
summary(RPVURmod3)
RPVURResid3 <- simulateResiduals(RPVURmod3, n = 1000)
plot(RPVURResid3)

# linear model
PAVURmod4 <- glmmTMB(quant.plant ~ highN * highP * role  * dpiUR + (1|set), 
                     data = PAVURdat, family = "lognormal")
summary(PAVURmod4)
PAVURResid4 <- simulateResiduals(PAVURmod4, n = 1000)
plot(PAVURResid4)

RPVURmod4 <- glmmTMB(quant.plant ~ highN * highP * role  * dpiUR + (1|set), 
                     data = RPVURdat, family = "lognormal")
summary(RPVURmod4)
RPVURResid4 <- simulateResiduals(RPVURmod4, n = 1000)
plot(RPVURResid4)

# compare models
AIC(PAVURmod3, PAVURmod4) # linear
AIC(PAVURmod3) - AIC(PAVURmod4)
AIC(RPVURmod3, RPVURmod4) # polynomial
AIC(RPVURmod3) - AIC(RPVURmod4)

# Table S9
mod_sum(PAVURmod4, "pav_established_per_plant_model")

# Table S10
mod_sum(RPVURmod3, "rpv_established_per_plant_model")


#### visualize ####

# dodge points
dodge_size = 2

### invasion figures ###

PAVISimFig <- PAVIdat %>%
  distinct(role, Nutrient, highN, highP) %>%
  expand_grid(dpiI = min(PAVIdat$dpiI):max(PAVIdat$dpiI)) %>%
  mutate(set = NA) %>%
  mutate(quant.mg = predict(PAVImod, newdata = ., re.form = NA, type = "response"),
         quant.se = predict(PAVImod, newdata = ., se.fit = T, re.form = NA, type = "response")$se.fit)

PAVIfig <- ggplot(PAVIdat, aes(x = dpiI, y = quant.mg, color = Nutrient, fill = Nutrient, linetype = Nutrient)) +
  geom_ribbon(data = PAVISimFig, aes(ymin = quant.mg - quant.se, ymax = quant.mg + quant.se), 
              color = NA, alpha = 0.1, show.legend = F) +
  geom_line(data = PAVISimFig) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, alpha = 0.8, 
               position = position_dodge(dodge_size), linetype = "solid", show.legend = F) +
  stat_summary(geom = "point", size = 2, fun = "mean", position = position_dodge(dodge_size),
               shape = 21, color = "black") +
  scale_color_viridis_d(direction = -1, name = "Nutrient\nsupply") +
  scale_fill_viridis_d(direction = -1, name = "Nutrient\nsupply") +
  scale_linetype(name = "Nutrient\nsupply") +
  scale_y_continuous(labels=function(x)x/1e3) +
  labs(y = bquote("BYDV-PAV titer (thousands "~mg^-1~")"), title = "Invading virus") +
  fig_theme +
  theme(legend.position = c(0.2, 0.6),
        legend.key.width = unit(1, "cm"),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

RPVISimFig <- RPVIdat %>%
  select(role, Nutrient, highN, highP) %>%
  expand_grid(dpiI = min(RPVIdat$dpiI):max(RPVIdat$dpiI)) %>%
  unique() %>%
  mutate(set = NA) %>%
  mutate(quant.mg = predict(RPVImod, newdata = ., re.form = NA, type = "response"),
         quant.se = predict(RPVImod, newdata = ., se.fit = T, re.form = NA, type = "response")$se.fit)

RPVIfig <- ggplot(RPVIdat, aes(x = dpiI, y = quant.mg, color = Nutrient, 
                               fill = Nutrient, linetype = Nutrient)) +
  geom_ribbon(data = RPVISimFig, aes(ymin = quant.mg - quant.se, ymax = quant.mg + quant.se), 
              color = NA, alpha = 0.1, show.legend = F) +
  geom_line(data = RPVISimFig) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, alpha = 0.8, 
               position = position_dodge(dodge_size), linetype = "solid") +
  stat_summary(geom = "point", size = 2, fun = "mean", position = position_dodge(dodge_size),
               shape = 21, color = "black") +
  scale_color_viridis_d(direction = -1) +
  scale_fill_viridis_d(direction = -1) +
  scale_y_continuous(labels=function(x)x/1e3) +
  labs(x = "Days post invader inoculation", 
       y = bquote("CYDV-RPV titer (thousands "~mg^-1~")")) +
  fig_theme +
  theme(legend.position = "none")

### resident figures ###

PAVRSimFig <- PAVURdat %>%
  filter(role == "resident") %>%
  select(role, Nutrient, highN, highP) %>%
  expand_grid(tibble(dpiR = min(PAVURdat$dpiR):max(PAVURdat$dpiR),
                     dpiUR = min(PAVURdat$dpiUR):max(PAVURdat$dpiUR))) %>%
  unique() %>%
  mutate(set = NA) %>%
  mutate(quant.mg = predict(PAVURmod, newdata = ., re.form = NA, type = "response"),
         quant.se = predict(PAVURmod, newdata = ., se.fit = T, re.form = NA, type = "response")$se.fit)

PAVRfig <- PAVURdat %>%
  filter(role == "resident") %>%
  ggplot(aes(dpiR, quant.mg, color = Nutrient, fill = Nutrient, linetype = Nutrient)) +
  geom_ribbon(data = PAVRSimFig, aes(ymin = quant.mg - quant.se, ymax = quant.mg + quant.se), 
              color = NA, alpha= 0.1) +
  geom_line(data = PAVRSimFig,) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, alpha = 0.8, 
               position = position_dodge(dodge_size), linetype = "solid") +
  stat_summary(geom = "point", size = 2, fun = "mean", position = position_dodge(dodge_size),
               shape = 21, color = "black") +
  scale_color_viridis_d(direction = -1) +
  scale_fill_viridis_d(direction = -1) +
  scale_y_continuous(labels=function(x)x/1e3) +
  labs(title = "Resident virus") +
  fig_theme +
  theme(legend.position = "none",
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

RPVRSimFig <- RPVURdat %>%
  filter(role == "resident") %>%
  select(role, Nutrient, highN, highP) %>%
  expand_grid(tibble(dpiR = min(RPVURdat$dpiR):max(RPVURdat$dpiR),
                     dpiUR = min(RPVURdat$dpiUR):max(RPVURdat$dpiUR))) %>%
  unique() %>%
  mutate(set = NA) %>%
  mutate(quant.mg = predict(RPVURmod, newdata = ., re.form = NA, type = "response"),
         quant.se = predict(RPVURmod, newdata = ., se.fit = T, re.form = NA, type = "response")$se.fit)

RPVRfig <- RPVdat %>%
  filter(role == "resident") %>%
  ggplot(aes(dpiR, quant.mg, color = Nutrient, fill = Nutrient, linetype = Nutrient)) +
  geom_ribbon(data = RPVRSimFig, aes(ymin = quant.mg - quant.se, ymax = quant.mg + quant.se), 
              color = NA, alpha = 0.1, show.legend = F) +
  geom_line(data = RPVRSimFig) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, alpha = 0.8, 
               position = position_dodge(dodge_size), linetype = "solid") +
  stat_summary(geom = "point", size = 2, fun = "mean", position = position_dodge(dodge_size),
               shape = 21, color = "black") +
  scale_color_viridis_d(direction = -1) +
  scale_fill_viridis_d(direction = -1) +
  scale_y_continuous(labels=function(x)x/1e3) +
  labs(x = "Days post resident inoculation") +
  fig_theme +
  theme(legend.position = "none",
        axis.title.y = element_blank())

# combine figures
comb_fig <- PAVIfig + theme(axis.text.x = element_blank()) + 
  PAVRfig + theme(axis.text.x = element_blank()) + 
  RPVIfig + RPVRfig + 
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")") & 
  theme(plot.tag = element_text(size = 8, face = "bold"))

# Figure 4
ggsave("output/PAV_RPV_virus_titer_figure.pdf", comb_fig,
       width = 18, height = 13.5, units = "cm")


#### sample size table ####

# PAV sample sizes
pav_size <- PAVdat %>%
  filter(role == "invader" | (role %in% c("only", "resident") & quant.mg > 0)) %>%
  count(role, dpiI, nutrient) %>%
  pivot_wider(names_from = "role",
              values_from = "n",
              names_prefix = "pav_")

# RPV sample sizes
rpv_size <- RPVdat %>%
  count(role, dpiI, nutrient) %>%
  pivot_wider(names_from = "role",
              values_from = "n",
              names_prefix = "rpv_")

# combine viruses
virus_size <- full_join(pav_size, rpv_size)

# can we condense columns?
filter(virus_size, pav_invader != rpv_resident)
filter(virus_size, rpv_invader != pav_resident) # no

# Table S4
write_csv(virus_size, "output/PAV_RPV_sample_sizes.csv")
