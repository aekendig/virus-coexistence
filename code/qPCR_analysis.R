## Goal: analyze qPCR data


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)
library(glmmTMB)

# import data
dat <- read_csv("intermediate-data/qPCR_expt_data_cleaned.csv")

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
  filter(invasion == "I" & resident_est == 0) 
# 32 failed establishments

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

# remove failed resident establishment
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
# little single PAV data
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
  filter(PAV_role == "only") %>%
  summarise(missing = sum(missing_quant),
            contaminated = sum(single_cont),
            total = n(),
            maxRPV = max(RPV_quant.mg, na.rm = T)) %>%
  mutate(log_maxRPV = log10(maxRPV)))
# 92 out of 100 had some RPV
# max RPV is 171937, consistent with density figure

# same for RPV?
dat2 %>%
  filter(RPV_role == "only") %>%
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

# remove failed resident establishment
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
  mutate(quant.mg = PAV_quant.mg,
         log_quant.mg = log(quant.mg + 1),
         log10_quant.mg = log10(quant.mg + 1),
         role = PAV_role,
         rel_quant.mg = quant.mg/max(quant.mg, na.rm = T),
         rel_log_quant.mg = log_quant.mg/max(log_quant.mg),
         quant.plant = quant.mg * shoot_mass.g,
         log_quant.plant = log(quant.plant + 1))

RPVdat <- dat5 %>%
  filter(!is.na(RPV_role) & !is.na(RPV_quant.mg)) %>%
  mutate(quant.mg = RPV_quant.mg,
         log_quant.mg = log(quant.mg + 1),
         log10_quant.mg = log10(quant.mg + 1),
         role = RPV_role,
         rel_quant.mg = quant.mg/max(quant.mg, na.rm = T),
         rel_log_quant.mg = log_quant.mg/max(log_quant.mg),
         quant.plant = quant.mg * shoot_mass.g,
         log_quant.plant = log(quant.plant + 1))

# set for nutrients
bioInvDat <- dat5 %>%
  filter(invasion == "I") %>%
  mutate(log_shoot_mass.g = log(shoot_mass.g))

# set for residents
bioResDat <- dat5 %>%
  mutate(inoculation = paste(invasion, first_inoculation, sep = "_"),
         log_shoot_mass.g = log(shoot_mass.g))

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


#### biomass model ####

# nutrients, established virus
bioInvMod <- glmmTMB(log_shoot_mass.g ~ highN * highP * first_inoculation + (1|time), data = bioInvDat)
# model won't converge with set as a random effect
summary(bioInvMod)
mod_sum(bioInvMod, "invasion_biomass_model")

bioResMod <- glmmTMB(log_shoot_mass.g ~ highN * highP * inoculation + (1|set) + (1|time), data = bioResDat)
summary(bioResMod)
mod_sum(bioResMod, "resident_biomass_model")

# values for text

# N effect
noNBio <- exp(fixef(bioInvMod)$cond[1])
NBio <- exp(fixef(bioInvMod)$cond[1] + fixef(bioInvMod)$cond[2])
100 * (NBio - noNBio) / noNBio

# P, inoculation interactions
intcpt <- exp(fixef(bioResMod)$cond[1])
PBio <- exp(fixef(bioResMod)$cond[1] + fixef(bioResMod)$cond[3])
100 * (PBio - intcpt) / intcpt # increase

noPRPVBio <- exp(fixef(bioResMod)$cond[1] + fixef(bioResMod)$cond[6])
PRPVBio <- exp(fixef(bioResMod)$cond[1] + fixef(bioResMod)$cond[3] + fixef(bioResMod)$cond[6] + fixef(bioResMod)$cond[13])
100 * (PRPVBio - noPRPVBio) / noPRPVBio # decrease

noPPAVBio <- exp(fixef(bioResMod)$cond[1] + fixef(bioResMod)$cond[5])
PPAVBio <- exp(fixef(bioResMod)$cond[1] + fixef(bioResMod)$cond[3] + fixef(bioResMod)$cond[5] + fixef(bioResMod)$cond[12])
100 * (PPAVBio - noPPAVBio) / noPPAVBio # decrease

noPInvBio <- exp(fixef(bioResMod)$cond[1] + fixef(bioResMod)$cond[4])
PInvBio <- exp(fixef(bioResMod)$cond[1] + fixef(bioResMod)$cond[3] + fixef(bioResMod)$cond[4] + fixef(bioResMod)$cond[11])
100 * (PInvBio - noPInvBio) / noPInvBio # smaller decrease


#### estimate r and N0 parameters ####

# subset data
PAVIdat <- PAVdat %>%
  filter(role == "invader")

RPVIdat <- RPVdat %>%
  filter(role == "invader")

# simulated data/starting values
PAVIsim <- PAVIdat %>%
  select(dpiI) %>%
  unique() %>%
  mutate(log_N0 = log(1),
         r = 0.3,
         log_quant.mg = log_N0 + r * dpiI)

# visualize
ggplot(PAVIdat, aes(dpiI, log_quant.mg)) +
  geom_point() +
  geom_line(data = PAVIsim) +
  facet_wrap(~nutrient)

ggplot(PAVIdat, aes(dpiI, log_quant.plant)) +
  geom_point() +
  facet_wrap(~nutrient)

# simulated data/starting values
RPVIsim <- RPVIdat %>%
  select(dpiI) %>%
  unique() %>%
  mutate(log_N0 = 10,
         r = 0.3,
         log_quant.mg = log_N0 + r * dpiI)

# visualize
ggplot(RPVIdat, aes(dpiI, log_quant.mg)) +
  geom_point() +
  geom_line(data = RPVIsim) +
  facet_wrap(~nutrient)

ggplot(RPVIdat, aes(dpiI, log_quant.plant)) +
  geom_point() +
  facet_wrap(~nutrient)

# model with titer concentration
PAVImod <- glmmTMB(log_quant.mg ~ dpiI + highN:dpiI + highP:dpiI + highN:highP:dpiI + (1|set), data = PAVIdat)
summary(PAVImod)
RPVImod <- glmmTMB(log_quant.mg ~ dpiI + highN:dpiI + highP:dpiI + highN:highP:dpiI + (1|set), data = RPVIdat)
summary(RPVImod)

# model with plant-scale titer
PAVImod2 <- glmmTMB(log_quant.plant ~ dpiI + highN:dpiI + highP:dpiI + highN:highP:dpiI + (1|set), data = PAVIdat)
summary(PAVImod2)
RPVImod2 <- glmmTMB(log_quant.plant ~ dpiI + highN:dpiI + highP:dpiI + highN:highP:dpiI + (1|set), data = RPVIdat)
summary(RPVImod2)

# # model 2: N and P each have a different growth rate, no interaction
# PAVImod2 <- lm(log_quant.mg ~ dpiI + highN:dpiI + highP:dpiI, data = PAVIdat)
# summary(PAVImod2)
# RPVImod2 <- lm(log_quant.mg ~ dpiI + highN:dpiI + highP:dpiI, data = RPVIdat)
# summary(RPVImod2)
# 
# # model 3: only N causes a different growth rate
# PAVImod3 <- lm(log_quant.mg ~ dpiI + highN:dpiI, data = PAVIdat)
# summary(PAVImod3)
# RPVImod3 <- lm(log_quant.mg ~ dpiI + highN:dpiI, data = RPVIdat)
# summary(RPVImod3)
# 
# # model 4: only P causes a different growth rate
# PAVImod4 <- lm(log_quant.mg ~ dpiI + highP:dpiI, data = PAVIdat)
# summary(PAVImod4)
# RPVImod4 <- lm(log_quant.mg ~ dpiI + highP:dpiI, data = RPVIdat)
# summary(RPVImod4)
# 
# # model 5: all growth rates are the same
# PAVImod5 <- lm(log_quant.mg ~ dpiI, data = PAVIdat)
# summary(PAVImod5)
# RPVImod5 <- lm(log_quant.mg ~ dpiI, data = RPVIdat)
# summary(RPVImod5)
# 
# # compare
# anova(PAVImod1, PAVImod2) # no effect of losing interaction
# anova(PAVImod3, PAVImod5) # no effect of losing N
# anova(PAVImod4, PAVImod5) # no effect of losing P
# 
# anova(RPVImod1, RPVImod2) # no effect of losing interaction
# anova(RPVImod3, RPVImod5) # no effect of losing N
# anova(RPVImod4, RPVImod5) # no effect of losing P

# visualize model 1
PAVIsim1 <- PAVIdat %>%
  select(nutrient, highN, highP, dpiI) %>%
  unique() %>%
  mutate(set = NA) %>%
  mutate(log_quant.mg = predict(PAVImod, newdata = ., re.form = NA),
         quant.mg = exp(log_quant.mg))

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
  mutate(quant.mg = exp(predict(RPVImod, newdata = ., re.form = NA)))

ggplot(RPVIdat, aes(dpiI, quant.mg)) +
  geom_point() +
  geom_line(data = RPVIsim1) +
  facet_wrap(~nutrient)

# visualize model 2
PAVIsim2 <- PAVIdat %>%
  select(nutrient, highN, highP, dpiI) %>%
  unique() %>%
  mutate(set = NA) %>%
  mutate(log_quant.plant = predict(PAVImod2, newdata = ., re.form = NA),
         quant.plant = exp(log_quant.plant))

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
  mutate(quant.plant = exp(predict(RPVImod2, newdata = ., re.form = NA)))

ggplot(RPVIdat, aes(dpiI, quant.plant)) +
  geom_point() +
  geom_line(data = RPVIsim2) +
  facet_wrap(~nutrient)

# r values
fixef(PAVImod)$cond[2]
fixef(RPVImod)$cond[2]

# N0 values
exp(fixef(PAVImod)$cond[1])
exp(fixef(RPVImod)$cond[1])

# marginal N effect
dpiN_RPV <- fixef(RPVImod2)$cond[2] + fixef(RPVImod2)$cond[3]
(dpiN_RPV - fixef(RPVImod2)$cond[2]) / fixef(RPVImod2)$cond[2]

# save
mod_sum(PAVImod, "pav_invasion_model")
mod_sum(RPVImod, "rpv_invasion_model")
mod_sum(PAVImod2, "pav_invasion_per_plant_model")
mod_sum(RPVImod2, "rpv_invasion_per_plant_model")


#### compare resident to single ####

# select data
PAVURdat <- PAVdat %>%
  filter(role %in% c("only", "resident")) %>%
  mutate(dpiUR = dpiR - min(dpiR))

RPVURdat <- RPVdat %>%
  filter(role %in% c("only", "resident")) %>%
  mutate(dpiUR = dpiR - min(dpiR))

# visualize
ggplot(PAVURdat, aes(dpiUR, log_quant.mg, color = nutrient)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ role)
# role affects intercept and linearity

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

# # linear model with no effects of nutrients
# PAVURmod1 <- lm(log_quant.mg ~ dpiUR, data = PAVURdat)
# summary(PAVURmod1)
# RPVURmod1 <- lm(log_quant.mg ~ dpiUR, data = RPVURdat)
# summary(RPVURmod1)
# 
# # polynomial model with no effects of nutrients
# PAVURmod2 <- lm(log_quant.mg ~ dpiUR + I(dpiUR^2), data = PAVURdat)
# summary(PAVURmod2)
# RPVURmod2 <- lm(log_quant.mg ~ dpiUR + I(dpiUR^2), data = RPVURdat)
# summary(RPVURmod2)
# 
# # compare shapes
# anova(PAVURmod1, PAVURmod2) # not different
# anova(RPVURmod1, RPVURmod2) # quadratic

# polynomial model
PAVURmod1 <- glmmTMB(log_quant.mg ~ highN * highP * role  * (dpiUR + I(dpiUR^2)) + (1|set), data = PAVURdat)
summary(PAVURmod1)
RPVURmod1 <- glmmTMB(log_quant.mg ~ highN * highP * role  * (dpiUR + I(dpiUR^2)) + (1|set), data = RPVURdat)
summary(RPVURmod1)

# linear model
PAVURmod2 <- glmmTMB(log_quant.mg ~ highN * highP * role  * dpiUR + (1|set), data = PAVURdat)
summary(PAVURmod2)
RPVURmod2 <- glmmTMB(log_quant.mg ~ highN * highP * role  * dpiUR + (1|set), data = RPVURdat)
summary(RPVURmod2)

# compare models
AIC(PAVURmod1, PAVURmod2) # linear
AIC(RPVURmod1, RPVURmod2) # polynomial

# save model
PAVURmod <- PAVURmod2
RPVURmod <- RPVURmod1

# export
mod_sum(PAVURmod, "pav_established_model")
mod_sum(RPVURmod, "rpv_established_model")

# polynomial model with total titer
PAVURmod3 <- glmmTMB(log_quant.plant ~ highN * highP * role  * (dpiUR + I(dpiUR^2)) + (1|set), data = PAVURdat)
summary(PAVURmod3)
RPVURmod3 <- glmmTMB(log_quant.plant ~ highN * highP * role  * (dpiUR + I(dpiUR^2)) + (1|set), data = RPVURdat)
summary(RPVURmod3)

# linear model
PAVURmod4 <- glmmTMB(log_quant.plant ~ highN * highP * role  * dpiUR + (1|set), data = PAVURdat)
summary(PAVURmod4)
RPVURmod4 <- glmmTMB(log_quant.plant ~ highN * highP * role  * dpiUR + (1|set), data = RPVURdat)
summary(RPVURmod4)

# compare models
AIC(PAVURmod3, PAVURmod4) # linear (by 3)
AIC(RPVURmod3, RPVURmod4) # polynomial

# export
mod_sum(PAVURmod4, "pav_established_per_plant_model")
mod_sum(RPVURmod3, "rpv_established_per_plant_model")


# drop1(PAVURmod3, test = "F")
# drop1(RPVURmod3, test = "F")
# 
# # remove 4-way interactions
# PAVURmod4 <- update(PAVURmod3, .~. -highN:highP:role:dpiUR)
# summary(PAVURmod4)
# RPVURmod4 <- update(RPVURmod3, .~. -highN:highP:role:dpiUR-highN:highP:role:I(dpiUR^2))
# summary(RPVURmod4)
# 
# drop1(PAVURmod4, test = "F")
# drop1(RPVURmod4, test = "F")
# 
# # remove 3-way interactions
# PAVURmod5 <- update(PAVURmod4, .~. -highN:highP:role-highN:highP:dpiUR-highN:role:dpiUR-highP:role:dpiUR)
# summary(PAVURmod5)
# RPVURmod5 <- update(RPVURmod4, .~. -highN:highP:role-highN:highP:dpiUR-highN:role:dpiUR-highP:role:dpiUR-
#                       highN:highP:I(dpiUR^2)-highN:role:I(dpiUR^2)-highP:role:I(dpiUR^2))
# summary(RPVURmod5)
# 
# drop1(PAVURmod5, test = "F")
# drop1(RPVURmod5, test = "F")
# 
# # remove 2-way interactions
# PAVURmod6 <- lm(log_quant.mg ~ highN + highP + role  + dpiUR, data = PAVURdat)
# summary(PAVURmod6)
# RPVURmod6 <- lm(log_quant.mg ~ highN + highP + role  + dpiUR + I(dpiUR^2), data = RPVURdat)
# summary(RPVURmod6)
# 
# drop1(PAVURmod6, test = "F") # none sig
# drop1(RPVURmod6, test = "F") # keep time
# 
# # simplified model
# PAVURmod7 <- lm(log_quant.mg ~ 1, data = PAVURdat)
# summary(PAVURmod7)
# summary(RPVURmod2)
# 
# add1(PAVURmod7, ~ highN + highP + role + dpiUR, test = "F") # no effect
# add1(RPVURmod2, ~ .+ highN + highP + role, test = "F") # no effect
# 
# #plot(PAVURmod7)
# #plot(RPVURmod2)
# 
# # PAV value
# exp(coef(PAVURmod7)[1])


#### visualize ####

# figure settings
fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing.x = unit(0,"line"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 10),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.position = "none",
        legend.key.height = unit(5, "mm"),
        legend.key.width = unit(8, "mm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        strip.placement = "outside",
        plot.title = element_text(size = 12, vjust = 0))

# dodge points
dodge_size = 2

# shapes
shape_pal = c(21, 24, 22, 23, 4, 1)

# max values
maxPAV_quant.mg <- max(PAVdat$quant.mg, na.rm = T)
maxRPV_quant.mg <- max(RPVdat$quant.mg, na.rm = T)

### invasion figure ###
legFig <- PAVdat %>%
  filter(role == "invader") %>%
  mutate(Nutrient = fct_recode(nutrient, "low" = "L", "+N+P" = "NP", "+N" = "N", "+P" = "P") %>%
           fct_relevel("low", "+N", "+P")) %>%
  ggplot(aes(dpiI, rel_quant.mg)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, alpha = 0.8, position = position_dodge(dodge_size), aes(color = Nutrient)) +
  stat_summary(geom = "point", fun = "mean", shape = 21, position = position_dodge(dodge_size), aes(fill = Nutrient)) +
  scale_color_viridis_d(end = 0.7) +
  scale_fill_viridis_d(end = 0.7) +
  theme_bw() +
  theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.background = element_blank(),
        legend.position = c(0.25, 0.6),
        legend.key.size = unit(5, "mm"))

PAVISimFig <- PAVdat %>%
  filter(role == "invader") %>%
  select(role, dpiI, nutrient, highN, highP) %>%
  unique() %>%
  mutate(set = NA,
         Nutrient = fct_recode(nutrient, "low" = "L", "+N+P" = "NP", "+N" = "N", "+P" = "P"),
         Virus = "PAV") %>%
  mutate(log_quant.mg = predict(PAVImod, newdata = ., re.form = NA),
         log_quant.se = predict(PAVImod, newdata = ., se.fit = T, re.form = NA)$se.fit,
         log10_quant.mg = log10(exp(log_quant.mg)),
         log10_quant.se = log10(exp(log_quant.se)),
         rel_quant.mg = exp(log_quant.mg)/maxPAV_quant.mg,
         rel_quant.max = exp(log_quant.mg + log_quant.se)/maxPAV_quant.mg,
         rel_quant.min = exp(log_quant.mg - log_quant.se)/maxPAV_quant.mg) %>%
  full_join(tibble(dpiI = seq(min(RPVURdat$dpiI), max(RPVURdat$dpiI), length.out = 50),
                   dpiUR = seq(min(RPVURdat$dpiUR), max(RPVURdat$dpiUR), length.out = 50)) %>%
              expand_grid(tibble(highN = c(0, 1, 0, 1),
                                 highP = c(0, 0, 1, 1),
                                 Nutrient = c("low", "+N", "+P", "+N+P"))) %>%
              mutate(role = "resident",
                     set = NA,
                     Virus = "RPV") %>%
              mutate(log_quant.mg = predict(RPVURmod, newdata = ., re.form = NA),
                     log_quant.se = predict(RPVURmod, newdata = ., se.fit = T, re.form = NA)$se.fit,
                     log10_quant.mg = log10(exp(log_quant.mg)),
                     log10_quant.se = log10(exp(log_quant.se)),
                     rel_quant.mg = exp(log_quant.mg)/maxRPV_quant.mg,
                     rel_quant.max = exp(log_quant.mg + log_quant.se)/maxRPV_quant.mg,
                     rel_quant.min = exp(log_quant.mg - log_quant.se)/maxRPV_quant.mg)) %>%
  mutate(Nutrient = fct_relevel(Nutrient, "low", "+N", "+P"))

PAVIfig <- PAVdat %>%
  filter(role == "invader") %>%
  mutate(Virus = "PAV") %>%
  full_join(RPVdat %>%
              filter(role == "resident") %>%
              mutate(Virus = "RPV")) %>%
  mutate(Nutrient = fct_recode(nutrient, "low" = "L", "+N+P" = "NP", "+N" = "N", "+P" = "P") %>%
           fct_relevel("low", "+N", "+P")) %>%
  ggplot(aes(dpiI, log10_quant.mg)) +
  geom_ribbon(data = PAVISimFig, aes(ymin = log10_quant.mg - log10_quant.se, ymax = log10_quant.mg + log10_quant.se, shape = Virus, fill = Nutrient), color = NA, alpha= 0.1) +
  geom_line(data = PAVISimFig, aes(linetype = Virus, color = Nutrient)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, alpha = 0.8, position = position_dodge(dodge_size), aes(shape = Virus, color = Nutrient)) +
  stat_summary(geom = "point", size = 2, fun = "mean", position = position_dodge(dodge_size), aes(shape = Virus, fill = Nutrient)) +
  scale_color_viridis_d(end = 0.7, guide = "none") +
  scale_fill_viridis_d(end = 0.7, guide = "none") +
  scale_shape_manual(values = shape_pal[1:2], labels = c("BYDV-PAV invader", "CYDV-RPV resident")) +
  scale_linetype_manual(values = c("solid", "dashed"), labels = c("BYDV-PAV invader", "CYDV-RPV resident")) +
  ggtitle("(A) BYDV-PAV invades") +
  xlab("Days post invader inoculation") +
  ylab(expression(paste("Virus titer (", log[10], ")", sep = ""))) +
  fig_theme +
  theme(legend.position = c(0.77, 0.75))

RPVISimFig <- RPVdat %>%
  filter(role == "invader") %>%
  select(role, dpiI, nutrient, highN, highP) %>%
  unique() %>%
  mutate(set = NA) %>%
  mutate(log_quant.mg = predict(RPVImod, newdata = ., re.form = NA),
         log_quant.se = predict(RPVImod, newdata = ., se.fit = T, re.form = NA)$se.fit,
         log10_quant.mg = log10(exp(log_quant.mg)),
         log10_quant.se = log10(exp(log_quant.se)),
         rel_quant.mg = exp(log_quant.mg)/maxRPV_quant.mg,
         rel_quant.max = exp(log_quant.mg + log_quant.se)/maxRPV_quant.mg,
         rel_quant.min = exp(log_quant.mg - log_quant.se)/maxRPV_quant.mg,
         Virus = "RPV") %>%
  full_join(PAVdat %>%
              filter(role == "resident") %>%
              mutate(dpiUR = dpiR - min(dpiR)) %>%
              select(role, dpiI, dpiUR, nutrient, highN, highP) %>%
              unique() %>%
              mutate(set = NA) %>%
              mutate(log_quant.mg = predict(PAVURmod, newdata = ., re.form = NA),
                     log_quant.se = predict(PAVURmod, newdata = ., se.fit = T, re.form = NA)$se.fit,
                     log10_quant.mg = log10(exp(log_quant.mg)),
                     log10_quant.se = log10(exp(log_quant.se)),
                     rel_quant.mg = exp(log_quant.mg)/maxPAV_quant.mg,
                     rel_quant.max = exp(log_quant.mg + log_quant.se)/maxPAV_quant.mg,
                     rel_quant.min = exp(log_quant.mg - log_quant.se)/maxPAV_quant.mg,
                     Virus = "PAV")) %>%
  mutate(Nutrient = fct_recode(nutrient, "low" = "L", "+N+P" = "NP", "+N" = "N", "+P" = "P") %>%
           fct_relevel("low", "+N", "+P"))

RPVIfig <- RPVdat %>%
  filter(role == "invader") %>%
  mutate(Virus = "RPV") %>%
  full_join(PAVdat %>%
              filter(role == "resident") %>%
              mutate(Virus = "PAV")) %>%
  mutate(Nutrient = fct_recode(nutrient, "low" = "L", "+N+P" = "NP", "+N" = "N", "+P" = "P") %>%
           fct_relevel("low", "+N", "+P")) %>%
  ggplot(aes(dpiI, log10_quant.mg)) +
  geom_ribbon(data = RPVISimFig, aes(ymin = log10_quant.mg - log10_quant.se, ymax = log10_quant.mg + log10_quant.se, shape = Virus, fill = Nutrient), color = NA, alpha= 0.1) +
  geom_line(data = RPVISimFig, aes(linetype = Virus, color = Nutrient)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, alpha = 0.8, position = position_dodge(dodge_size), aes(shape = Virus, color = Nutrient)) +
  stat_summary(geom = "point", size = 2, fun = "mean", position = position_dodge(dodge_size), aes(shape = Virus, fill = Nutrient)) +
  scale_color_viridis_d(end = 0.7, guide = "none") +
  scale_fill_viridis_d(end = 0.7, guide = "none") +
  scale_shape_manual(values = shape_pal[3:4], labels = c("BYDV-PAV resident", "CYDV-RPV invader")) +
  scale_linetype_manual(values = c("solid", "dashed"), labels = c("BYDV-PAV resident", "CYDV-RPV invader")) +
  ggtitle("(B) CYDV-RPV invades") +
  xlab("Days post invader inoculation") +
  fig_theme +
  theme(axis.title.y = element_blank(),
        legend.position = c(0.77, 0.65))


### resident figure ###

RPVRSimFig <- tibble(dpiR = seq(min(RPVURdat$dpiR), max(RPVURdat$dpiR), length.out = 50),
                     dpiUR = seq(min(RPVURdat$dpiUR), max(RPVURdat$dpiUR), length.out = 50)) %>%
  expand_grid(tibble(highN = c(0, 1, 0, 1),
                     highP = c(0, 0, 1, 1),
                     Nutrient = c("low", "+N", "+P", "+N+P") %>%
                       fct_relevel("low", "+N", "+P")) %>%
                expand_grid(tibble(role = c("resident", "only")))) %>%
  mutate(set = NA) %>%
  mutate(log_quant.mg = predict(RPVURmod, newdata = ., re.form = NA),
         log_quant.se = predict(RPVURmod, newdata = ., se.fit = T, re.form = NA)$se.fit,
         log10_quant.mg = log10(exp(log_quant.mg)),
         log10_quant.se = log10(exp(log_quant.se)),
         rel_quant.mg = exp(log_quant.mg)/maxRPV_quant.mg,
         rel_quant.max = exp(log_quant.mg + log_quant.se)/maxRPV_quant.mg,
         rel_quant.min = exp(log_quant.mg - log_quant.se)/maxRPV_quant.mg)

RPVRfig <- RPVdat %>%
  filter(role != "invader") %>%
  mutate(Nutrient = fct_recode(nutrient, "low" = "L", "+N+P" = "NP", "+N" = "N", "+P" = "P") %>%
           fct_relevel("low", "+N", "+P")) %>%
  ggplot(aes(dpiR, log10_quant.mg, group = interaction(role, Nutrient), shape = role, linetype = role)) +
  geom_ribbon(data = RPVRSimFig, aes(ymin = log10_quant.mg - log10_quant.se, ymax = log10_quant.mg + log10_quant.se, fill = Nutrient), color = NA, alpha= 0.1, show.legend = F) +
  geom_line(data = RPVRSimFig, aes(color = Nutrient)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, alpha = 0.8, position = position_dodge(dodge_size), linetype = "solid", aes(color = Nutrient)) +
  stat_summary(geom = "point", size = 2, fun = "mean", position = position_dodge(dodge_size), aes(fill = Nutrient)) +
  scale_color_viridis_d(end = 0.7, guide = "none") +
  scale_fill_viridis_d(end = 0.7, guide = "none") +
  scale_shape_manual(values = shape_pal[c(2,5)], name = "Invasion", labels = c("CYDV-RPV resident", "CYDV-RPV only")) +
  scale_linetype_manual(values = c("dashed", "dotted"), name = "Invasion", labels = c("CYDV-RPV resident", "CYDV-RPV only")) +
  ggtitle("(C) Established CYDV-RPV") +
  xlab("Days post resident inoculation") +
  ylab(expression(paste("Virus titer (", log[10], ")", sep = ""))) +
  fig_theme +
  theme(legend.position = c(0.3, 0.15))

PAVRSimFig <- PAVURdat %>%
  select(role, dpiR, dpiUR, nutrient, highN, highP) %>%
  unique() %>%
  mutate(set = NA) %>%
  mutate(log_quant.mg = predict(PAVURmod, newdata = ., re.form = NA),
         log_quant.se = predict(PAVURmod, newdata = ., se.fit = T, re.form = NA)$se.fit,
         log10_quant.mg = log10(exp(log_quant.mg)),
         log10_quant.se = log10(exp(log_quant.se)),
         rel_quant.mg = exp(log_quant.mg)/maxPAV_quant.mg,
         rel_quant.max = exp(log_quant.mg + log_quant.se)/maxPAV_quant.mg,
         rel_quant.min = exp(log_quant.mg - log_quant.se)/maxPAV_quant.mg) %>%
  mutate(Nutrient = fct_recode(nutrient, "low" = "L", "+N+P" = "NP", "+N" = "N", "+P" = "P") %>%
           fct_relevel("low", "+N", "+P"))

# figures
PAVRfig <- PAVdat %>%
  filter(role != "invader") %>%
  mutate(Nutrient = fct_recode(nutrient, "low" = "L", "+N+P" = "NP", "+N" = "N", "+P" = "P") %>%
           fct_relevel("low", "+N", "+P")) %>%
  ggplot(aes(dpiR, log10_quant.mg, group = interaction(role, Nutrient), shape = role, linetype = role)) +
  geom_ribbon(data = PAVRSimFig, aes(ymin = log10_quant.mg - log10_quant.se, ymax = log10_quant.mg + log10_quant.se, fill = Nutrient), color = NA, alpha= 0.1, show.legend = F) +
  geom_line(data = PAVRSimFig, aes(color = Nutrient)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, alpha = 0.8, position = position_dodge(dodge_size), linetype = "solid", aes(color = Nutrient)) +
  stat_summary(geom = "point", size = 2, fun = "mean", position = position_dodge(dodge_size), aes(fill = Nutrient)) +
  scale_color_viridis_d(end = 0.7, guide = "none") +
  scale_fill_viridis_d(end = 0.7, guide = "none") +
  scale_shape_manual(values = shape_pal[c(3,6)], name = "Invasion", labels = c("BYDV-PAV resident", "BYDV-PAV only")) +
  scale_linetype_manual(values = c("solid", "dotdash"), name = "Invasion", labels = c("BYDV-PAV resident", "BYDV-PAV only")) +
  ggtitle("(D) Established BYDV-PAV") +
  xlab("Days post resident inoculation") +
  fig_theme +
  theme(legend.position = c(0.75, 0.95),
        axis.title.y = element_blank())

# add nutrient legend
nut_leg <- get_legend(legFig)
PAVIfig2 <- align_plots(PAVIfig, nut_leg, axis = "l")
PAVIfig3 <- ggdraw(PAVIfig2[[1]]) + draw_plot(PAVIfig2[[2]])

# combine
tiff("output/PAV_RPV_virus_titer_figure.tiff", width = 6.5, height = 6.5, units = "in", res = 300)
plot_grid(PAVIfig3, RPVIfig, RPVRfig, PAVRfig,
          nrow = 2,
          rel_widths = c(1, 0.9, 1, 0.9))
dev.off()

# outlier
filter(PAVdat, role == "only" & nutrient == "L" & time == 3) %>%
  select(tube_label, quant.mg, expt_notes, extraction_notes)
# looked at raw data from qPCR_data_processing
# one technical replicate was > 1000, the other two were not, their average was not
# averages below 1000 (or minimum standard) were converted to zero
