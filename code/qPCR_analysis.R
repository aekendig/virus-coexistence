## Goal: analyze qPCR data


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(patchwork)
library(glmmTMB)

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
  filter(single_cont == 1) 
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

# total uninoculated controls
dat2 %>%
  filter(PAV_role == "only" | RPV_role == "only")
# 197

# RPV with NA titer
(init_RPV_NA <- filter(dat2, is.na(RPV_quant.mg)) %>% 
  count(RPV_role, PAV_role))

# RPV with NA titer after lowering the threshold
filter(dat2, is.na(RPV_quant.mg) | RPV_quant.mg < PAV_single_summary$maxRPV) %>% 
  count(RPV_role, PAV_role) %>%
  left_join(init_RPV_NA %>%
              rename(n_init = n)) %>%
  mutate(n_new = n - n_init)

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
# rename nutrient levels
dat5 <- dat4 %>%
  filter(!(invasion == "I" & resident_est == 0)) %>%
  mutate(Nutrient = fct_recode(nutrient, "low" = "L", "+N+P" = "NP", "+N" = "N", "+P" = "P") %>%
           fct_relevel("low", "+N", "+P"))
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

# block-level replication
dat5 %>%
  distinct(set, tube_label) %>%
  count(set)

# divide datasets by virus
PAVdat <- dat5 %>%
  filter(!is.na(PAV_role) & !is.na(PAV_quant.mg)) %>%
  mutate(quant.mg = PAV_quant.mg,
         log_quant.mg = log(quant.mg + 1),
         log10_quant.mg = log10(quant.mg + 1),
         role = PAV_role,
         quant.plant = quant.mg * shoot_mass.g,
         log_quant.plant = log(quant.plant + 1),
         log10_quant.plant = log10(quant.plant + 1))

RPVdat <- dat5 %>%
  filter(!is.na(RPV_role) & !is.na(RPV_quant.mg)) %>%
  mutate(quant.mg = RPV_quant.mg,
         log_quant.mg = log(quant.mg + 1),
         log10_quant.mg = log10(quant.mg + 1),
         role = RPV_role,
         quant.plant = quant.mg * shoot_mass.g,
         log_quant.plant = log(quant.plant + 1),
         log10_quant.plant = log10(quant.plant + 1))

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
AIC(PAVURmod3, PAVURmod4) # linear (by 6)
AIC(RPVURmod3, RPVURmod4) # polynomial

# export
mod_sum(PAVURmod4, "pav_established_per_plant_model")
mod_sum(RPVURmod3, "rpv_established_per_plant_model")

# effect sizes

# function to exponentiate effects
eff_exp_fun <- function(mod, coef1, coef2){
  
  titer1 <- exp(fixef(mod)$cond[coef1])
  titer2 <- exp(fixef(mod)$cond[coef1] + fixef(mod)$cond[coef2])
  return(100 * (titer2 - titer1) / titer1)
  
}

# PAV establishment lower for single than resident
# titer
eff_exp_fun(PAVURmod2, 1, 4)

# virions/plant
eff_exp_fun(PAVURmod4, 1, 4)

# N reduces RPV's squared term
# titer
100* fixef(RPVURmod1)$cond[11]/fixef(RPVURmod1)$cond[6]

# PAV establishment declines with time
# virions/plant
eff_exp_fun(PAVURmod4, 1, 5)


#### visualize ####

# dodge points
dodge_size = 2

### invasion figures ###

PAVISimFig <- PAVIdat %>%
  select(role, dpiI, Nutrient, highN, highP) %>%
  unique() %>%
  mutate(set = NA) %>%
  mutate(log_quant.mg = predict(PAVImod, newdata = ., re.form = NA),
         log_quant.se = predict(PAVImod, newdata = ., se.fit = T, re.form = NA)$se.fit)

PAVIfig <- ggplot(PAVIdat, aes(x = dpiI, y = log_quant.mg, color = Nutrient, fill = Nutrient, linetype = Nutrient)) +
  geom_ribbon(data = PAVISimFig, aes(ymin = log_quant.mg - log_quant.se, ymax = log_quant.mg + log_quant.se), 
              color = NA, alpha = 0.1) +
  geom_line(data = PAVISimFig) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, alpha = 0.8, 
               position = position_dodge(dodge_size), linetype = "solid") +
  stat_summary(geom = "point", size = 2, fun = "mean", position = position_dodge(dodge_size),
               shape = 21, color = "black") +
  scale_color_viridis_d(direction = -1) +
  scale_fill_viridis_d(direction = -1) +
  scale_y_continuous(limits  = c(-1, 9)) +
  labs(y = "BYDV-PAV titer (ln[x + 1])", title = "Invading virus") +
  fig_theme +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

RPVISimFig <- RPVIdat %>%
  select(role, dpiI, Nutrient, highN, highP) %>%
  unique() %>%
  mutate(set = NA) %>%
  mutate(log_quant.mg = predict(RPVImod, newdata = ., re.form = NA),
         log_quant.se = predict(RPVImod, newdata = ., se.fit = T, re.form = NA)$se.fit)

RPVIfig <- ggplot(RPVIdat, aes(x = dpiI, y = log_quant.mg, color = Nutrient, fill = Nutrient, linetype = Nutrient)) +
  geom_ribbon(data = RPVISimFig, aes(ymin = log_quant.mg - log_quant.se, ymax = log_quant.mg + log_quant.se), 
              color = NA, alpha = 0.1, show.legend = F) +
  geom_line(data = RPVISimFig) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, alpha = 0.8, 
               position = position_dodge(dodge_size), linetype = "solid") +
  stat_summary(geom = "point", size = 2, fun = "mean", position = position_dodge(dodge_size),
               shape = 21, color = "black") +
  scale_color_viridis_d(direction = -1) +
  scale_fill_viridis_d(direction = -1) +
  scale_y_continuous(limits  = c(12, 17)) +
  labs(x = "Days post invader inoculation", y = "CYDV-RPV titer (ln[x + 1])") +
  fig_theme +
  theme(legend.position = "none")

### resident figures ###

PAVRSimFig <- PAVURdat %>%
  filter(role == "resident") %>%
  select(role, dpiR, dpiUR, Nutrient, highN, highP) %>%
  unique() %>%
  mutate(set = NA) %>%
  mutate(log_quant.mg = predict(PAVURmod, newdata = ., re.form = NA),
         log_quant.se = predict(PAVURmod, newdata = ., se.fit = T, re.form = NA)$se.fit)

PAVRfig <- PAVdat %>%
  filter(role == "resident") %>%
  ggplot(aes(dpiR, log_quant.mg, color = Nutrient, fill = Nutrient, linetype = Nutrient)) +
  geom_ribbon(data = PAVRSimFig, aes(ymin = log_quant.mg - log_quant.se, ymax = log_quant.mg + log_quant.se), 
              color = NA, alpha= 0.1, show.legend = F) +
  geom_line(data = PAVRSimFig,) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, alpha = 0.8, 
               position = position_dodge(dodge_size), linetype = "solid", show.legend = F) +
  stat_summary(geom = "point", size = 2, fun = "mean", position = position_dodge(dodge_size),
               shape = 21, color = "black") +
  scale_color_viridis_d(direction = -1, name = "Nutrient\nsupply") +
  scale_fill_viridis_d(direction = -1, name = "Nutrient\nsupply") +
  scale_linetype(name = "Nutrient\nsupply") +
  scale_y_continuous(limits  = c(-1, 9)) +
  labs(title = "Resident virus") +
  fig_theme +
  theme(axis.title = element_blank(),
        legend.key.width = unit(1, "cm"),
        plot.title = element_text(hjust = 0.5))

RPVRSimFig <- tibble(dpiR = seq(min(RPVURdat$dpiR), max(RPVURdat$dpiR), length.out = 50),
                     dpiUR = seq(min(RPVURdat$dpiUR), max(RPVURdat$dpiUR), length.out = 50),
                     role = "resident") %>%
  expand_grid(tibble(highN = c(0, 1, 0, 1),
                     highP = c(0, 0, 1, 1),
                     Nutrient = c("low", "+N", "+P", "+N+P") %>%
                       fct_relevel("low", "+N", "+P"))) %>%
  mutate(set = NA) %>%
  mutate(log_quant.mg = predict(RPVURmod, newdata = ., re.form = NA),
         log_quant.se = predict(RPVURmod, newdata = ., se.fit = T, re.form = NA)$se.fit)

RPVRfig <- RPVdat %>%
  filter(role == "resident") %>%
  ggplot(aes(dpiR, log_quant.mg, color = Nutrient, fill = Nutrient, linetype = Nutrient)) +
  geom_ribbon(data = RPVRSimFig, aes(ymin = log_quant.mg - log_quant.se, ymax = log_quant.mg + log_quant.se), 
              color = NA, alpha = 0.1, show.legend = F) +
  geom_line(data = RPVRSimFig) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, alpha = 0.8, 
               position = position_dodge(dodge_size), linetype = "solid") +
  stat_summary(geom = "point", size = 2, fun = "mean", position = position_dodge(dodge_size),
               shape = 21, color = "black") +
  scale_color_viridis_d(direction = -1) +
  scale_fill_viridis_d(direction = -1) +
  scale_y_continuous(limits  = c(12, 17)) +
  labs(x = "Days post resident inoculation") +
  fig_theme +
  theme(legend.position = "none",
        axis.title.y = element_blank())

# combine figures
comb_fig <- PAVIfig + theme(axis.text.x = element_blank()) + 
  PAVRfig + theme(axis.text.x = element_blank(),
                  axis.text.y = element_blank()) + 
  RPVIfig + RPVRfig + theme(axis.text.y = element_blank()) +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 8, face = "bold"))

ggsave("output/PAV_RPV_virus_titer_figure.pdf", comb_fig,
       width = 6, height = 5, units = "in")

# outlier
filter(PAVdat, role == "only" & nutrient == "L" & time == 3) %>%
  select(tube_label, quant.mg, expt_notes, extraction_notes)
# looked at raw data from qPCR_data_processing
# one technical replicate was > 1000, the other two were not, their average was not
# averages below 1000 (or minimum standard) were converted to zero
