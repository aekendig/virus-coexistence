## Goal: analyze plant growth data from exp2

## moved this code to qPCR_analysis.R - don't need this script anymore


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)
library(glmmTMB)

# import data
dat <- read_csv("intermediate-data/qPCR_analysis_script_data_cleaned.csv")

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

# add inoculation treatment
dat2 <- dat %>%
  mutate(inoculation = paste(invasion, first_inoculation, sep = "_"),
         log_shoot_mass.g = log(shoot_mass.g),
         harvest_day = dpiI - 5)


#### visualize ####

ggplot(dat2, aes(x = harvest_day, y = shoot_mass.g, color = nutrient)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_se") +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ inoculation)

ggplot(dat2, aes(x = harvest_day, y = log_shoot_mass.g, color = nutrient)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_se") +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ inoculation)


#### model ####

bio_mod <- glmmTMB(log_shoot_mass.g ~ highN * highP * inoculation + (1|set) + (1|time), data = dat2)
summary(bio_mod)
mod_sum(bio_mod, "biomass_model")
