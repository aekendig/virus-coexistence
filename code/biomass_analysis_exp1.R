## Goal: analyze plant growth data from exp1


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)

# import data
dat <- read_csv("data/concentration_analysis_all_data_exp1.csv")
# Google Drive/Within Host Coexistence/nutrients-plant-viruses/output/concentration_analysis_all_data.csv
sdat <- read_csv("data/sample_exp_molc_data_exp1.csv")
# Google Drive/Within Host Coexistence/nutrients-plant-viruses/data/sample_exp_molc_data.csv


#### edit data ####

# remove failed infections and infected healthy plants
# add sdat
qdat2 <- qdat %>%
  filter((quant_zero == 0 & inoc %in% c("PAV", "coinfection", "RPV")) | inoc == "healthy") %>%
  left_join(sdat %>%
              filter(material == "shoot") %>%
              select(round, time, inoc, nutrient, replicate, RTPCR_PAV, RTPCR_RPV)) %>%
  mutate(conc = quant_ul * 50 / mass_ext_mg,
         tot_mass_g = shoot_mass_g + root_mass_g)

# accidental RPV inoculations
(acc.r <- qdat2 %>%
    filter(inoc %in% c("PAV", "healthy") & ((target == "RPV" & quant_zero == 0) | RTPCR_RPV == 1)) %>%
    select(sample, round, time, inoc, nutrient, replicate, quant_zero, conc, RTPCR_RPV))

# accidental PAV inoculations  
(acc.p <- qdat2 %>%
    filter(inoc %in% c("RPV", "healthy") & ((target == "PAV" & quant_zero == 0) | RTPCR_PAV == 1)) %>%
    select(sample, round, time, inoc, nutrient, replicate, quant_zero, conc, RTPCR_PAV)) 

# remove accidental inoculations
qdat3 <- qdat2 %>%
  filter(!(sample %in% c(acc.r$sample, acc.p$sample)))

# missing data
sum(is.na(qdat3$shoot_mass_g))


#### visualize ####

ggplot(qdat3, aes(x = time, y = shoot_mass_g, color = nutrient)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_se") +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~inoc)
# a lot of healthy data missing


ggplot(qdat3, aes(x = time, y = root_mass_g, color = nutrient)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_se") +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~inoc)


ggplot(qdat3, aes(x = time, y = tot_mass_g, color = nutrient)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_se") +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~inoc)
