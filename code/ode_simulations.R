## Goal: simulate mutual invasions


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(deSolve)
library(cowplot)
library(lemon)

# load model and settings
source("code/model_settings.R")


#### long-term plant ####

# simulation
long_term_plant <- virus2_model_sim(params_def2, "PAV", V0_b = 0, V0_c = 0,
                                    plant_time = 250, res_time = 12, 
                                    inv_time = 1000-250-12) %>%
  virus2_model_format(params_def2)

# figure
pdf("output/long_term_plant_simulation_figure.pdf", width = 6.5, height = 3.75)
plant_fig_fun(long_term_plant, params_def2, -3e-3, -1e-3)
dev.off()


#### short-term plant ####

# settings for short-term simulations
plant_days <- 11
res_days <- 12
inv_days <- 1000-plant_days-res_days

# simulation
short_term_plant <- virus2_model_sim(params_def2, "PAV", V0_b = 0, V0_c = 0, 
                                     plant_time = plant_days, res_time = res_days, 
                                     inv_time = inv_days) %>%
  virus2_model_format(params_def2)

# figure
pdf("output/short_term_plant_simulation_figure.pdf", width = 6.5, height = 3.75)
plant_fig_fun(short_term_plant, params_def2, -2e-3, -5e-4)
dev.off()


#### single virus simulations ####

# PAV with plant
pav_first_sim <- virus2_model_sim(params_def2, "PAV", V0_b = V0, V0_c = 0,
                                  plant_time = plant_days, res_time = res_days, 
                                  inv_time = inv_days) %>%
  virus2_model_format(params_def2)

pav_second_sim <- virus2_model_sim(params_def2, "RPV", V0_b = V0, V0_c = 0,
                                   plant_time = plant_days, res_time = res_days, 
                                   inv_time = inv_days) %>%
  virus2_model_format(params_def2)

# RPV with plant
rpv_first_sim <- virus2_model_sim(params_def2, "RPV", V0_b = 0, V0_c = V0,
                                  plant_time = plant_days, res_time = res_days, 
                                  inv_time = inv_days) %>%
  virus2_model_format(params_def2)

rpv_second_sim <- virus2_model_sim(params_def2, "PAV", V0_b = 0, V0_c = V0,
                                   plant_time = plant_days, res_time = res_days, 
                                   inv_time = inv_days) %>%
  virus2_model_format(params_def2)


#### virus invasion simulations ####

# PAV invades
pav_inv_sim <- virus2_model_sim(params_def2, "RPV", V0_b = V0, V0_c = V0,
                                plant_time = plant_days, res_time = res_days, 
                                inv_time = inv_days) %>%
  virus2_model_format(params_def2)

# RPV invades
rpv_inv_sim <- virus2_model_sim(params_def2, "PAV", V0_b = V0, V0_c = V0,
                                plant_time = plant_days, res_time = res_days, 
                                inv_time = inv_days) %>%
  virus2_model_format(params_def2)


#### edit simulation data ####

# edit simulations
pav_inv_sim2 <- pav_inv_sim %>%
  filter(time > (plant_days + res_days) & variable2 == "PAV") %>%
  mutate(scenario = "(A) Invader")

rpv_inv_sim2 <- rpv_inv_sim %>%
  filter(time > (plant_days + res_days) & variable2 == "RPV") %>%
  mutate(scenario = "(A) Invader")

pav_second_sim2 <- pav_second_sim %>%
  filter(time > (plant_days + res_days) & variable2 == "PAV") %>%
  mutate(scenario = "(B) Invader alone")

rpv_second_sim2 <- rpv_second_sim %>%
  filter(time > (plant_days + res_days) & variable2 == "RPV") %>%
  mutate(scenario = "(B) Invader alone")

pav_res_sim2 <- rpv_inv_sim %>%
  filter(time > plant_days & variable2 == "PAV") %>%
  mutate(scenario = "(C) Resident")

rpv_res_sim2 <- pav_inv_sim %>%
  filter(time > plant_days & variable2 == "RPV") %>%
  mutate(scenario = "(C) Resident")

pav_first_sim2 <- pav_first_sim %>%
  filter(time > plant_days & variable2 == "PAV") %>%
  mutate(scenario = "(D) Resident alone")

rpv_first_sim2 <- rpv_first_sim %>%
  filter(time > plant_days & variable2 == "RPV") %>%
  mutate(scenario = "(D) Resident alone")

# combine
virus_sim <- pav_inv_sim2 %>%
  full_join(rpv_inv_sim2) %>%
  full_join(pav_second_sim2) %>%
  full_join(rpv_second_sim2) %>%
  full_join(pav_res_sim2) %>%
  full_join(rpv_res_sim2) %>%
  full_join(pav_first_sim2) %>%
  full_join(rpv_first_sim2) %>%
  mutate(virus_conc = log10(value + 1),
         virus = variable2,
         fac_lab_y = case_when(virus == "PAV" ~ "BYDV-PAV~titer~(log[10])",
                               virus == "RPV" ~ "CYDV-RPV~titer~(log[10])"),
         fac_lab_x = str_replace_all(scenario, " ", "~"),
         fac_lab_x = paste0("bold(", fac_lab_x, ")"))


#### virus figure ####

# tried to align strip x labels left, 
# but resident alone only shows full text with hjust = 0.5

pdf("output/invasion_simulation_figure.pdf", width = 6.5, height = 3.75)
ggplot(virus_sim, 
       aes(x = time, y = virus_conc, color = nutrient, linetype = nutrient)) +
  geom_line(aes(size = lim_nut_H)) +
  facet_rep_grid(rows = vars(fac_lab_y),
                 cols = vars(fac_lab_x),
                 labeller = label_parsed,
                 switch = "y",
                 scales = "free_y") +
  scale_color_viridis_d(end = 0.7, name = "Nutrient\nsupply") +
  scale_linetype(name = "Nutrient\nsupply") +
  scale_size_manual(values = c(1.2, 0.6), name = "Limiting\nnutrient",
                    labels = c("N", "P")) +
  labs(x = "Time (days)") +
  fig_theme +
  theme(axis.title.y = element_blank()) +
  guides(color = guide_legend(order = 1, override.aes = list(size = 1.2)), linetype = guide_legend(order = 1))
dev.off()


#### plant figures with viruses ####

pdf("output/pav_only_plant_simulation_figure.pdf", width = 6.5, height = 3.75)
plant_fig_fun(pav_first_sim, params_def2, -2e-3, -5e-4)
dev.off()

pdf("output/rpv_only_plant_simulation_figure.pdf", width = 6.5, height = 3.75)
plant_fig_fun(rpv_first_sim, params_def2, -2e-3, -5e-4)
dev.off()

pdf("output/pav_inv_plant_simulation_figure.pdf", width = 6.5, height = 3.75)
plant_fig_fun(pav_inv_sim, params_def2, -2e-3, -5e-4)
dev.off()

pdf("output/rpv_inv_plant_simulation_figure.pdf", width = 6.5, height = 3.75)
plant_fig_fun(rpv_inv_sim, params_def2, -2e-3, -5e-4)
dev.off()