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

# figure labels
Qn_lab <- tibble(time = 0, 
                 Q_n = as.numeric(params_def2["Qmin_n"]), 
                 label = "Q['min,N']")
Qp_lab <- tibble(time = 0, 
                 Q_p = as.numeric(params_def2["Qmin_p"]), 
                 label = "Q['min,P']")


#### long-term plant ####

# viruses
V0_b <- V0_c <- 0

# simulation
long_term_plant <- virus2_model_sim(params_def2, "PAV", 
                                    plant_time = 250, res_time = 12, 
                                    inv_time = 1000-250-12) %>%
  virus2_model_format(params_def2)

#### start here: check that above works as expected by making figure ###

# figure
pdf("output/long_term_plant_simulation_figure.pdf", width = 6.5, height = 3.75)
plant_fig_fun(low_plant_long_sim, n_plant_long_sim, p_plant_long_sim, np_plant_long_sim, -3e-3, -1e-3)
dev.off()


#### short-term plant ####

# time intervals
plant_days <- 11
res_days <- 12
inv_days <- 100-plant_days-res_days

# viruses
V0_b <- V0_c <- 0

# simulations
low_plant_short_sim <- sim_fun("low", "low", "low", "PAV", plant_days, res_days, inv_days)
n_plant_short_sim <- sim_fun("high", "low", "N", "PAV", plant_days, res_days, inv_days)
p_plant_short_sim <- sim_fun("low", "high", "P", "PAV", plant_days, res_days, inv_days)
np_plant_short_sim <- sim_fun("high", "high", "N+P", "PAV", plant_days, res_days, inv_days)

# figure
pdf("output/short_term_plant_simulation_figure.pdf", width = 6.5, height = 3.75)
plant_fig_fun(low_plant_short_sim, n_plant_short_sim, p_plant_short_sim, np_plant_short_sim, -2e-3, -5e-4)
dev.off()


#### single virus simulations ####

# time intervals
plant_days <- 11
res_days <- 12
inv_days <- 100-plant_days-res_days

# PAV with plant
V0_b <- V0_init
V0_c <- 0

low_pav_first_sim <- sim_fun("low", "low", "low", "PAV", plant_days, res_days, inv_days)
n_pav_first_sim <- sim_fun("high", "low", "N", "PAV", plant_days, res_days, inv_days)
p_pav_first_sim <- sim_fun("low", "high", "P", "PAV", plant_days, res_days, inv_days)
np_pav_first_sim <- sim_fun("high", "high", "N+P", "PAV", plant_days, res_days, inv_days)

low_pav_second_sim <- sim_fun("low", "low", "low", "RPV", plant_days, res_days, inv_days)
n_pav_second_sim <- sim_fun("high", "low", "N", "RPV", plant_days, res_days, inv_days)
p_pav_second_sim <- sim_fun("low", "high", "P", "RPV", plant_days, res_days, inv_days)
np_pav_second_sim <- sim_fun("high", "high", "N+P", "RPV", plant_days, res_days, inv_days)

# RPV with plant
V0_b <- 0
V0_c <- V0_init

low_rpv_first_sim <- sim_fun("low", "low", "low", "RPV", plant_days, res_days, inv_days)
n_rpv_first_sim <- sim_fun("high", "low", "N", "RPV", plant_days, res_days, inv_days)
p_rpv_first_sim <- sim_fun("low", "high", "P", "RPV", plant_days, res_days, inv_days)
np_rpv_first_sim <- sim_fun("high", "high", "N+P", "RPV", plant_days, res_days, inv_days)

low_rpv_second_sim <- sim_fun("low", "low", "low", "PAV", plant_days, res_days, inv_days)
n_rpv_second_sim <- sim_fun("high", "low", "N", "PAV", plant_days, res_days, inv_days)
p_rpv_second_sim <- sim_fun("low", "high", "P", "PAV", plant_days, res_days, inv_days)
np_rpv_second_sim <- sim_fun("high", "high", "N+P", "PAV", plant_days, res_days, inv_days)


#### virus invasion simulations ####

# time intervals
plant_days <- 11
res_days <- 12
inv_days <- 100-plant_days-res_days

# viruses
V0_b <- V0_c <- V0_init

# simulations
low_pav_inv_sim <- sim_fun("low", "low", "low", "RPV", plant_days, res_days, inv_days)
n_pav_inv_sim <- sim_fun("high", "low", "N", "RPV", plant_days, res_days, inv_days)
p_pav_inv_sim <- sim_fun("low", "high", "P", "RPV", plant_days, res_days, inv_days)
np_pav_inv_sim <- sim_fun("high", "high", "N+P", "RPV", plant_days, res_days, inv_days)

low_rpv_inv_sim <- sim_fun("low", "low", "low", "PAV", plant_days, res_days, inv_days)
n_rpv_inv_sim <- sim_fun("high", "low", "N", "PAV", plant_days, res_days, inv_days)
p_rpv_inv_sim <- sim_fun("low", "high", "P", "PAV", plant_days, res_days, inv_days)
np_rpv_inv_sim <- sim_fun("high", "high", "N+P", "PAV", plant_days, res_days, inv_days)


#### edit simulation data ####

pav_first_sim <- sim_dat_fun(low_pav_first_sim, n_pav_first_sim, 
                             p_pav_first_sim, np_pav_first_sim) %>%
  filter(time > plant_days) %>%
  mutate(virus_conc = PAV_log10,
         virus = resident,
         scenario = "(D) Resident alone")

pav_second_sim <- sim_dat_fun(low_pav_second_sim, n_pav_second_sim, 
                              p_pav_second_sim, np_pav_second_sim) %>%
  filter(time > (plant_days + res_days)) %>%
  mutate(virus_conc = PAV_log10,
         virus = invader,
         scenario = "(B) Invader alone")

rpv_first_sim <- sim_dat_fun(low_rpv_first_sim, n_rpv_first_sim, 
                             p_rpv_first_sim, np_rpv_first_sim) %>%
  filter(time > plant_days) %>%
  mutate(virus_conc = RPV_log10,
         virus = resident,
         scenario = "(D) Resident alone")

rpv_second_sim <- sim_dat_fun(low_rpv_second_sim, n_rpv_second_sim, 
                              p_rpv_second_sim, np_rpv_second_sim) %>%
  filter(time > (plant_days + res_days)) %>%
  mutate(virus_conc = RPV_log10,
         virus = invader,
         scenario = "(B) Invader alone")

pav_res_sim <- sim_dat_fun(low_rpv_inv_sim, n_rpv_inv_sim, 
                           p_rpv_inv_sim, np_rpv_inv_sim) %>%
  filter(time > plant_days) %>%
  mutate(virus_conc = PAV_log10,
         virus = resident,
         scenario = "(C) Resident")

rpv_inv_sim <- pav_res_sim %>%
  filter(time > (plant_days + res_days)) %>%
  mutate(virus_conc = RPV_log10,
         virus = invader,
         scenario = "(A) Invader")

rpv_res_sim <- sim_dat_fun(low_pav_inv_sim, n_pav_inv_sim, 
                           p_pav_inv_sim, np_pav_inv_sim) %>%
  filter(time > plant_days) %>%
  mutate(virus_conc = RPV_log10,
         virus = resident,
         scenario = "(C) Resident")

pav_inv_sim <- rpv_res_sim %>%
  filter(time > (plant_days + res_days)) %>%
  mutate(virus_conc = PAV_log10,
         virus = invader,
         scenario = "(A) Invader")

# combine
virus_sim <- pav_inv_sim %>%
  full_join(rpv_inv_sim) %>%
  full_join(pav_second_sim) %>%
  full_join(rpv_second_sim) %>%
  full_join(pav_res_sim) %>%
  full_join(rpv_res_sim) %>%
  full_join(pav_first_sim) %>%
  full_join(rpv_first_sim) %>%
  mutate(fac_lab_y = case_when(virus == "PAV" ~ "BYDV-PAV~titer~(log[10])",
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
plant_fig_fun(low_pav_first_sim, n_pav_first_sim, 
              p_pav_first_sim, np_pav_first_sim, -2e-3, -5e-4)
dev.off()

pdf("output/rpv_only_plant_simulation_figure.pdf", width = 6.5, height = 3.75)
plant_fig_fun(low_rpv_first_sim, n_rpv_first_sim, 
              p_rpv_first_sim, np_rpv_first_sim, -2e-3, -5e-4)
dev.off()

pdf("output/pav_inv_plant_simulation_figure.pdf", width = 6.5, height = 3.75)
plant_fig_fun(low_pav_inv_sim, n_pav_inv_sim, 
              p_pav_inv_sim, np_pav_inv_sim, -2e-3, -5e-4)
dev.off()

pdf("output/rpv_inv_plant_simulation_figure.pdf", width = 6.5, height = 3.75)
plant_fig_fun(low_rpv_inv_sim, n_rpv_inv_sim, 
              p_rpv_inv_sim, np_rpv_inv_sim, -2e-3, -5e-4)
dev.off()