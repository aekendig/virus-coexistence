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
Qn_lab <- tibble(time = 0, Q_n = Qmin_n, label = "Q['min,N']")
Qp_lab <- tibble(time = 0, Q_p = Qmin_p, label = "Q['min,P']")


#### figure functions ####

#### start here: formatting function below in model_settings ####

# simulation data formatting
sim_dat_fun <- function(low_sim, n_sim, p_sim, np_sim){
  
  dat_out <- low_sim %>%
    full_join(n_sim) %>%
    full_join(p_sim) %>%
    full_join(np_sim) %>%
    mutate(Qlim_n = Qmin_n / Q_n,
           Qlim_p = Qmin_p / Q_p) %>%
    rowwise() %>%
    mutate(Qlim = max(Qlim_n, Qlim_p)) %>%
    ungroup() %>%
    mutate(nutrient = fct_recode(nutrient, "+N" = "N",
                                 "+P" = "P",
                                 "+N+P" = "N+P") %>%
             fct_relevel("low", "+N", "+P"),
           lim_nut_H = case_when(Qlim == Qlim_n ~ "Q[N]",
                                 Qlim == Qlim_p ~ "Q[P]"),
           PAV_log10 = log10(V_b + 1),
           RPV_log10 = log10(V_c + 1)) %>%
    rename(PAV = V_b, RPV = V_c)
  
  return(dat_out)
  
}

# plant_figure function
plant_fig_fun <- function(low_sim, n_sim, p_sim, np_sim, q_adj_n, q_adj_p){
  
  # combine
  plant_dat <- sim_dat_fun(low_sim, n_sim, p_sim, np_sim)
  
  # panels
  plant_H_fig <- ggplot(plant_dat, aes(time, H, linetype = nutrient, color = nutrient)) +
    geom_line(aes(size = lim_nut_H)) +
    scale_color_viridis_d(end = 0.7, name = "Nutrient\nsupply") +
    scale_linetype(name = "Nutrient\nsupply") +
    scale_size_manual(values = c(1.2, 0.6), guide= "none") +
    labs(x = "Time (days)", y = "Plant biomass (g)", title = "(A)") +
    fig_theme +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 1.2)))
  
  plant_lim_fig <- ggplot(plant_dat, aes(time, Qlim, linetype = nutrient, color = nutrient)) +
    geom_line(aes(size = lim_nut_H)) +
    scale_color_viridis_d(end = 0.7, guide = "none") +
    scale_linetype(guide = "none") +
    scale_size_manual(values = c(1.2, 0.6), name = "Limiting\nnutrient",
                      labels = c("N", "P")) +
    labs(x = "Time (days)", title = "(B)", 
         y = expression(paste("Limiting nutrient ratio (", Q[min], "/Q)", sep = ""))) +
    fig_theme +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank())
  
  plant_Rn_fig <- ggplot(plant_dat, aes(time, R_n, linetype = nutrient, color = nutrient)) +
    geom_line(aes(size = lim_nut_H)) +
    scale_color_viridis_d(end = 0.7) +
    scale_size_manual(values = c(1.2, 0.6)) +
    labs(x = "Time (days)", y = "Environment N (g)", title = "(C)") +
    fig_theme +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_blank())
  
  plant_Rp_fig <- ggplot(plant_dat, aes(time, R_p, linetype = nutrient, color = nutrient)) +
    geom_line(aes(size = lim_nut_H)) +
    scale_color_viridis_d(end = 0.7) +
    scale_size_manual(values = c(1.2, 0.6)) +
    labs(x = "", y = "Environment P (g)", title = "(D)") +
    fig_theme +
    theme(legend.position = "none")
  
  plant_Qn_fig <- ggplot(plant_dat, aes(time, Q_n)) +
    geom_hline(yintercept = Qmin_n, color = "black") +
    geom_text(data = Qn_lab, aes(label = label), parse = T, color = "black", fontface = "italic",
              size = 3, hjust = 0, vjust = 0, nudge_y = q_adj_n) +
    geom_line(aes(linetype = nutrient, color = nutrient, size = lim_nut_H)) +
    scale_color_viridis_d(end = 0.7) +
    scale_size_manual(values = c(1.2, 0.6)) +
    labs(x = "Time (days)", y = "Plant N concentration", title = "(E)") +
    fig_theme +
    theme(legend.position = "none")
  
  plant_Qp_fig <- ggplot(plant_dat, aes(time, Q_p)) +
    geom_hline(yintercept = Qmin_p, color = "black") +
    geom_text(data = Qp_lab, aes(label = label), parse = T, color = "black", fontface = "italic", 
              size = 3, hjust = 0, vjust = 0, nudge_y = q_adj_p) +
    geom_line(aes(linetype = nutrient, color = nutrient, size = lim_nut_H)) +
    scale_color_viridis_d(end = 0.7) +
    scale_size_manual(values = c(1.2, 0.6)) +
    labs(x = "", y = "Plant P concentration", title = "(F)") +
    fig_theme +
    theme(legend.position = "none")
  
  # extract legends
  nut_leg <- get_legend(plant_H_fig)
  lim_leg <- get_legend(plant_lim_fig)
  
  # combine figures
  plant_top_fig <- plot_grid(plant_H_fig + theme(legend.position = "none"), 
                             plant_lim_fig + theme(legend.position = "none"), 
                             plant_Rn_fig,
                             nut_leg,
                             nrow = 1,
                             rel_widths = c(1, 1, 1, 0.35))
  plant_bot_fig <- plot_grid(plant_Rp_fig, plant_Qn_fig, plant_Qp_fig, 
                             lim_leg,
                             nrow = 1,
                             rel_widths = c(1, 1, 1, 0.35))
  
  # output
  return(plot_grid(plant_top_fig, plant_bot_fig,
                   nrow = 2, rel_heights = c(0.85, 1)))
  
}


#### long-term plant ####

# time intervals
plant_days <- 250
res_days <- 12
inv_days <- 1000-plant_days-res_days

# viruses
V0_b <- V0_c <- 0

# simulations
low_plant_long_sim <- sim_fun("low", "low", "low", "PAV", plant_days, res_days, inv_days)
n_plant_long_sim <- sim_fun("high", "low", "N", "PAV", plant_days, res_days, inv_days)
p_plant_long_sim <- sim_fun("low", "high", "P", "PAV", plant_days, res_days, inv_days)
np_plant_long_sim <- sim_fun("high", "high", "N+P", "PAV", plant_days, res_days, inv_days)

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