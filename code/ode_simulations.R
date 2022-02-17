## Goal: simulate mutual invasions


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(deSolve)
library(cowplot)
library(lemon)
library(janitor)
library(ggpattern)

# load model and settings
source("code/model_settings.R")

# settings for short-term simulations
plant_days <- 11
res_days <- 12
inv_days <- 19


#### long-term plant ####

# simulation
long_term_plant <- virus2_model_sim(params_def2, "PAV", V0_b = 0, V0_c = 0,
                                    plant_time = 250, res_time = 12, 
                                    inv_time = 7500-250-12) %>%
  virus2_model_format(params_def2)

# figure
pdf("output/long_term_plant_simulation_figure.pdf", width = 6.5, height = 3.75)
plant_fig_fun(long_term_plant, params_def2, -8e-4, -5e-4)
dev.off()


#### short-term plant ####

# simulation
short_term_plant <- virus2_model_sim(params_def2, "PAV", V0_b = 0, V0_c = 0, 
                                     plant_time = plant_days, res_time = res_days, 
                                     inv_time = inv_days) %>%
  virus2_model_format(params_def2)

# figure
pdf("output/short_term_plant_simulation_figure.pdf", width = 6.5, height = 3.75)
plant_fig_fun(short_term_plant, params_def2, -7e-4, -4e-4)
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
  filter(time > (plant_days + res_days) & str_starts(variable2, "PAV") == T) %>%
  mutate(scenario = "(A) Invader")

rpv_inv_sim2 <- rpv_inv_sim %>%
  filter(time > (plant_days + res_days) & str_starts(variable2, "RPV") == T) %>%
  mutate(scenario = "(A) Invader")

pav_second_sim2 <- pav_second_sim %>%
  filter(time > (plant_days + res_days) & str_starts(variable2, "PAV") == T) %>%
  mutate(scenario = "(B) Invader alone")

rpv_second_sim2 <- rpv_second_sim %>%
  filter(time > (plant_days + res_days) & str_starts(variable2, "RPV") == T) %>%
  mutate(scenario = "(B) Invader alone")

pav_res_sim2 <- rpv_inv_sim %>%
  filter(time > plant_days & str_starts(variable2, "PAV") == T) %>%
  mutate(scenario = "(C) Resident")

rpv_res_sim2 <- pav_inv_sim %>%
  filter(time > plant_days & str_starts(variable2, "RPV") == T) %>%
  mutate(scenario = "(C) Resident")

pav_first_sim2 <- pav_first_sim %>%
  filter(time > plant_days & str_starts(variable2, "PAV") == T) %>%
  mutate(scenario = "(D) Resident alone")

rpv_first_sim2 <- rpv_first_sim %>%
  filter(time > plant_days & str_starts(variable2, "RPV") == T) %>%
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
  mutate(virus_abund = log10(value + 1),
         virus = str_sub(variable2, 1, 3),
         abund_type = str_sub(variable2, 5, 7),
         fac_lab_y = case_when(virus == "PAV" & abund_type == "con" ~ "Log[10]~BYDV-PAV~(g^-1)",
                               virus == "RPV" & abund_type == "con" ~ "Log[10]~CYDV-RPV~(g^-1)",
                               virus == "PAV" & abund_type == "pop" ~ "Log[10]~BYDV-PAV~(plant^-1)",
                               virus == "RPV" & abund_type == "pop" ~ "Log[10]~CYDV-RPV~(plant^-1)"),
         fac_lab_x = str_replace_all(scenario, " ", "~"),
         fac_lab_x = paste0("bold(", fac_lab_x, ")"))


#### virus figure ####

# tried to align strip x labels left, 
# but resident alone only shows full text with hjust = 0.5

pdf("output/invasion_simulation_concentration_figure.pdf", width = 6.5, height = 3.75)
ggplot(filter(virus_sim, abund_type == "con"), 
       aes(x = time, y = virus_abund, color = nutrient, linetype = nutrient)) +
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

pdf("output/invasion_simulation_population_figure.pdf", width = 6.5, height = 3.75)
ggplot(filter(virus_sim, abund_type == "pop"), 
       aes(x = time, y = virus_abund, color = nutrient, linetype = nutrient)) +
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
plant_fig_fun(pav_first_sim, params_def2, -7e-4, -4e-4)
dev.off()

pdf("output/rpv_only_plant_simulation_figure.pdf", width = 6.5, height = 3.75)
plant_fig_fun(rpv_first_sim, params_def2, -7e-4, -4e-4)
dev.off()

pdf("output/pav_inv_plant_simulation_figure.pdf", width = 6.5, height = 3.75)
plant_fig_fun(pav_inv_sim, params_def2, -7e-4, -4e-4)
dev.off()

pdf("output/rpv_inv_plant_simulation_figure.pdf", width = 6.5, height = 3.75)
plant_fig_fun(rpv_inv_sim, params_def2, -7e-4, -4e-4)
dev.off()


#### virus growth rates ####

pav_inv_gr <- pav_inv_sim %>%
  virus2_growth_rate("RPV", plant_days, res_days) %>%
  mutate(virus = str_sub(variable2, 1, 3),
         DPP = if_else(virus == "PAV", plant_days + res_days, plant_days))

rpv_inv_gr <- rpv_inv_sim %>%
  virus2_growth_rate("PAV", plant_days, res_days) %>%
  mutate(virus = str_sub(variable2, 1, 3),
         DPP = if_else(virus == "RPV", plant_days + res_days, plant_days))

pav_second_gr <- pav_second_sim %>%
  virus2_growth_rate("RPV", plant_days, res_days) %>%
  mutate(virus = str_sub(variable2, 1, 3),
         DPP = plant_days + res_days) %>%
  filter(virus == "PAV")

rpv_second_gr <- rpv_second_sim %>%
  virus2_growth_rate("PAV", plant_days, res_days) %>%
  mutate(virus = str_sub(variable2, 1, 3),
         DPP = plant_days + res_days) %>%
  filter(virus == "RPV")

pav_first_gr <- pav_first_sim %>%
  virus2_growth_rate("PAV", plant_days, res_days) %>%
  mutate(virus = str_sub(variable2, 1, 3),
         DPP = plant_days) %>%
  filter(virus == "PAV")

rpv_first_gr <- rpv_first_sim %>%
  virus2_growth_rate("RPV", plant_days, res_days) %>%
  mutate(virus = str_sub(variable2, 1, 3),
         DPP = plant_days) %>%
  filter(virus == "RPV")

# combine
inv_gr <- pav_inv_gr %>%
  full_join(rpv_inv_gr) %>%
  mutate(role = if_else(DPP == 11, "resident", "invader") %>%
           fct_relevel("resident"))

alone_gr <- pav_first_gr %>%
  full_join(rpv_first_gr) %>%
  full_join(pav_second_gr) %>%
  full_join(rpv_second_gr) %>%
  mutate(role = if_else(DPP == 11, "resident", "invader") %>%
           fct_relevel("resident"))


#### growth rate figures ####

pdf("output/invasion_growth_rate_figure.pdf", width = 4.5, height = 2.5)
ggplot(inv_gr, aes(x = DPP, y = growth)) +
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  geom_col_pattern(aes(pattern = role, fill = nutrient),
                   position = "dodge",
                   color = "black", 
                   pattern_fill = "black",
                   pattern_color = "black",
                   pattern_density = 0.1,
                   pattern_key_scale_factor = 0.6) +
  facet_wrap(~ virus) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf, size = 1) +
  scale_pattern_manual(values = c(resident = "none", invader = "stripe"),
                       name = "Virus\nrole") +
  scale_x_continuous(breaks = c(plant_days, plant_days + res_days)) +
  scale_fill_viridis_d(end = 0.7, name = "Nutrient\nsupply") +
  labs(y = "Growth rate when rare", x = "Days post planting") +
  fig_theme + 
  theme(panel.spacing.x = unit(1,"line"),
        axis.line.y.left = element_blank()) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none")))
dev.off()

pdf("output/alone_growth_rate_figure.pdf", width = 4.5, height = 2.5)
ggplot(alone_gr, aes(x = DPP, y = growth)) +
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  geom_col_pattern(aes(pattern = role, fill = nutrient),
                   position = "dodge",
                   color = "black", 
                   pattern_fill = "black",
                   pattern_color = "black",
                   pattern_density = 0.1,
                   pattern_key_scale_factor = 0.6) +
  facet_wrap(~ virus) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf, size = 1) +
  scale_pattern_manual(values = c(resident = "none", invader = "none"),
                       guide = "none") +
  scale_x_continuous(breaks = c(plant_days, plant_days + res_days)) +
  scale_fill_viridis_d(end = 0.7, name = "Nutrient\nsupply") +
  labs(y = "Growth rate when rare", x = "Days post planting") +
  fig_theme +
  theme(panel.spacing.x = unit(1,"line"),
        axis.line.y.left = element_blank()) +
  guides(fill = guide_legend(override.aes = list(pattern = "none")))
dev.off()