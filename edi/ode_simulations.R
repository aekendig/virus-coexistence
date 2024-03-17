#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse) # version 2.0.0
library(deSolve) # version 1.38
library(cowplot) # version 1.1.1
library(lemon) # version 0.4.6
library(janitor) # version 2.2.0
library(ggpattern) # version 1.0.1
library(patchwork) # version 1.1.2
library(pals) # version 1.7

# load model and settings
source("code/model_settings.R")

# settings for short-term simulations
plant_days <- 11
res_days <- 12
inv_days <- 60-11-12


#### long-term plant ####

# simulation
long_term_plant <- virus2_model_sim(params_def2, "PAV", V0_b = 0, V0_c = 0,
                                    plant_time = 250, res_time = 12, 
                                    inv_time = 7500-250-12) %>%
  virus2_model_format(params_def2)

# Appendix S2: Figure S1
ggsave("output/long_term_plant_simulation_figure.pdf", 
       plant_fig_fun(long_term_plant, params_def2, -3e-3, -1e-3),
       width = 18, height = 10.5, units = "cm")


#### short-term plant ####

# simulation
short_term_plant <- virus2_model_sim(params_def2, "PAV", V0_b = 0, V0_c = 0, 
                                     plant_time = plant_days, res_time = res_days, 
                                     inv_time = inv_days) %>%
  virus2_model_format(params_def2)

# figure 1
ggsave("output/short_term_plant_simulation_figure.pdf",
       plant_fig_fun(short_term_plant, params_def2, -1e-5, -6e-4),
       width = 18, height = 10.5, units = "cm")

# save simulation
write_csv(short_term_plant, "output/short_term_plant_simulation.csv")


#### simulation function function ####

# change parameter value and simulate
param_fun <- function(params_in, 
                      param_foc1 = NA_character_, param_val1 = NA_real_, 
                      param_foc2 = NA_character_, param_val2 = NA_real_,
                      param_foc3 = NA_character_, param_val3 = NA_real_, 
                      param_foc4 = NA_character_, param_val4 = NA_real_,
                      first_virus, output_type = "growth",
                      V0_b = V0, V0_c = V0, inits_in = init_def2,
                      plant_time = plant_days, res_time = res_days, inv_time = inv_days){
  
  # update parameter value
  if(!is.na(param_foc1)){
    params_in[param_foc1] <- param_val1
  }
  
  if(!is.na(param_foc2)){
    params_in[param_foc2] <- param_val2
  }
  
  if(!is.na(param_foc3)){
    params_in[param_foc3] <- param_val3
  }
  
  if(!is.na(param_foc4)){
    params_in[param_foc4] <- param_val4
  }
  
  # run model
  mod_out <- virus2_model_sim(params = params_in, first_virus = first_virus, 
                              V0_b = V0_b, V0_c = V0_c,
                              plant_time = plant_time, res_time = res_time, inv_time = inv_time, 
                              inits = inits_in)
  
  # edit output
  mod_out2 <- virus2_model_format(mod_out, params_in, Qlim = F)
  
  # figure of time series
  print(ggplot(mod_out2, aes(x = time, y = value2, 
                             color = nutrient, linetype = nutrient)) +
          geom_line() +
          facet_wrap(~ variable2, scales = "free_y") +
          scale_color_viridis_d() +
          fig_theme +
          labs(title = paste0(param_foc1, " = ", param_val1, "; ", param_foc2, " = ", param_val2)))
  
  # outputs
  if(output_type == "growth"){
    
    # calculate virus growth rates
    mod_out3 <- virus2_growth_rate(mod_out2, first_virus,
                                   plant_time = plant_time, res_time = res_time)
    
  } else if(is.numeric(output_type) == T) {
    
    # extract specified time point
    mod_out3 <- mod_out2 %>%
      filter(time == output_type)
    
  } else {
    
    # calculate virus growth rates
    mod_out3a <- virus2_growth_rate(mod_out2, first_virus,
                                    plant_time = plant_time, res_time = res_time) %>%
      mutate(variable2 = str_replace(variable2, "conc", "inv_gr")) %>%
      rename(value2 = growth)
    
    # extract specified time point
    mod_out3b <- mod_out2 %>%
      filter(time == as.numeric(output_type))
    
    # combine
    mod_out3 <- full_join(mod_out3a, mod_out3b)
    
  }
  
  
  # return output
  return(mod_out3)
}

# test function
virus2_model_sim(params_def2, "RPV", V0_b = V0, V0_c = V0,
                 plant_time = plant_days, res_time = res_days, 
                 inv_time = inv_days) %>%
  virus2_model_format(params_def2) %>%
  filter(time == max(time) & variable2 %in% c("PAV_conc", "PAV_pop")) %>%
  mutate(virus_conc = log10(value + 1))

param_fun(params_def2, 
          param_foc1 = "q_nb", param_val1 = 1.1e-3, 
          param_foc2 = "q_pb", param_val2 = 7.4e-5, 
          first_virus = "RPV", output_type = 60) %>%
  filter(variable2 %in% c("PAV_conc", "PAV_pop"))
# output should match output of above

param_fun(params_def2, 
          param_foc1 = "a_n_lo", param_val1 = 1.1e-10, 
          first_virus = "RPV", output_type = 60) %>%
  filter(variable2 %in% c("PAV_conc", "PAV_pop"))
# output should match output of above unless N supply is low

# growth rate when rare
param_fun(params_def2, 
          param_foc1 = "q_nb", param_val1 = 1.1e-3, 
          param_foc2 = "q_pb", param_val2 = 7.4e-5, 
          first_virus = "RPV", output_type = "growth")

# growth rate and values
param_fun(params_def2, 
          param_foc1 = "q_nb", param_val1 = 1.1e-3, 
          param_foc2 = "q_pb", param_val2 = 7.4e-5, 
          first_virus = "RPV", output_type = "23")


#### days post planting ####

# day ranges
dpp_in <- tibble(plant_days = 1:60) %>% # days before resident
  mutate(res_days = 1, # days before invasion
         inv_days = 1) # days after invasion

# PAV alone
pdf("output/ode_simulations_pav_1st_dpp.pdf")
pav_1st_dpp <- dpp_in %>%
  mutate(sim_out = pmap(., function(plant_days, res_days, inv_days) param_fun(params_in = params_def2, plant_time = plant_days, res_time = res_days, inv_time = inv_days, first_virus = "PAV", output_type = "growth", V0_c = 0))) %>%
  unnest(cols = c(sim_out))
dev.off()

# RPV alone
pdf("output/ode_simulations_rpv_1st_dpp.pdf")
rpv_1st_dpp <- dpp_in %>%
  mutate(sim_out = pmap(., function(plant_days, res_days, inv_days) param_fun(params_in = params_def2, plant_time = plant_days, res_time = res_days, inv_time = inv_days, first_virus = "RPV", output_type = "growth", V0_b = 0))) %>%
  unnest(cols = c(sim_out))
dev.off()

# save datasets
write_csv(pav_1st_dpp, "output/ode_simulations_pav_1st_dpp.csv")
write_csv(rpv_1st_dpp, "output/ode_simulations_rpv_1st_dpp.csv")

# re-import datasets
pav_1st_dpp <- read_csv("output/ode_simulations_pav_1st_dpp.csv") %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P", "+N+P"))
rpv_1st_dpp <- read_csv("output/ode_simulations_rpv_1st_dpp.csv") %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P", "+N+P"))

# figures
pav_dpp_fig <- pav_1st_dpp %>%
  filter(variable2 == "PAV_conc") %>%
  ggplot(aes(x = plant_days, y = growth, color = nutrient, linetype = nutrient)) +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  geom_line(linewidth = 0.9) +
  scale_color_viridis_d(direction = -1, name = "Nutrient\nsupply") +
  scale_linetype(name = "Nutrient\nsupply") +
  labs(x = "Infection time (DPP)", y = "BYDV-PAV growth rate") +
  fig_theme

rpv_dpp_fig <- rpv_1st_dpp %>%
  filter(variable2 == "RPV_conc") %>%
  ggplot(aes(x = plant_days, y = growth, color = nutrient, linetype = nutrient)) +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  geom_line(linewidth = 0.9) +
  scale_color_viridis_d(direction = -1, name = "Nutrient\nsupply") +
  scale_linetype(name = "Nutrient\nsupply") +
  labs(x = "Infection time (DPP)", y = "CYDV-RPV growth rate") +
  fig_theme


#### virus nutrient use ####

# virus N:P
np_zb_ratio_vir <- as.numeric(params_def2["z_nb"]) / as.numeric(params_def2["z_pb"])
np_zc_ratio_vir <- as.numeric(params_def2["z_nc"]) / as.numeric(params_def2["z_pc"])

# values
# too high Z values make Q's go negative
z_pb_vals <- sort(c(as.numeric(params_def2["z_pb"]),
                    10^seq(-18, -7, length.out = 50)))

z_pc_vals <- sort(c(as.numeric(params_def2["z_pc"]),
                    10^seq(-18, -7, length.out = 50)))

z_nb_vals <- z_pb_vals * np_zb_ratio_vir
z_nc_vals <- z_pc_vals * np_zc_ratio_vir

# data frame
z_b_in <- tibble(param_foc1 = "z_nb",
                 param_val1 = z_nb_vals,
                 param_foc2 = "z_pb",
                 param_val2 = z_pb_vals)

z_c_in <- tibble(param_foc1 = "z_nc",
                 param_val1 = z_nc_vals,
                 param_foc2 = "z_pc",
                 param_val2 = z_pc_vals)

# PAV simulation
pdf("output/ode_simulations_pav_1st_z.pdf")
pav_1st_z <- z_b_in %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2) param_fun(params_def2, param_foc1 = param_foc1, param_val1 = param_val1, param_foc2 = param_foc2, param_val2 = param_val2, first_virus = "PAV", output_type = 23, inv_time = 1, V0_c = 0))) %>%
  unnest(cols = c(sim_out))
dev.off()

# RPV simulation
pdf("output/ode_simulations_rpv_1st_z.pdf")
rpv_1st_z <- z_c_in %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2) param_fun(params_def2, param_foc1 = param_foc1, param_val1 = param_val1, param_foc2 = param_foc2, param_val2 = param_val2, first_virus = "RPV", output_type = 23, inv_time = 1, V0_b = 0))) %>%
  unnest(cols = c(sim_out))
dev.off()

# save datasets
write_csv(pav_1st_z, "output/ode_simulations_pav_1st_z.csv")
write_csv(rpv_1st_z, "output/ode_simulations_rpv_1st_z.csv")

# re-import datasets
pav_1st_z <- read_csv("output/ode_simulations_pav_1st_z.csv") %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P", "+N+P"))
rpv_1st_z <- read_csv("output/ode_simulations_rpv_1st_z.csv") %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P", "+N+P"))

# extract nutrient info
pav_z_Q <- pav_1st_z %>%
  filter(variable2 == "Q_n") %>%
  select(param_val1, value2, nutrient) %>%
  rename(virus_z = param_val1,
         host_q = value2) %>%
  mutate(internal_nutrient = "N") %>%
  full_join(pav_1st_z %>%
              filter(variable2 == "Q_p") %>%
              select(param_val1, value2, nutrient) %>%
              rename(virus_z = param_val1,
                     host_q = value2) %>%
              mutate(internal_nutrient = "P"))

rpv_z_Q <- rpv_1st_z %>%
  filter(variable2 == "Q_n") %>%
  select(param_val1, value2, nutrient) %>%
  rename(virus_z = param_val1,
         host_q = value2) %>%
  mutate(internal_nutrient = "N") %>%
  full_join(rpv_1st_z %>%
              filter(variable2 == "Q_p") %>%
              select(param_val1, value2, nutrient) %>%
              rename(virus_z = param_val1,
                     host_q = value2) %>%
              mutate(internal_nutrient = "P"))

# figure
pav_z_fig <- ggplot(pav_z_Q, aes(x = virus_z, y = host_q, 
                                 color = nutrient, linetype = nutrient,
                                 size = internal_nutrient)) +
  geom_line() +
  scale_x_log10() +
  scale_color_viridis_d(direction = -1, name = "Nutrient\nsupply") +
  scale_linetype(name = "Nutrient\nsupply") +
  scale_size_manual(values = c(1.2, 0.6), name = "Internal\nnutrient")  +
  labs(x = "BYDV-PAV N use", 
       y = "Plant nutrient concentration") +
  fig_theme

rpv_z_fig <- ggplot(rpv_z_Q, aes(x = virus_z, y = host_q, 
                                 color = nutrient, linetype = nutrient,
                                 size = internal_nutrient)) +
  geom_line() +
  scale_x_log10() +
  scale_color_viridis_d(direction = -1, name = "Nutrient\nsupply") +
  scale_linetype(name = "Nutrient\nsupply") +
  scale_size_manual(values = c(1.2, 0.6), name = "Internal\nnutrient")  +
  labs(x = "CYDV-RPV N use", 
       y = "Plant nutrient concentration") +
  fig_theme

# repeat simulations for virus growth rate
# PAV simulation
pdf("output/ode_simulations_pav_growth_z.pdf")
pav_growth_z <- z_b_in %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2) param_fun(params_def2, param_foc1 = param_foc1, param_val1 = param_val1, param_foc2 = param_foc2, param_val2 = param_val2, first_virus = "PAV", output_type = "growth", inv_time = 1, V0_c = 0))) %>%
  unnest(cols = c(sim_out))
dev.off()

# RPV simulation
pdf("output/ode_simulations_rpv_growth_z.pdf")
rpv_growth_z <- z_c_in %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2) param_fun(params_def2, param_foc1 = param_foc1, param_val1 = param_val1, param_foc2 = param_foc2, param_val2 = param_val2, first_virus = "RPV", output_type = "growth", inv_time = 1, V0_b = 0))) %>%
  unnest(cols = c(sim_out))
dev.off()

# save datasets
write_csv(pav_growth_z, "output/ode_simulations_pav_growth_z.csv")
write_csv(rpv_growth_z, "output/ode_simulations_rpv_growth_z.csv")

# re-import datasets
pav_growth_z <- read_csv("output/ode_simulations_pav_growth_z.csv") %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P", "+N+P"))
rpv_growth_z <- read_csv("output/ode_simulations_rpv_growth_z.csv") %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P", "+N+P"))

# positive growth rates (for other simulation)
pav_growth_pos <- pav_growth_z  %>%
  filter(variable2 == "PAV_conc" & growth > 0)
rpv_growth_pos <- rpv_growth_z  %>%
  filter(variable2 == "RPV_conc" & growth > 0)

# extract growth info
pav_z_growth <- pav_growth_z %>%
  filter(variable2 == "PAV_conc") %>%
  select(param_val1, growth, nutrient) %>%
  rename(virus_z = param_val1)

rpv_z_growth <- rpv_growth_z %>%
  filter(variable2 == "RPV_conc") %>%
  select(param_val1, growth, nutrient) %>%
  rename(virus_z = param_val1)

# figure
pav_z_growth_fig <- ggplot(pav_z_growth, aes(x = virus_z, y = growth)) +
  geom_hline(yintercept = 0, size = 0.25) +
  geom_line(aes(color = nutrient, linetype = nutrient),
            size = 0.9) +
  scale_x_log10() +
  scale_color_viridis_d(direction = -1, name = "Nutrient\nsupply") +
  scale_linetype(name = "Nutrient\nsupply") +
  scale_size_manual(values = c(1.2, 0.6), name = "Internal\nnutrient")  +
  labs(x = "BYDV-PAV N use", 
       y = "BYDV-PAV growth rate") +
  fig_theme

rpv_z_growth_fig <- ggplot(rpv_z_growth, aes(x = virus_z, y = growth)) +
  geom_hline(yintercept = 0, size = 0.25) +
  geom_line(aes(color = nutrient, linetype = nutrient),
            size = 0.9) +
  scale_x_log10() +
  scale_color_viridis_d(direction = -1, name = "Nutrient\nsupply") +
  scale_linetype(name = "Nutrient\nsupply") +
  scale_size_manual(values = c(1.2, 0.6), name = "Internal\nnutrient")  +
  labs(x = "CYDV-RPV N use", 
       y = "CYDV-RPV growth rate") +
  fig_theme

z_growth_fig <- pav_z_growth_fig + theme(legend.position = "none") +
  rpv_z_growth_fig + theme(legend.spacing.y = unit(0.01, "cm")) +
  plot_layout(nrow = 1) + 
  plot_annotation(tag_levels = 'a',
                  tag_prefix = "(",
                  tag_suffix = ")") & 
  theme(plot.tag = element_text(size = 8, face = "bold"))

# Figure S3
ggsave("output/single_virus_z_growth.pdf", z_growth_fig,
       width = 6, height = 2.5, units = "in")


#### minimum nutrient concentrations ####

# import default simulation
plant_def <- read_csv("output/short_term_plant_simulation.csv") %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P", "+N+P"))

# plant nutrient content at time of invasion
q_thresh <- plant_def %>%
  filter(time == plant_days + res_days & 
           variable2 %in% c("Q_n", "Q_p", "Qlim")) %>%
  select(-variable) %>%
  pivot_wider(values_from = "value",
              names_from = "variable2") %>%
  mutate(q_nb_thresh = Q_n * (1 - ((1-Qlim) * 
                                     as.numeric(params_def2["g"]) + 
                                     as.numeric(params_def2["c_b"])) /
                                as.numeric(params_def2["r_b"])),
         q_nc_thresh = Q_n * (1 - ((1-Qlim) * 
                                     as.numeric(params_def2["g"]) + 
                                     as.numeric(params_def2["c_c"])) /
                                as.numeric(params_def2["r_c"])),
         q_pb_thresh = Q_p * (1 - ((1-Qlim) * 
                                     as.numeric(params_def2["g"]) + 
                                     as.numeric(params_def2["c_b"])) /
                                as.numeric(params_def2["r_b"])),
         q_pc_thresh = Q_p * (1 - ((1-Qlim) * 
                                     as.numeric(params_def2["g"]) + 
                                     as.numeric(params_def2["c_c"])) /
                                as.numeric(params_def2["r_c"])),
         growth = 0)

# thresholds calculated under the assumption the focal
# nutrient is limiting virus replication
q_thresh %>%
  mutate(Q_nb_lim = if_else(q_nb_thresh/Q_n > 1e-6/Q_p, "yes", "no"),
         Q_pb_lim = if_else(q_pb_thresh/Q_p > 1e-6/Q_n, "yes", "no"),
         Q_nc_lim = if_else(q_nc_thresh/Q_n > 1e-6/Q_p, "yes", "no"),
         Q_pc_lim = if_else(q_pc_thresh/Q_p > 1e-6/Q_n, "yes", "no"))
# all are "yes"

# values
q_nb_vals <- q_pb_vals <- sort(c(as.numeric(params_def2["q_nb"]),
                                 as.numeric(params_def2["q_pb"]),
                                 unique(q_thresh$q_nb_thresh),
                                 unique(q_thresh$q_pb_thresh),
                                 1e-6, 5e-6, 1e-5, 2.5e-5, 5e-5, 3e-4, 5e-4, 7e-4, 3e-3, 5e-3, 7e-3, 1e-2))

q_nc_vals <- q_pc_vals <- sort(c(as.numeric(params_def2["q_nc"]),
                                 as.numeric(params_def2["q_pc"]),
                                 unique(q_thresh$q_nc_thresh),
                                 unique(q_thresh$q_pc_thresh),
                                 1e-6, 5e-6, 1e-5, 2.5e-5, 5e-5, 3e-4, 6e-4, 9e-4, 3e-3, 5e-3, 7e-3, 1e-2))

# data frame
q_b_in <- tibble(param_foc1 = "q_nb",
                 param_val1 = q_nb_vals,
                 param_foc2 = "q_pb",
                 param_val2 = 1e-6) %>%
  full_join(tibble(param_foc2 = "q_pb",
                   param_val2 = q_pb_vals,
                   param_foc1 = "q_nb",
                   param_val1 = 1e-6))

q_c_in <- tibble(param_foc1 = "q_nc",
                 param_val1 = q_nc_vals,
                 param_foc2 = "q_pc",
                 param_val2 = 1e-6) %>%
  full_join(tibble(param_foc2 = "q_pc",
                   param_val2 = q_pc_vals,
                   param_foc1 = "q_nc",
                   param_val1 = 1e-6))

# PAV simulation
pdf("output/ode_simulations_pav_2nd_q.pdf")
pav_2nd_q <- q_b_in %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2) param_fun(params_def2, param_foc1 = param_foc1, param_val1 = param_val1, param_foc2 = param_foc2, param_val2 = param_val2, first_virus = "RPV", output_type = "growth", V0_c = 0))) %>%
  unnest(cols = c(sim_out))
dev.off()

# RPV simulation
pdf("output/ode_simulations_rpv_2nd_q.pdf")
rpv_2nd_q <- q_c_in %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2) param_fun(params_def2, param_foc1 = param_foc1, param_val1 = param_val1, param_foc2 = param_foc2, param_val2 = param_val2, first_virus = "PAV", output_type = "growth", V0_b = 0))) %>%
  unnest(cols = c(sim_out))
dev.off()

# save datasets
write_csv(pav_2nd_q, "output/ode_simulations_pav_2nd_q.csv")
write_csv(rpv_2nd_q, "output/ode_simulations_rpv_2nd_q.csv")

# re-import datasets
pav_2nd_q <- read_csv("output/ode_simulations_pav_2nd_q.csv") %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P", "+N+P"))
rpv_2nd_q <- read_csv("output/ode_simulations_rpv_2nd_q.csv") %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P", "+N+P"))

# preliminary figure
pav_2nd_q %>%
  filter(variable2 == "PAV_conc") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = growth)) +
  geom_point(size = 10, shape = 15) +
  facet_wrap(~ nutrient) +
  scale_color_viridis_c(name = "Growth rate", direction = -1) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Min N conc.", y = "Min P conc.")

rpv_2nd_q %>%
  filter(variable2 == "RPV_conc") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = growth)) +
  geom_point(size = 10, shape = 15) +
  facet_wrap(~ nutrient) +
  scale_color_viridis_c(name = "Growth rate", direction = -1) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Min N conc.", y = "Min P conc.")

# extract growth info
pav_q_growth <- pav_2nd_q %>%
  filter(variable2 == "PAV_conc" & param_val2 == min(param_val2)) %>%
  select(param_val1, growth, nutrient) %>%
  rename(virus_q = param_val1) %>%
  mutate(internal_nutrient = "N") %>%
  full_join(pav_2nd_q %>%
              filter(variable2 == "PAV_conc" & param_val1 == min(param_val1)) %>%
              select(param_val2, growth, nutrient) %>%
              rename(virus_q = param_val2) %>%
              mutate(internal_nutrient = "P"))

rpv_q_growth <- rpv_2nd_q %>%
  filter(variable2 == "RPV_conc" & param_val2 == min(param_val2)) %>%
  select(param_val1, growth, nutrient) %>%
  rename(virus_q = param_val1) %>%
  mutate(internal_nutrient = "N") %>%
  full_join(rpv_2nd_q %>%
              filter(variable2 == "RPV_conc" & param_val1 == min(param_val1)) %>%
              select(param_val2, growth, nutrient) %>%
              rename(virus_q = param_val2) %>%
              mutate(internal_nutrient = "P"))

# make q_thresh long
q_thresh2 <- q_thresh %>%
  pivot_longer(cols = q_nb_thresh:q_pc_thresh,
               names_to = c("internal_nutrient", "virus"),
               names_pattern = "q_(.)(.)_thresh",
               values_to = "q_thresh") %>%
  mutate(internal_nutrient = toupper(internal_nutrient))

pav_q_thresh <- q_thresh2 %>%
  filter(virus == "b")

rpv_q_thresh <- q_thresh2 %>%
  filter(virus == "c")

# figure
pav_q_fig <- ggplot(pav_q_growth, aes(x = virus_q, y = growth)) +
  geom_hline(yintercept = 0, size = 0.25) +
  geom_line(aes(color = nutrient, linetype = nutrient,
                size = internal_nutrient)) +
  geom_point(data = pav_q_thresh,
             aes(x = q_thresh, color = nutrient, shape = internal_nutrient),
             size = 2) +
  scale_color_viridis_d(direction = -1, name = "Nutrient\nsupply") +
  scale_fill_viridis_d(direction = -1, name = "Nutrient\nsupply") +
  scale_linetype(name = "Nutrient\nsupply") +
  scale_size_manual(values = c(1.2, 0.6), name = "Internal\nnutrient")  +
  scale_shape_manual(values = c(19, 15), name = "Internal\nnutrient")  +
  scale_x_log10(labels = scientific_10) +
  labs(x = "Min. nutrient conc. for replication",
       y = "BYDV-PAV growth rate") +
  fig_theme +
  guides(shape = guide_legend(override.aes = list(color = "black")),
         color = guide_legend(override.aes = list(shape = NA)))

rpv_q_fig <- ggplot(rpv_q_growth, aes(x = virus_q, y = growth)) +
  geom_hline(yintercept = 0, size = 0.25) +
  geom_line(aes(color = nutrient, linetype = nutrient,
                size = internal_nutrient)) +
  geom_point(data = rpv_q_thresh,
             aes(x = q_thresh, color = nutrient, shape = internal_nutrient),
             size = 2) +
  scale_color_viridis_d(direction = -1, name = "Nutrient\nsupply") +
  scale_fill_viridis_d(direction = -1, name = "Nutrient\nsupply") +
  scale_linetype(name = "Nutrient\nsupply") +
  scale_size_manual(values = c(1.2, 0.6), name = "Internal\nnutrient")  +
  scale_shape_manual(values = c(19, 15), name = "Internal\nnutrient")  +
  scale_x_log10(labels = scientific_10) +
  labs(x = "Min. nutrient conc. for replication",
       y = "CYDV-RPV growth rate") +
  fig_theme +
  guides(shape = guide_legend(override.aes = list(color = "black")),
         color = guide_legend(override.aes = list(shape = NA)))


#### combine figures ####

# get legend
leg1 <- get_legend(pav_dpp_fig)
leg2 <- get_legend(pav_z_fig +
                     guides(color = "none",
                            fill = "none",
                            linetype = "none"))
leg3 <- get_legend(pav_q_fig +
                     guides(color = "none",
                            fill = "none",
                            linetype = "none",
                            size = "none"))

# remove legends
pav_dpp_fig2 <- pav_dpp_fig + theme(legend.position = "none")
rpv_dpp_fig2 <- rpv_dpp_fig + theme(legend.position = "none")
pav_z_fig2 <- pav_z_fig + theme(legend.position = "none")
rpv_z_fig2 <- rpv_z_fig + theme(legend.position = "none")
pav_q_fig2 <- pav_q_fig + theme(legend.position = "none")
rpv_q_fig2 <- rpv_q_fig + theme(legend.position = "none")

comb_fig <- pav_dpp_fig2 + rpv_dpp_fig2 + leg1 + 
  pav_z_fig2 + rpv_z_fig2 + leg2 + 
  pav_q_fig2 + rpv_q_fig2 + leg3 +
  plot_layout(ncol = 3, widths = c(1, 1, 0.2), tag_level = 'new') + 
  plot_annotation(tag_levels = list(c("(a)", "(b)", "", "(c)", "(d)", "", "(e)", "(f)", ""))) & 
  theme(plot.tag = element_text(size = 8, face = "bold"))

# Figure S2
ggsave("output/single_virus_ode_simulations.pdf", comb_fig,
       width = 18, height = 18, units = "cm")

# one virus figure
rpv_dpp_fig3 <- rpv_dpp_fig2 %+% labs(y = "Virus growth rate")
rpv_z_fig3 <- rpv_z_fig2 %+% labs(x = "Virus N use")
rpv_q_fig3 <- rpv_q_fig2 %+% labs(y = "Virus growth rate")

sing_comb_fig <- rpv_dpp_fig3 + leg1 + 
  rpv_z_fig3 + leg2 + 
  rpv_q_fig3 + leg3 +
  plot_layout(ncol = 2, widths = c(1, 0.2), tag_level = 'new') + 
  plot_annotation(tag_levels = list(c("(a)", "", "(b)", "", "(c)", ""))) & 
  theme(plot.tag = element_text(size = 8, face = "bold"))

# Figure 2
ggsave("output/one_col_single_virus_ode_simulations.pdf", sing_comb_fig,
       width = 8.5, height = 18, units = "cm")


#### virus nutrient use and minimum nutrient concentrations ####

# extract z and q values
# only use z values where resident growth is positive

# PAV is resident, RPV invader
z_q_b_in <- pav_growth_pos %>%
  select(param_foc1, param_val1, param_foc2, param_val2) %>%
  unique() %>%
  expand_grid(q_c_in %>%
                rename(param_foc3 = param_foc1,
                       param_val3 = param_val1,
                       param_foc4 = param_foc2,
                       param_val4 = param_val2)) %>%
  filter(param_val1 >= 1e-13)

# RPV is resident, PAV invader
z_q_c_in <- rpv_growth_pos %>%
  select(param_foc1, param_val1, param_foc2, param_val2) %>%
  unique() %>%
  expand_grid(q_b_in %>%
                rename(param_foc3 = param_foc1,
                       param_val3 = param_val1,
                       param_foc4 = param_foc2,
                       param_val4 = param_val2)) %>%
  filter(param_val1 >= 1e-13)

# PAV simulation
pdf("output/ode_simulations_pav_res_q_z.pdf")
pav_res_qz <- z_q_b_in %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2, param_foc3, param_val3, param_foc4, param_val4) param_fun(params_def2, param_foc1 = param_foc1, param_val1 = param_val1, param_foc2 = param_foc2, param_val2 = param_val2, param_foc3 = param_foc3, param_val3 = param_val3, param_foc4 = param_foc4, param_val4 = param_val4, first_virus = "PAV", output_type = "23", inv_time = 3))) %>%
  unnest(cols = c(sim_out))
dev.off()

# RPV simulation
pdf("output/ode_simulations_rpv_res_q_z.pdf")
rpv_res_qz <- z_q_c_in %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2, param_foc3, param_val3, param_foc4, param_val4) param_fun(params_def2, param_foc1 = param_foc1, param_val1 = param_val1, param_foc2 = param_foc2, param_val2 = param_val2, param_foc3 = param_foc3, param_val3 = param_val3, param_foc4 = param_foc4, param_val4 = param_val4, first_virus = "RPV", output_type = "23", inv_time = 3))) %>%
  unnest(cols = c(sim_out))
dev.off()

# save datasets
write_csv(pav_res_qz, "output/ode_simulations_pav_res_q_z.csv")
write_csv(rpv_res_qz, "output/ode_simulations_rpv_res_q_z.csv")

# re-import datasets
pav_res_qz <- read_csv("output/ode_simulations_pav_res_q_z.csv") %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P", "+N+P"))
rpv_res_qz <- read_csv("output/ode_simulations_rpv_res_q_z.csv") %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P", "+N+P"))

# add threshold qs
pav_res_qz_thresh <- pav_res_qz %>%
  left_join(q_thresh %>%
              select(nutrient, ends_with("_thresh"))) %>%
  mutate(thresh_q = if_else(variable2 == "RPV_inv_gr" & 
                              (round_half_up(param_val3, 6) == round_half_up(q_nc_thresh, 6) | 
                                 round_half_up(param_val4, 6) == round_half_up(q_pc_thresh, 6)),
                            1, 0))

rpv_res_qz_thresh <- rpv_res_qz %>%
  left_join(q_thresh %>%
              select(nutrient, ends_with("_thresh"))) %>%
  mutate(thresh_q = if_else(variable2 == "PAV_inv_gr" &
                              (round_half_up(param_val3, 6) == round_half_up(q_nb_thresh, 6) | 
                                 round_half_up(param_val4, 6) == round_half_up(q_pb_thresh, 6)),
                            1, 0))

# figures
rpv_inv_n_fig <- pav_res_qz_thresh %>%
  filter(variable2 == "RPV_inv_gr" & param_val4 == min(param_val4)) %>%
  ggplot(aes(x = param_val1, y = value2)) +
  geom_line(aes(color = param_val3,
                group = as.factor(param_val3),
                linewidth = as.factor(thresh_q))) +
  scale_color_viridis_c(name = "Invader min.\nN conc.", trans = "log",
                        breaks = c(0.01, 1e-4, 1e-6),
                        direction = -1, begin = 0.1) +
  scale_linewidth_manual(values = c(0.5, 1.5),
                         guide = "none") +
  facet_wrap(~ nutrient) +
  scale_x_log10() +
  labs(x = "BYDV-PAV (resident) N use", y = "CYDV-RPV (invader) growth rate") +
  fig_theme

rpv_inv_p_fig <- pav_res_qz_thresh %>%
  filter(variable2 == "RPV_inv_gr" & param_val3 == min(param_val3)) %>%
  ggplot(aes(x = param_val2, y = value2)) +
  geom_line(aes(color = param_val4,
                group = as.factor(param_val4),
                linewidth = as.factor(thresh_q))) +
  scale_color_viridis_c(name = "Invader min.\nP conc.", trans = "log",
                        breaks = c(0.01, 1e-4, 1e-6),
                        direction = -1, begin = 0.1) +
  scale_linewidth_manual(values = c(0.5, 1.5),
                         guide = "none") +
  facet_wrap(~ nutrient) +
  scale_x_log10() +
  labs(x = "BYDV-PAV (resident) P use", y = "CYDV-RPV (invader) growth rate") +
  fig_theme

pav_inv_n_fig <- rpv_res_qz_thresh %>%
  filter(variable2 == "PAV_inv_gr" & param_val4 == min(param_val4)) %>%
  ggplot(aes(x = param_val1, y = value2)) +
  geom_line(aes(color = param_val3,
                group = as.factor(param_val3),
                linewidth = as.factor(thresh_q))) +
  scale_color_viridis_c(name = "Invader min.\nN conc.", trans = "log",
                        breaks = c(0.01, 1e-4, 1e-6),
                        direction = -1, begin = 0.1) +
  scale_linewidth_manual(values = c(0.5, 1.5),
                         guide = "none") +
  facet_wrap(~ nutrient) +
  scale_x_log10() +
  labs(x = "CYDV-RPV (resident) N use", y = "BYDV-PAV (invader) growth rate") +
  fig_theme

pav_inv_p_fig <- rpv_res_qz_thresh %>%
  filter(variable2 == "PAV_inv_gr" & param_val3 == min(param_val3)) %>%
  ggplot(aes(x = param_val2, y = value2)) +
  geom_line(aes(color = param_val4,
                group = as.factor(param_val4),
                linewidth = as.factor(thresh_q))) +
  scale_color_viridis_c(name = "Invader min.\nP conc.", trans = "log",
                        breaks = c(0.01, 1e-4, 1e-6),
                        direction = -1, begin = 0.1) +
  scale_linewidth_manual(values = c(0.5, 1.5),
                         guide = "none") +
  facet_wrap(~ nutrient) +
  scale_x_log10() +
  labs(x = "CYDV-RPV (resident) P use", y = "BYDV-PAV (invader) growth rate") +
  fig_theme

# Figure S4
ggsave("output/pav_invasion_n_ode_simulations.pdf", pav_inv_n_fig,
       width = 6.5, height = 5, units = "in")

# Figure S5
ggsave("output/pav_invasion_p_ode_simulations.pdf", pav_inv_p_fig,
       width = 6.5, height = 5, units = "in")

# Figure S6
ggsave("output/rpv_invasion_n_ode_simulations.pdf", rpv_inv_n_fig,
       width = 6.5, height = 5, units = "in")

# Figure S7
ggsave("output/rpv_invasion_p_ode_simulations.pdf", rpv_inv_p_fig,
       width = 6.5, height = 5, units = "in")

# figure with plant nutrients
rpv_inv_n_gr_fig <- pav_res_qz_thresh %>%
  filter(variable2 == "RPV_inv_gr" & thresh_q == 1 & 
           param_val4 == min(param_val4)) %>%
  ggplot(aes(x = param_val1, y = value2)) +
  geom_line(aes(color = nutrient,
                linetype = nutrient),
            linewidth = 0.9) +
  scale_color_viridis_d(direction = -1, name = "Nutrient\nsupply") +
  scale_linetype(name = "Nutrient\nsupply") +
  scale_x_log10() +
  labs(x = "BYDV-PAV (resident) N use", y = "CYDV-RPV (invader) growth rate") +
  fig_theme

rpv_inv_n_plant_fig <- rpv_inv_n_gr_fig %+%
  distinct(filter(pav_res_qz_thresh, variable2 == "Q_n" & 
                    param_val4 == min(param_val4)),
           param_val1, value2, nutrient) +
  labs(y = "Plant N concentration") +
  fig_theme +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())

rpv_inv_p_gr_fig <- pav_res_qz_thresh %>%
  filter(variable2 == "RPV_inv_gr" & thresh_q == 1 & 
           param_val3 == min(param_val3)) %>%
  ggplot(aes(x = param_val2, y = value2)) +
  geom_line(aes(color = nutrient,
                linetype = nutrient),
            linewidth = 0.9) +
  scale_color_viridis_d(direction = -1, name = "Nutrient\nsupply") +
  scale_linetype(name = "Nutrient\nsupply") +
  scale_x_log10() +
  labs(x = "BYDV-PAV (resident) P use", y = "CYDV-RPV (invader) growth rate") +
  fig_theme

rpv_inv_p_plant_fig <- rpv_inv_p_gr_fig %+%
  distinct(filter(pav_res_qz_thresh, variable2 == "Q_p" & 
                    param_val3 == min(param_val3)),
           param_val2, value2, nutrient) +
  labs(y = "Plant P concentration") +
  fig_theme +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())

pav_inv_n_gr_fig <- rpv_res_qz_thresh %>%
  filter(variable2 == "PAV_inv_gr" & thresh_q == 1 & 
           param_val4 == min(param_val4)) %>%
  ggplot(aes(x = param_val1, y = value2)) +
  geom_line(aes(color = nutrient,
                linetype = nutrient),
            linewidth = 0.9) +
  scale_color_viridis_d(direction = -1, name = "Nutrient\nsupply") +
  scale_linetype(name = "Nutrient\nsupply") +
  scale_x_log10() +
  labs(x = "CYDV-RPV (resident) N use", y = "BYDV-PAV (invader) growth rate") +
  fig_theme

pav_inv_n_plant_fig <- pav_inv_n_gr_fig %+%
  distinct(filter(rpv_res_qz_thresh, variable2 == "Q_n" & 
                    param_val4 == min(param_val4)),
           param_val1, value2, nutrient) +
  labs(y = "Plant N concentration") +
  fig_theme +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())

pav_inv_p_gr_fig <- rpv_res_qz_thresh %>%
  filter(variable2 == "PAV_inv_gr" & thresh_q == 1 & 
           param_val3 == min(param_val3)) %>%
  ggplot(aes(x = param_val2, y = value2)) +
  geom_line(aes(color = nutrient,
                linetype = nutrient),
            linewidth = 0.9) +
  scale_color_viridis_d(direction = -1, name = "Nutrient\nsupply") +
  scale_linetype(name = "Nutrient\nsupply") +
  scale_x_log10() +
  labs(x = "CYDV-RPV (resident) P use", y = "BYDV-PAV (invader) growth rate") +
  fig_theme

pav_inv_p_plant_fig <- pav_inv_p_gr_fig %+%
  distinct(filter(rpv_res_qz_thresh, variable2 == "Q_p" & 
                    param_val3 == min(param_val3)),
           param_val2, value2, nutrient) +
  labs(y = "Plant P concentration") +
  fig_theme +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())

# combine
rpv_inv_fig <- rpv_inv_n_plant_fig +
  rpv_inv_p_plant_fig +
  rpv_inv_n_gr_fig +
  rpv_inv_p_gr_fig +
  plot_layout(ncol = 2, guides = "collect") + 
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")") & 
  theme(plot.tag = element_text(size = 8, face = "bold"))

pav_inv_fig <- pav_inv_n_plant_fig + theme(axis.text.x = element_text(hjust = 0.7)) +
  pav_inv_p_plant_fig +
  pav_inv_n_gr_fig + theme(axis.text.x = element_text(hjust = 0.7)) +
  pav_inv_p_gr_fig +
  plot_layout(ncol = 2, guides = "collect") + 
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")") & 
  theme(plot.tag = element_text(size = 8, face = "bold"))

# Figure 3
ggsave("output/rpv_invasion_ode_simulations.pdf", rpv_inv_fig,
       width = 18, height = 13.5, units = "cm")

# Figure S8
ggsave("output/pav_invasion_ode_simulations.pdf", pav_inv_fig,
       width = 18, height = 13.5, units = "cm")
