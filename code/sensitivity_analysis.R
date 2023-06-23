## Goal: sensitivity analysis of virus parameters

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(deSolve)
library(cowplot)
library(lemon)
library(janitor)
library(patchwork)
library(pals)

# load model and settings
source("code/model_settings.R")


#### model functions ####

# change parameter value
param_fun <- function(params_in, 
                      param_foc1 = NA_character_, param_val1 = NA_real_, 
                      param_foc2 = NA_character_, param_val2 = NA_real_,
                      param_foc3 = NA_character_, param_val3 = NA_real_, 
                      param_foc4 = NA_character_, param_val4 = NA_real_,
                      first_virus, output_type = "growth",
                      V0_b = V0, V0_c = V0, inits_in = init_def2,
                      plant_time = 11, res_time = 12, inv_time = 37){
  
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
  
  # save time series
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


#### test functions ####

# test function
virus2_model_sim(params_def2, "RPV", V0_b = V0, V0_c = V0,
                 plant_time = 11, res_time = 12, 
                 inv_time = 37) %>%
  virus2_model_format(params_def2) %>%
  filter(time == max(time) & variable2 %in% c("PAV_conc", "PAV_pop")) %>%
  mutate(virus_conc = log10(value + 1))

pdf("output/temp_sensitivity_analysis_fig.pdf")
param_fun(params_def2, 
          param_foc1 = "q_nb", param_val1 = 1.1e-3, 
          param_foc2 = "q_pb", param_val2 = 7.4e-5, 
          first_virus = "RPV", output_type = 60) %>%
  filter(variable2 %in% c("PAV_conc", "PAV_pop"))
dev.off()
# output should match output of above

pdf("output/temp_sensitivity_analysis_fig.pdf")
param_fun(params_def2, 
          param_foc1 = "a_n_lo", param_val1 = 1.1e-10, 
          first_virus = "RPV", output_type = 60) %>%
  filter(variable2 %in% c("PAV_conc", "PAV_pop"))
dev.off()
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
pdf("output/sensitivity_analysis_pav_1st_dpp.pdf")
pav_1st_dpp <- dpp_in %>%
  mutate(sim_out = pmap(., function(plant_days, res_days, inv_days) param_fun(params_in = params_def2, plant_time = plant_days, res_time = res_days, inv_time = inv_days, first_virus = "PAV", output_type = "growth", V0_c = 0))) %>%
  unnest(cols = c(sim_out))
dev.off()

# RPV alone
pdf("output/sensitivity_analysis_rpv_1st_dpp.pdf")
rpv_1st_dpp <- dpp_in %>%
  mutate(sim_out = pmap(., function(plant_days, res_days, inv_days) param_fun(params_in = params_def2, plant_time = plant_days, res_time = res_days, inv_time = inv_days, first_virus = "RPV", output_type = "growth", V0_b = 0))) %>%
  unnest(cols = c(sim_out))
dev.off()

# save datasets
write_csv(pav_1st_dpp, "output/sensitivity_analysis_pav_1st_dpp.csv")
write_csv(rpv_1st_dpp, "output/sensitivity_analysis_rpv_1st_dpp.csv")

# re-import datasets
pav_1st_dpp <- read_csv("output/sensitivity_analysis_pav_1st_dpp.csv") %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P", "+N+P"))
rpv_1st_dpp <- read_csv("output/sensitivity_analysis_rpv_1st_dpp.csv") %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P", "+N+P"))

# figures
pav_dpp_fig <- pav_1st_dpp %>%
  filter(variable2 == "PAV_conc") %>%
  ggplot(aes(x = plant_days, y = growth, color = nutrient, linetype = nutrient)) +
  geom_hline(yintercept = 0, size = 0.25) +
  geom_line(size = 0.9) +
  scale_color_viridis_d(direction = -1, name = "Nutrient\nsupply") +
  scale_linetype(name = "Nutrient\nsupply") +
  labs(x = "Infection time (DPP)", y = "BYDV-PAV growth rate") +
  fig_theme

rpv_dpp_fig <- rpv_1st_dpp %>%
  filter(variable2 == "RPV_conc") %>%
  ggplot(aes(x = plant_days, y = growth, color = nutrient, linetype = nutrient)) +
  geom_hline(yintercept = 0, size = 0.25) +
  geom_line(size = 0.9) +
  scale_color_viridis_d(direction = -1, name = "Nutrient\nsupply") +
  scale_linetype(name = "Nutrient\nsupply") +
  labs(x = "Infection time (DPP)", y = "CYDV-RPV growth rate") +
  fig_theme

# plant by itself
short_term_plant <- virus2_model_sim(params_def2, "PAV", V0_b = 0, V0_c = 0, 
                                     plant_time = 11, res_time = 12, 
                                     inv_time = 60-11-12) %>%
  virus2_model_format(params_def2)

# combine grwr with limiting nutrient concentration
# this doesn't quite work because time is a continuous variable in short_term_plant
pav_1st_lim <- pav_1st_dpp %>%
  filter(variable2 == "PAV_conc") %>%
  left_join(short_term_plant %>%
              filter(variable2 == "Qlim") %>%
              rename(plant_days = time,
                     Qlim = value) %>%
              select(plant_days, nutrient, Qlim))

pav_1st_lim %>%
  ggplot(aes(x = Qlim, y = growth, color = nutrient, linetype = nutrient)) +
  geom_hline(yintercept = 0, size = 0.25) +
  geom_line() +
  scale_color_viridis_d(direction = -1, name = "Nutrient\nsupply") +
  scale_linetype(name = "Nutrient\nsupply") +
  labs(x = expression(paste("Limiting nutrient ratio (", Q[min], "/Q)", sep = "")), 
       y = "Growth rate") +
  fig_theme
# this makes sense because this is the virus's growth rate too


#### virus nutrient use ####

# virus N:P
np_zb_ratio_vir <- as.numeric(params_def2["z_nb"]) / as.numeric(params_def2["z_pb"])
np_zc_ratio_vir <- as.numeric(params_def2["z_nc"]) / as.numeric(params_def2["z_pc"])

# values
# too high Z values make Q's go negative
z_pb_vals <- sort(c(as.numeric(params_def2["z_pb"]),
                  10^seq(-18, -12, length.out = 7),
                  10^seq(-11, -7, length.out = 15)))

z_pc_vals <- sort(c(as.numeric(params_def2["z_pc"]),
                    10^seq(-18, -12, length.out = 7),
                    10^seq(-11, -7, length.out = 15)))

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
pdf("output/sensitivity_analysis_pav_1st_z.pdf")
pav_1st_z <- z_b_in %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2) param_fun(params_def2, param_foc1 = param_foc1, param_val1 = param_val1, param_foc2 = param_foc2, param_val2 = param_val2, first_virus = "PAV", output_type = 23, inv_time = 1, V0_c = 0))) %>%
  unnest(cols = c(sim_out))
dev.off()

# RPV simulation
pdf("output/sensitivity_analysis_rpv_1st_z.pdf")
rpv_1st_z <- z_c_in %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2) param_fun(params_def2, param_foc1 = param_foc1, param_val1 = param_val1, param_foc2 = param_foc2, param_val2 = param_val2, first_virus = "RPV", output_type = 23, inv_time = 1, V0_b = 0))) %>%
  unnest(cols = c(sim_out))
dev.off()

# save datasets
write_csv(pav_1st_z, "output/sensitivity_analysis_pav_1st_z.csv")
write_csv(rpv_1st_z, "output/sensitivity_analysis_rpv_1st_z.csv")

# re-import datasets
pav_1st_z <- read_csv("output/sensitivity_analysis_pav_1st_z.csv") %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P", "+N+P"))
rpv_1st_z <- read_csv("output/sensitivity_analysis_rpv_1st_z.csv") %>%
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
pdf("output/sensitivity_analysis_pav_growth_z.pdf")
pav_growth_z <- z_b_in %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2) param_fun(params_def2, param_foc1 = param_foc1, param_val1 = param_val1, param_foc2 = param_foc2, param_val2 = param_val2, first_virus = "PAV", output_type = "growth", inv_time = 1, V0_c = 0))) %>%
  unnest(cols = c(sim_out))
dev.off()

# RPV simulation
pdf("output/sensitivity_analysis_rpv_growth_z.pdf")
rpv_growth_z <- z_c_in %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2) param_fun(params_def2, param_foc1 = param_foc1, param_val1 = param_val1, param_foc2 = param_foc2, param_val2 = param_val2, first_virus = "RPV", output_type = "growth", inv_time = 1, V0_b = 0))) %>%
  unnest(cols = c(sim_out))
dev.off()

# save datasets
write_csv(pav_growth_z, "output/sensitivity_analysis_pav_growth_z.csv")
write_csv(rpv_growth_z, "output/sensitivity_analysis_rpv_growth_z.csv")

# re-import datasets
pav_growth_z <- read_csv("output/sensitivity_analysis_pav_growth_z.csv") %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P", "+N+P"))
rpv_growth_z <- read_csv("output/sensitivity_analysis_rpv_growth_z.csv") %>%
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
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 8, face = "bold"))

ggsave("output/single_virus_z_growth_sensitivity_analysis.pdf", z_growth_fig,
       width = 5, height = 2, units = "in")


#### minimum nutrient concentrations ####

# import default simulation
plant_def <- read_csv("output/short_term_plant_simulation.csv") %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P", "+N+P"))

# plant nutrient content at time of invasion
q_thresh <- plant_def %>%
  filter(time == 11 + 12 & variable2 %in% c("Q_n", "Q_p", "Qlim")) %>%
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

# values (simplified)
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
pdf("output/sensitivity_analysis_pav_2nd_q.pdf")
pav_2nd_q <- q_b_in %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2) param_fun(params_def2, param_foc1 = param_foc1, param_val1 = param_val1, param_foc2 = param_foc2, param_val2 = param_val2, first_virus = "RPV", output_type = "growth", V0_c = 0))) %>%
  unnest(cols = c(sim_out))
dev.off()

pdf("output/sensitivity_analysis_rpv_2nd_q.pdf")
rpv_2nd_q <- q_c_in %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2) param_fun(params_def2, param_foc1 = param_foc1, param_val1 = param_val1, param_foc2 = param_foc2, param_val2 = param_val2, first_virus = "PAV", output_type = "growth", V0_b = 0))) %>%
  unnest(cols = c(sim_out))
dev.off()

# save datasets
write_csv(pav_2nd_q, "output/sensitivity_analysis_pav_2nd_q.csv")
write_csv(rpv_2nd_q, "output/sensitivity_analysis_rpv_2nd_q.csv")

# re-import datasets
pav_2nd_q <- read_csv("output/sensitivity_analysis_pav_2nd_q.csv") %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P", "+N+P"))
rpv_2nd_q <- read_csv("output/sensitivity_analysis_rpv_2nd_q.csv") %>%
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
  plot_annotation(tag_levels = list(c("A", "B", "", "C", "D", "", "E", "F", ""))) & 
  theme(plot.tag = element_text(size = 8, face = "bold"))

ggsave("output/single_virus_sensitivity_analysis.pdf", comb_fig,
       width = 6, height = 6.5, units = "in")

# one resident figure
rpv_dpp_fig3 <- rpv_dpp_fig2 %+% labs(y = "Virus growth rate")
rpv_z_fig3 <- rpv_z_fig2 %+% labs(x = "Virus N use")
rpv_q_fig3 <- rpv_q_fig2 %+% labs(y = "Virus growth rate")

sing_comb_fig <- rpv_dpp_fig3 + leg1 + 
  rpv_z_fig3 + leg2 + 
  rpv_q_fig3 + leg3 +
  plot_layout(ncol = 2, widths = c(1, 0.2), tag_level = 'new') + 
  plot_annotation(tag_levels = list(c("A", "", "B", "", "C", ""))) & 
  theme(plot.tag = element_text(size = 8, face = "bold"))

ggsave("output/one_col_single_virus_sensitivity_analysis.pdf", sing_comb_fig,
       width = 3, height = 6.5, units = "in")


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
                       param_val4 = param_val2))

z_q_c_in <- rpv_growth_pos %>%
  select(param_foc1, param_val1, param_foc2, param_val2) %>%
  unique() %>%
  expand_grid(q_b_in %>%
                rename(param_foc3 = param_foc1,
                       param_val3 = param_val1,
                       param_foc4 = param_foc2,
                       param_val4 = param_val2))

# PAV simulation
pdf("output/sensitivity_analysis_pav_res_q_z.pdf")
pav_res_qz <- z_q_b_in %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2, param_foc3, param_val3, param_foc4, param_val4) param_fun(params_def2, param_foc1 = param_foc1, param_val1 = param_val1, param_foc2 = param_foc2, param_val2 = param_val2, param_foc3 = param_foc3, param_val3 = param_val3, param_foc4 = param_foc4, param_val4 = param_val4, first_virus = "PAV", output_type = "23", inv_time = 3))) %>%
  unnest(cols = c(sim_out))
dev.off()

# RPV simulation
pdf("output/sensitivity_analysis_rpv_res_q_z.pdf")
rpv_res_qz <- z_q_c_in %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2, param_foc3, param_val3, param_foc4, param_val4) param_fun(params_def2, param_foc1 = param_foc1, param_val1 = param_val1, param_foc2 = param_foc2, param_val2 = param_val2, param_foc3 = param_foc3, param_val3 = param_val3, param_foc4 = param_foc4, param_val4 = param_val4, first_virus = "RPV", output_type = "23", inv_time = 3))) %>%
  unnest(cols = c(sim_out))
dev.off()

# save datasets -- DO THIS
write_csv(pav_res_qz, "output/sensitivity_analysis_pav_res_q_z.csv")
write_csv(rpv_res_qz, "output/sensitivity_analysis_rpv_res_q_z.csv")

# re-import datasets
pav_res_qz <- read_csv("output/sensitivity_analysis_pav_res_q_z.csv") %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P", "+N+P"))
rpv_res_qz <- read_csv("output/sensitivity_analysis_rpv_res_q_z.csv") %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P", "+N+P"))

# add threshold qs
pav_res_qz_thresh <- pav_res_qz %>%
  left_join(q_thresh %>%
              select(nutrient, ends_with("_thresh"))) %>%
  mutate(thresh_q = if_else(variable2 == "RPV_inv_gr" & 
                              (param_val3 == q_nc_thresh | param_val4 == q_pc_thresh),
                            1, 0)) %>%
  filter(param_val1 >= 1e-13)

rpv_res_qz_thresh <- rpv_res_qz %>%
  left_join(q_thresh %>%
              select(nutrient, ends_with("_thresh"))) %>%
  mutate(thresh_q = if_else(variable2 == "PAV_inv_gr" &
                              (param_val3 == q_nb_thresh | param_val4 == q_pb_thresh),
                            1, 0)) %>%
  filter(param_val1 >= 1e-13)

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

ggsave("output/rpv_invasion_n_sensitivity_analysis.pdf", rpv_inv_n_fig,
       width = 6.5, height = 5, units = "in")
ggsave("output/rpv_invasion_p_sensitivity_analysis.pdf", rpv_inv_p_fig,
       width = 6.5, height = 5, units = "in")
ggsave("output/pav_invasion_n_sensitivity_analysis.pdf", pav_inv_n_fig,
       width = 6.5, height = 5, units = "in")
ggsave("output/pav_invasion_p_sensitivity_analysis.pdf", pav_inv_p_fig,
       width = 6.5, height = 5, units = "in")

#### start here ####
# four panel figure, plant nutrients on top, virus growth rate on bottom
# one resident, one invader
# other resident and invader in supplement
# four figures above in supplement (provide trajectories when q != threshold q)
# figure with plant nutrients
rpv_inv_n_conc_fig <- pav_res_qz_thresh %>%
  filter(variable2 == "RPV_inv_gr" & thresh_q == 1 & 
           param_val4 == min(param_val4)) %>%
  ggplot(aes(x = param_val1, y = value2)) +
  geom_line(aes(color = nutrient,
                linetype = nutrient)) +
  scale_color_viridis_d(direction = -1, name = "Nutrient\nsupply") +
  scale_linetype(name = "Nutrient\nsupply") +
  scale_x_log10() +
  labs(x = "BYDV-PAV (resident) N use", y = "CYDV-RPV (invader) growth rate") +
  fig_theme

rpv_inv_n_plant_fig <- rpv_inv_n_conc_fig %+%
 filter(pav_res_qz_thresh, variable2 == "Q_n" & 
  param_val4 == min(param_val4)) + # lots of reps by q_nc, but all the same
  labs(y = "Plant N concentration") +
  fig_theme +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())

pav_inv_n_conc_fig <- rpv_res_qz_thresh %>%
  filter(variable2 == "PAV_inv_gr" & thresh_q == 1 & 
           param_val4 == min(param_val4)) %>%
  ggplot(aes(x = param_val1, y = value2)) +
  geom_line(aes(color = nutrient,
                linetype = nutrient)) +
  scale_color_viridis_d(direction = -1, name = "Nutrient\nsupply") +
  scale_linetype(name = "Nutrient\nsupply") +
  scale_x_log10() +
  labs(x = "CYDV-RPV (resident) N use", y = "BYDV-PAV (invader) growth rate") +
  fig_theme

pav_inv_n_plant_fig <- pav_inv_n_conc_fig %+%
  filter(rpv_res_qz_thresh, variable2 == "Q_n" & 
           param_val4 == min(param_val4)) + # lots of reps by q_nc, but all the same
  labs(y = "Plant N concentration") +
  fig_theme +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())



#### older code ####

# warm/cold figures
# figures
rpv_inv_n_low_fig <- ggplot(pav_res_qz_n_low, aes(x = param_val1, y = param_val3, fill = value2)) +
  geom_tile(width = 1, height = 1) +
  scale_fill_gradientn(colours = coolwarm(100), limits = c(-0.7, 0.7),
                       name = "CYDV-RPV\ngrowth rate") +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "BYDV-PAV N (resident) use", y = "CYDV-RPV (invader) min. N conc.") +
  fig_theme

rpv_inv_n_n_fig <- rpv_inv_n_low_fig %+% pav_res_qz_n_n
rpv_inv_n_p_fig <- rpv_inv_n_low_fig %+% pav_res_qz_n_p
rpv_inv_n_np_fig <- rpv_inv_n_low_fig %+% pav_res_qz_n_np

rpv_inv_p_low_fig <- ggplot(pav_res_qz_p_low, aes(x = param_val2, y = param_val4, fill = value2)) +
  geom_tile(width = 1, height = 1) +
  scale_fill_gradientn(colours = coolwarm(100), limits = c(-0.7, 0.7),
                       name = "CYDV-RPV\ngrowth rate") +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "BYDV-PAV (resident) P use", y = "CYDV-RPV (invader) min. P conc.") +
  fig_theme

rpv_inv_p_n_fig <- rpv_inv_p_low_fig %+% pav_res_qz_p_n
rpv_inv_p_p_fig <- rpv_inv_p_low_fig %+% pav_res_qz_p_p
rpv_inv_p_np_fig <- rpv_inv_p_low_fig %+% pav_res_qz_p_np

pav_inv_n_low_fig <- ggplot(rpv_res_qz_n_low, aes(x = param_val1, y = param_val3, fill = value2)) +
  geom_tile(width = 1, height = 1) +
  scale_fill_gradientn(colours = coolwarm(100), limits = c(-0.7, 0.7),
                       name = "BYDV-PAV\ngrowth rate") +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "CYDV-RPV (resident) N use", y = "BYDV-PAV (invader) min. N conc.") +
  fig_theme

pav_inv_n_n_fig <- pav_inv_n_low_fig %+% pav_res_qz_n_n
pav_inv_n_p_fig <- pav_inv_n_low_fig %+% pav_res_qz_n_p
pav_inv_n_np_fig <- pav_inv_n_low_fig %+% pav_res_qz_n_np

pav_inv_p_low_fig <- ggplot(rpv_res_qz_p_low, aes(x = param_val2, y = param_val4, fill = value2)) +
  geom_tile(width = 1, height = 1) +
  scale_fill_gradientn(colours = coolwarm(100), limits = c(-0.7, 0.7),
                       name = "BYDV-PAV\ngrowth rate") +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "CYDV-RPV (resident) P use", y = "BYDV-PAV (invader) min. P conc.") +
  fig_theme

pav_inv_p_n_fig <- pav_inv_p_low_fig %+% pav_res_qz_p_n
pav_inv_p_p_fig <- pav_inv_p_low_fig %+% pav_res_qz_p_p
pav_inv_p_np_fig <- pav_inv_p_low_fig %+% pav_res_qz_p_np


#### resource supply rates ####

# values for low nutrient supply
a_n_vals <- c(1.1e-13, 1.1e-10, 1.1e-6, 5.6e-5)
a_p_vals <- c(1.6e-14, 1.6e-11, 1.6e-7, 8.2e-6, 5.6e-5)

# data frame
a_in <- tibble(param_foc1 = "a_n_lo",
               param_val1 = a_n_vals) %>%
  expand_grid(tibble(param_foc2 = "a_p_lo",
                     param_val2 = a_p_vals))

# PAV invades RPV
pdf("output/sensitivity_analysis_pav_inv_a.pdf")
pav_inv_a <- a_in %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2) param_fun(params_def2, param_foc1 = param_foc1, param_val1 = param_val1, param_foc2 = param_foc2, param_val2 = param_val2, first_virus = "RPV", output_type = "first"))) %>%
  unnest(cols = c(sim_out))
dev.off()

# figure
pav_inv_a %>%
  filter(variable2 %in% c("PAV_conc", "RPV_conc") & nutrient == "low") %>%
  mutate(growth_sign = ifelse(growth <= 0, "<= 0", "> 0")) %>%
  ggplot(aes(x = param_val1, y = param_val2, color = growth, shape = growth_sign)) +
  geom_point(size = 10) +
  facet_wrap(~ variable2) +
  scale_color_viridis_c(name = "Growth rate", direction = -1) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "N supply", y = "P supply")
# growth rate increases with N, but not P


#### minimum nutrient concentrations and virus nutrient content ####

# virus N:P
np_z_ratio_vir <- as.numeric(params_def2["z_nb"]) / as.numeric(params_def2["z_pb"])
np_q_ratio_vir <- as.numeric(params_def2["q_nb"]) / as.numeric(params_def2["q_pb"])

# values
z_pb_vals2 <- sort(c(10^c(13, 12, 11, 10, 7),
                     as.numeric(params_def2["z_pb"])))

q_pb_vals2 <- sort(c(10^seq(-6, -2, length.out = 9),
                     as.numeric(params_def2["q_pb"])))

z_nb_vals2 <- z_pb_vals2 * np_z_ratio_vir
q_nb_vals2 <- q_pb_vals2 * np_q_ratio_vir

# data frame
z_q_b_in <- tibble(param_foc1 = "z_nb",
                   param_val1 = z_nb_vals2,
                   param_foc3 = "z_pb",
                   param_val3 = z_pb_vals2) %>%
  expand_grid(tibble(param_foc2 = "q_nb",
                     param_val2 = q_nb_vals2,
                     param_foc4 = "q_pb",
                     param_val4 = q_pb_vals2))

# PAV simulation
pdf("output/sensitivity_analysis_pav_1st_z_q.pdf")
pav_1st_z_q <- z_q_b_in %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2, param_foc3, param_val3, param_foc4, param_val4) param_fun(params_def2, param_foc1 = param_foc1, param_val1 = param_val1, param_foc2 = param_foc2, param_val2 = param_val2, param_foc3 = param_foc3, param_val3 = param_val3, param_foc4 = param_foc4, param_val4 = param_val4, first_virus = "PAV", V0_c = 0, inv_time = 0))) %>%
  unnest(cols = c(sim_out))
dev.off()

write_csv(pav_1st_z_q, "output/sensitivity_analysis_pav_1st_z_q.csv")
pav_1st_z_q <- read_csv("output/sensitivity_analysis_pav_1st_z_q.csv") %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P", "+N+P"))

# figure
pav_1st_z_q %>%
  filter(variable2 == "Q_n") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = value2)) +
  geom_point(size = 10, shape = 15) +
  facet_wrap(~ nutrient) +
  scale_color_viridis_c(name = "Plant N", direction = -1) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "N content", y = "Min N conc.")

# RPV invasion simulation
pdf("output/sensitivity_analysis_rpv_inv_z_q.pdf")
rpv_inv_z_q <- z_q_b_in %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2, param_foc3, param_val3, param_foc4, param_val4) param_fun(params_def2, param_foc1 = param_foc1, param_val1 = param_val1, param_foc2 = param_foc2, param_val2 = param_val2, param_foc3 = param_foc3, param_val3 = param_val3, param_foc4 = param_foc4, param_val4 = param_val4, first_virus = "PAV", output_type = "first"))) %>%
  unnest(cols = c(sim_out))
dev.off()

write_csv(rpv_inv_z_q, "output/sensitivity_analysis_rpv_inv_z_q.csv")
rpv_inv_z_q <- read_csv("output/sensitivity_analysis_rpv_inv_z_q.csv") %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P", "+N+P"))

# figure
rpv_inv_z_q %>%
  filter(variable2 == "RPV_conc") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = growth)) +
  geom_point(size = 10, shape = 15) +
  facet_wrap(~ nutrient) +
  scale_color_viridis_c(name = "Growth rate", direction = -1) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "N content", y = "Min N conc.")

rpv_inv_z_q %>%
  filter(variable2 == "RPV_conc") %>%
  ggplot(aes(x = param_val2, y = growth, color = as.factor(param_val1))) +
  geom_line() +
  facet_wrap(~ nutrient) +
  scale_color_viridis_d(name = "Resident N content", direction = -1) +
  scale_x_log10() +
  labs(x = "Resident min. N conc.", y = "Invader growth rate when rare")

rpv_inv_z_q %>%
  filter(variable2 == "RPV_conc") %>%
  ggplot(aes(x = param_val4, y = growth, color = as.factor(param_val3))) +
  geom_line() +
  facet_wrap(~ nutrient) +
  scale_color_viridis_d(name = "Resident P content", direction = -1) +
  scale_x_log10() +
  labs(x = "Resident min. P conc.", y = "Invader growth rate when rare")

rpv_inv_z_q %>%
  filter(variable2 == "PAV_conc") %>%
  ggplot(aes(x = param_val2, y = growth, color = as.factor(param_val1))) +
  geom_line() +
  facet_wrap(~ nutrient) +
  scale_color_viridis_d(name = "Resident N content", direction = -1) +
  scale_x_log10() +
  labs(x = "Resident min. N conc.", y = "Resident growth rate when rare")


#### minimum nutrient concentration (q_n)/nutrient content of virus (z_n) ####

# values for z and q
z_nc_vals <- 1.7*10^c(-18, -12, -6, -1)
q_nc_vals2 <- 1.1*10^c(-5, -4, -3, -2, -1) 

# data frame
qz_nc_in <- tibble(param_foc1 = "q_nc",
                   param_val1 = q_nc_vals2) %>%
  expand_grid(tibble(param_foc2 = "z_nc",
                   param_val2 = z_nc_vals))

# PAV invades RPV
pdf("output/sensitivity_analysis_pav_inv_qz_nc.pdf")
pav_inv_qz_nc <- qz_nc_in %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2) param_fun(params_def2, param_foc1, param_val1, param_foc2, param_val2, first_virus = "RPV", output_type = "last"))) %>%
  unnest(cols = c(sim_out))
dev.off()

# figure
pav_inv_qz_nc %>%
  filter(variable2 == "RPV_pop") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = value2)) +
  geom_point(size = 10) +
  facet_wrap(~ nutrient) +
  scale_color_viridis_c(name = "RPV pop", direction = -1) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "RPV min N (q)", y = "RPV N conc (z)")
# increasing q decreases growth regardless of z
# increasing z decreases growth
# this is because at low z, the plant controls the nutrients levels, so there are excess
# nutrients for the virus and it can grow a lot, but as z increases, the nutrient levels
# meet the virus's requirements and there aren't as much excess for fast growth
# growth increases when z is set at its highest level and lowest q level, 
# but this is a computation error -- Q_n skips to a negative value despite max functions

pav_inv_qz_nc %>%
  filter(variable2 == "PAV_pop") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = value2)) +
  geom_point(size = 10) +
  facet_wrap(~ nutrient) +
  scale_color_viridis_c(name = "PAV pop", direction = -1) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "RPV min N (q)", y = "RPV N conc (z)")
# increasing RPV's q doesn't affect PAV much, if at all
# increasing RPV's z decreases PAV's growth (until last value, see above)
# increasing RPV's z decreases plant growth
# for high z, decreasing q below the plant's decreases plant and PAV growth

# invasion growth rate
pav_inv_qz_nc2 <- qz_nc_in %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2) param_fun(params_def2, param_foc1, param_val1, param_foc2, param_val2, first_virus = "RPV", output_type = "first"))) %>%
  unnest(cols = c(sim_out))

pav_inv_qz_nc2 %>%
  filter(!(param_val1 == min(param_val1) & param_val2 == max(param_val2))) %>%
  mutate(growth_sign = ifelse(growth <= 0, "<= 0", "> 0")) %>%
  ggplot(aes(x = param_val1, y = param_val2, color = growth, shape = growth_sign)) +
  geom_point(size = 10) +
  facet_grid(variable2 ~ nutrient) +
  scale_color_viridis_c(name = "growth rate", direction = -1) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "RPV min N (q)", y = "RPV N conc (z)")
# PAV has negative growth rate at high z and low q
# increasing RPV's z can cause PAV go to negative before RPV


#### minimum N concentrations (q_n) with higher z_nc ####

# values for z and q
q_nc_vals2 <- 1.1*10^c(-5, -4, -3, -2, -1) # same as above
q_nb_vals2 <- 1.1*10^c(-5, -4, -3, -2, -1)

# data frame
q_n_in2 <- tibble(param_foc1 = "q_nc",
                   param_val1 = q_nc_vals2) %>%
    expand_grid(tibble(param_foc2 = "q_nb",
                     param_val2 = q_nb_vals2))

# set RPV (resident's) z higher
params_q_n <- params_def2
params_q_n["z_nc"] <- 1.7 * 10^-6

# PAV invades RPV
pdf("output/sensitivity_analysis_pav_inv_q_n_high_z.pdf")
pav_inv_q_n2 <- q_n_in2 %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2, first_virus) param_fun(params_q_n, param_foc1, param_val1, param_foc2, param_val2, first_virus = "RPV", output_type = "last"))) %>%
  unnest(cols = c(sim_out))
dev.off()

# figure
pav_inv_q_n2 %>%
  filter(variable2 == "RPV_pop") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = value2)) +
  geom_point(size = 10) +
  facet_wrap(~ nutrient) +
  scale_color_viridis_c(name = "RPV pop", direction = -1) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "RPV min N (q)", y = "PAV min N (q)")
# RPV grows slightly faster with lower q (more excess nutrients)

pav_inv_q_n2 %>%
  filter(variable2 == "PAV_pop") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = value2)) +
  geom_point(size = 10) +
  facet_wrap(~ nutrient) +
  scale_color_viridis_c(name = "PAV pop", direction = -1) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "RPV min N (q)", y = "PAV min N (q)")
# when PAV's q is lower than RPV's, it can grow fast because of excess nutrients
# it grows slowly when its q is equal to or greater than RPV's
# because resources are limited (it's actually only dying here, but c is low)


#### minimum N concentration (q_n)/mortality (c) with higher z_nc ####

# hypothesis: when PAV's q is greater than RPV's and it's c is high, it will decline

# values for z and q
q_nc_vals2 <- 1.1*10^c(-5, -4, -3, -2, -1) # same as above
c_b_vals <- c(0.1, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19)

# data frame
qnc_cb_in <- tibble(param_foc1 = "q_nc",
                       param_val1 = q_nc_vals2) %>%
  expand_grid(tibble(param_foc2 = "c_b",
                     param_val2 = c_b_vals))

# set RPV (resident's) z higher
params_qnc_cb <- params_def2
params_qnc_cb["z_nc"] <- 1.7 * 10^-6

# PAV invades RPV
pdf("output/sensitivity_analysis_pav_inv_qnc_cb_high_znc.pdf")
pav_inv_qnc_cb <- qnc_cb_in %>%
  mutate(first_virus = "RPV") %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2, first_virus) param_fun(params_qnc_cb, param_foc1, param_val1, param_foc2, param_val2, first_virus))) %>%
  unnest(cols = c(sim_out))
dev.off()

# figure
pav_inv_qnc_cb %>%
  filter(variable2 == "RPV_pop") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = value2)) +
  geom_point(size = 10) +
  facet_wrap(~ nutrient) +
  scale_color_viridis_c(name = "RPV pop", direction = -1) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "RPV min N (q)", y = "PAV mortality (c)")
# RPV grows slightly faster with lower q (more excess nutrients)

pav_inv_qnc_cb %>%
  filter(variable2 == "PAV_pop") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = value2)) +
  geom_point(size = 10) +
  facet_wrap(~ nutrient) +
  scale_color_viridis_c(name = "PAV pop", direction = -1) +
  scale_x_log10() +
  labs(x = "RPV min N (q)", y = "PAV mortality (c)")
# a high enough c_b makes PAV decline rather than grow
# decreasing RPV's q can pull down resources more,
# which accelerates PAV's decline
# RPV's q alone can make PAV go from positive to negative growth
# but the change is relatively small with low c_b
# could measure PAV's growth rate next
# growth rate in response to RPV's q (and z)
# should go from positive to negative
# may need to ramp up c_b if the changes are too subtle


#### older code ####

#### q_nc/z_nc -> plant N conc. ####

# values for z and q
z_nc_vals2 <- 10^c(-18, -12, -11, -10, -9, -8, -6, -2)
q_nc_vals3 <- 1.1e-3 * seq(1, 10, length.out = 100)

# data frame
qz_nc_in2 <- tibble(param_foc1 = "q_nc",
                   param_val1 = q_nc_vals3) %>%
  expand_grid(tibble(param_foc2 = "z_nc",
                     param_val2 = z_nc_vals2))

# RPV alone
# pdf("output/sensitivity_analysis_rpv_res_qz_nc.pdf")
# rpv_res_qz_nc <- qz_nc_in2 %>%
#   mutate(first_virus = "RPV") %>%
#   mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2, first_virus) param_fun(params_def2, param_foc1, param_val1, param_foc2, param_val2, first_virus, V0_b = 0))) %>%
#   unnest(cols = c(sim_out))
# dev.off()
# 
# # save simulation output (large)
# write_csv(rpv_res_qz_nc, "output/sensitivity_analysis_rpv_res_qz_nc.csv")
rpv_res_qz_nc <- read_csv("output/sensitivity_analysis_rpv_res_qz_nc.csv")

# figure
rpv_res_qz_nc %>%
  filter(variable2 == "Q_n") %>%
  ggplot(aes(x = param_val1, y = value, color = as.factor(param_val2))) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  geom_hline(yintercept = as.numeric(params_def1["Qmin_n"]), color = "black", linetype = "dashed") +
  geom_line() +
  facet_wrap(~ nutrient) +
  scale_color_viridis_d(name = "RPV N conc (z)") +
  labs(x = "RPV min N (q)", y = "Plant N concentration") +
  fig_theme
#### use N+P panel from this figure ####

# plant biomass
rpv_res_qz_nc %>%
  filter(variable2 == "H") %>%
  ggplot(aes(x = param_val1, y = value, color = as.factor(param_val2))) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  geom_line() +
  facet_wrap(~ nutrient) +
  scale_color_viridis_d(name = "RPV N conc (z)") +
  labs(x = "RPV min N (q)", y = "Plant biomass") +
  fig_theme
# increasing plant N concentration increases plant biomass


#### q_nb/z_nb -> plant N conc. ####

# values for z and q
z_nb_vals <- 10^c(-18, -12, -11, -10, -9, -8, -6, -2)
q_nb_vals2 <- 1.1e-3 * seq(1, 10, length.out = 100)

# data frame
qz_nb_in <- tibble(param_foc1 = "q_nb",
                    param_val1 = q_nb_vals2) %>%
  expand_grid(tibble(param_foc2 = "z_nb",
                     param_val2 = z_nb_vals))

# PAV alone
pdf("output/sensitivity_analysis_pav_res_qz_nb.pdf")
pav_res_qz_nb <- qz_nb_in %>%
  mutate(first_virus = "PAV") %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2, first_virus) param_fun(params_def2, param_foc1, param_val1, param_foc2, param_val2, first_virus, V0_c = 0))) %>%
  unnest(cols = c(sim_out))
dev.off()

# save simulation output (large)
write_csv(pav_res_qz_nb, "output/sensitivity_analysis_pav_res_qz_nb.csv")
pav_res_qz_nb <- read_csv("output/sensitivity_analysis_pav_res_qz_nb.csv")

# plant nutrient conc.
pav_res_qz_nb %>%
  filter(variable2 == "Q_n") %>%
  ggplot(aes(x = param_val1, y = value, color = as.factor(param_val2))) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  geom_line() +
  facet_wrap(~ nutrient) +
  scale_color_viridis_d(name = "PAV N conc (z)") +
  labs(x = "PAV min N (q)", y = "Plant N concentration") +
  fig_theme

# plant biomass
pav_res_qz_nb %>%
  filter(variable2 == "H") %>%
  ggplot(aes(x = param_val1, y = value, color = as.factor(param_val2))) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  geom_line() +
  facet_wrap(~ nutrient) +
  scale_color_viridis_d(name = "PAV N conc (z)") +
  labs(x = "PAV min N (q)", y = "Plant biomass") +
  fig_theme
# increasing plant N concentration increases plant biomass


#### q_pc/z_pc -> plant P conc. ####

# values for z and q
z_pc_vals <- 10^c(-19, -13, -12, -11, -10, -9, -7, -3)
q_pc_vals <- 7.4e-5 * seq(1, 10, length.out = 100)

# data frame
qz_pc_in <- tibble(param_foc1 = "q_pc",
                    param_val1 = q_pc_vals) %>%
  expand_grid(tibble(param_foc2 = "z_pc",
                     param_val2 = z_pc_vals))

# RPV alone
# pdf("output/sensitivity_analysis_rpv_res_qz_pc.pdf")
# rpv_res_qz_pc <- qz_pc_in %>%
#   mutate(first_virus = "RPV") %>%
#   mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2, first_virus) param_fun(params_def2, param_foc1, param_val1, param_foc2, param_val2, first_virus, V0_b = 0))) %>%
#   unnest(cols = c(sim_out))
# dev.off()
# 
# # save simulation output (large)
# write_csv(rpv_res_qz_pc, "output/sensitivity_analysis_rpv_res_qz_pc.csv")
rpv_res_qz_pc <- read_csv("output/sensitivity_analysis_rpv_res_qz_pc.csv")

# Q_p values become negative in one simulation: remove
rpv_res_qz_pc2 <- rpv_res_qz_pc %>%
  filter(!(variable2 == "Q_p" & value < 0))

# figure
rpv_res_qz_pc2 %>%
  filter(variable2 == "Q_p") %>%
  ggplot(aes(x = param_val1, y = value, color = as.factor(param_val2))) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  geom_line() +
  facet_wrap(~ nutrient) +
  scale_color_viridis_d(name = "RPV P conc (z)") +
  labs(x = "RPV min P (q)", y = "Plant P concentration") +
  fig_theme


#### z_n -> biomass ####

# values for z
z_nc_vals3 <- 10^seq(-18, 0, length.out = 10)
z_nb_vals2 <- 10^seq(-18, 0, length.out = 10)

# data frame
z_n_in <- tibble(param_foc1 = "z_nb",
                 param_val1 = z_nb_vals2) %>%
  expand_grid(tibble(param_foc2 = "z_nc",
                     param_val2 = z_nc_vals3))

# set q values
params_q_n <- params_def2
params_q_n["q_nc"] <- 5.5e-3 # q_nc x 5
params_q_n["q_nb"] <- 1.1e-2 # q_nb x 10

# # PAV invasion
# pdf("output/sensitivity_analysis_pav_inv_z_n.pdf")
# pav_inv_z_n <- z_n_in %>%
#   mutate(first_virus = "RPV") %>%
#   mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2, first_virus) param_fun(params_q_n, param_foc1, param_val1, param_foc2, param_val2, first_virus, inv_time = (500-11-12)))) %>%
#   unnest(cols = c(sim_out))
# dev.off()
# 
# # save simulation output (large)
# write_csv(pav_inv_z_n, "output/sensitivity_analysis_pav_inv_z_n.csv")
pav_inv_z_n <- read_csv("output/sensitivity_analysis_pav_inv_z_n.csv")

# figure
pav_inv_z_n %>%
  filter(variable2 == "Q_n") %>%
  ggplot(aes(x = param_val1, y = value, color = as.factor(param_val2))) +
  geom_hline(yintercept = as.numeric(params_def1["Qmin_n"]), # plant's min
             color = "black", linetype = "dashed") +
  geom_hline(yintercept = 5.5e-3, # RPV's min
             color = "red", linetype = "dashed") +
  geom_hline(yintercept = 1.1e-2, # PAV's min
             color = "pink", linetype = "dashed") +
  geom_line() +
  facet_wrap(~ nutrient) +
  scale_x_log10() +
  scale_color_viridis_d(name = "RPV N conc (z)") +
  labs(x = "PAV N conc (z)", y = "Plant N conc") +
  fig_theme
# plant N conc settles on one of the three minimums
# depends on virus N conc

pav_inv_z_n %>%
  filter(variable2 == "H") %>%
  ggplot(aes(x = param_val1, y = value, color = as.factor(param_val2))) +
  geom_line() +
  facet_wrap(~ nutrient) +
  scale_x_log10() +
  scale_color_viridis_d(name = "RPV N conc (z)") +
  labs(x = "PAV N conc (z)", y = "Plant size (g)") +
  fig_theme
# when N is limiting, higher plant N conc from higher virus N conc (above)
# increases plant biomass


#### z_nc -> invasion ####

# values for z
z_nc_vals4 <- 10^seq(-18, 0, length.out = 100)
z_nb_vals3 <- 10^-4

# data frame
z_nc_in <- tibble(param_foc1 = "z_nb",
                 param_val1 = z_nb_vals3) %>%
  expand_grid(tibble(param_foc2 = "z_nc",
                     param_val2 = z_nc_vals4))

# set q values
params_q_n <- params_def2
params_q_n["q_nc"] <- 5.5e-3 # q_nc x 5
params_q_n["q_nb"] <- 1.1e-2 # q_nb x 10

# # PAV invasion
# pdf("output/sensitivity_analysis_pav_inv_z_nc.pdf")
# pav_inv_z_nc <- z_nc_in %>%
#   mutate(first_virus = "RPV") %>%
#   mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2, first_virus) param_fun(params_q_n, param_foc1, param_val1, param_foc2, param_val2, first_virus, inv_time = (500-11-12)))) %>%
#   unnest(cols = c(sim_out))
# dev.off()
# 
# # save simulation output (large)
# write_csv(pav_inv_z_nc, "output/sensitivity_analysis_pav_inv_z_nc.csv")
pav_inv_z_nc <- read_csv("output/sensitivity_analysis_pav_inv_z_nc.csv")

# scale virus values relative to initial
pav_inv_z_nc2 <- pav_inv_z_nc %>%
  mutate(value_scale = case_when(variable2 == "PAV_conc" ~ value / log10(V0 + 1),
                                 variable2 == "RPV_conc" ~ value / log10(V0 + 1),
                                 TRUE ~ NA_real_))

# figure
pav_inv_z_nc2 %>%
  filter(str_detect(variable2, "conc") == T) %>%
  ggplot(aes(x = param_val2, y = value_scale, color = variable2)) +
  geom_hline(yintercept = 1) +
  geom_line() +
  facet_wrap(~ nutrient) +
  scale_x_log10() +
  scale_color_viridis_d(name = "Virus") +
  labs(x = "RPV P conc (z)", y = "Relative virus concentration") +
  fig_theme

# plant biomass
pav_inv_z_nc2 %>%
  filter(variable2 == "H") %>%
  ggplot(aes(x = param_val2, y = value, color = nutrient, linetype = nutrient)) +
  geom_line() +
  scale_color_viridis_d(name = "Nutrient") +
  scale_linetype(name = "Nutrient") +
  fig_theme
# plants with high P reach biomass > 40 g


#### z_nc/z_pc -> invasion ####

# values for z
z_nc_vals4 <- 10^seq(-18, 0, length.out = 19)
z_pb_vals <- 10^seq(-18, 0, length.out = 19)

z_nc_vals4 <- 10^-18
z_pb_vals <- 10^-18

# data frame
z_np_in <- tibble(param_foc1 = "z_pb",
                 param_val1 = z_pb_vals) %>%
  expand_grid(tibble(param_foc2 = "z_nc",
                     param_val2 = z_nc_vals4))

# set q values
params_qz_np <- params_def2
params_qz_np["q_nc"] <- 5.5e-3 # q_nc x 5
params_qz_np["q_nb"] <- 1.1e-2 # q_nb x 10
params_qz_np["z_nb"] <- 10^-4
params_qz_np["q_pc"] <- 7.4e-4 # q_pc x 10
params_qz_np["q_pb"] <- 3.7e-4 # q_pb x 5
params_qz_np["z_pc"] <- 10^-4

# PAV invasion
pdf("output/sensitivity_analysis_pav_inv_z_np.pdf")
pav_inv_z_np <- z_np_in %>%
  mutate(first_virus = "RPV") %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2, first_virus) param_fun(params_qz_np, param_foc1, param_val1, param_foc2, param_val2, first_virus, inv_time = (500-11-12)))) %>%
  unnest(cols = c(sim_out))
dev.off()

# these simulations crash and the H values go very high - need to resolve