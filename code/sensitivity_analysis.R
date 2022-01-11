## Goal: sensitivity analysis of virus parameters


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


#### model functions ####

# change parameter value
param_fun <- function(params_in, param_foc1, param_val1, param_foc2, param_val2, 
                      first_virus, V0_b = V0, V0_c = V0, inits_in = init_def2,
                      plant_time = 11, res_time = 12, inv_time = 100-11-12){
  
  # update parameter value
  params_in[param_foc1] <- param_val1
  params_in[param_foc2] <- param_val2
  
  # update resource supply
  if(param_foc1 == "a_n_lo"){
    # inits_in["R_n_low"] <- 5.6e-5 * 7 + param_val1 * 3
    # inits_in["R_n_p"] <- 5.6e-5 * 7 + param_val1 * 3
    inits_in["R_n_low"] <- param_val1 * 10
    inits_in["R_n_p"] <- param_val1 * 10
  }
  if(param_foc2 == "a_p_lo"){
    # inits_in["R_p_low"] <- 8.2e-6 * 7 + param_val2 * 3
    # inits_in["R_p_n"] <- 8.2e-6 * 7 + param_val2 * 3
    inits_in["R_p_low"] <- param_val2 * 10
    inits_in["R_p_n"] <- param_val2 * 10
  }
  
  # run model
  mod_out <- virus2_model_sim(params = params_in, first_virus = first_virus, 
                              V0_b = V0_b, V0_c = V0_c,
                              plant_time = plant_time, res_time = res_time, inv_time = inv_time, 
                              inits = inits_in)
  
  # edit output
  mod_out2 <- mod_out %>%
    filter(time == max(time)) %>%
    mutate(VbH_low = V_b_low * H_low,
           VbH_n = V_b_n * H_n,
           VbH_p = V_b_p * H_p,
           VbH_np = V_b_np * H_np,
           VcH_low = V_c_low * H_low,
           VcH_n = V_c_n * H_n,
           VcH_p = V_c_p * H_p,
           VcH_np = V_c_np * H_np) %>%
    pivot_longer(cols = everything(),
                 names_to = "variable",
                 values_to = "value") %>%
    mutate(nutrient = case_when(str_ends(variable, "low") == T ~ "low",
                                str_ends(variable, "np") == T ~ "+N+P",
                                str_ends(variable, "n") == T ~ "+N",
                                str_ends(variable, "p") == T ~ "+P") %>%
             fct_relevel("low", "+N", "+P"),
           variable2 = case_when(str_starts(variable, "R_n") == T ~ "R_n",
                                 str_starts(variable, "R_p") == T ~ "R_p",
                                 str_starts(variable, "Q_n") == T ~ "Q_n",
                                 str_starts(variable, "Q_p") == T ~ "Q_p",
                                 str_starts(variable, "H") == T ~ "H",
                                 str_starts(variable, "V_b") == T ~ "PAV_conc",
                                 str_starts(variable, "V_c") == T ~ "RPV_conc",
                                 str_starts(variable, "VbH") == T ~ "PAV_pop",
                                 str_starts(variable, "VcH") == T ~ "RPV_pop"),
           virus_abund = if_else(str_starts(variable, "V") == T, 
                                 log10(value + 1), NA_real_),
           abund_type = if_else(str_starts(variable, "V") == T, 
                                str_sub(variable2, 5, 7), NA_character_))
  
  # return output
  return(mod_out2)

}

# test function
virus2_model_sim(params_def2, "RPV", V0_b = V0, V0_c = V0,
                 plant_time = 11, res_time = 12, 
                 inv_time = 100-11-12) %>%
  virus2_model_format(params_def2) %>%
  filter(time == max(time) & variable2 %in% c("PAV_conc", "PAV_pop")) %>%
  mutate(virus_conc = log10(value + 1))

param_fun(params_def2, "r_b", 0.0961, "r_c", 0.221, "RPV") %>%
  filter(variable2 %in% c("PAV_conc", "PAV_pop"))

param_fun(params_def2, "a_n_lo", 1.1e-10, "r_c", 0.221, "RPV") %>%
  filter(variable2 %in% c("PAV_conc", "PAV_pop"))


#### resource supply rates ####

# values for low nutrient supply
a_n_vals <- 10^(-13:-4)
a_p_vals <- 10^(-13:-4)

# data frame
a_in <- tibble(param_foc1 = "a_n_lo",
               param_val1 = a_n_vals) %>%
  expand_grid(tibble(param_foc2 = "a_p_lo",
                     param_val2 = a_p_vals))

# PAV invades RPV
pav_inv_a <- a_in %>%
  mutate(first_virus = "RPV") %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2, first_virus) param_fun(params_def2, param_foc1, param_val1, param_foc2, param_val2, first_virus))) %>%
  unnest(cols = c(sim_out))

# figure
pav_inv_a %>%
  filter(variable2 == "PAV_conc") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = virus_abund)) +
  geom_point(size = 10) +
  facet_wrap(~ nutrient) +
  scale_color_viridis_c() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "N supply", y = "P supply")

pav_inv_a %>%
  filter(variable2 == "PAV_pop") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = virus_abund)) +
  geom_point(size = 10) +
  facet_wrap(~ nutrient) +
  scale_color_viridis_c() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "N supply", y = "P supply")

pav_inv_a %>%
  filter(variable2 == "H") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = value)) +
  geom_point(size = 10) +
  facet_wrap(~ nutrient) +
  scale_color_viridis_c() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "N supply", y = "P supply")

# plant dynamics with very low nutrients
n <- 1

params_a_low <- params_def2
params_a_low["a_n_low"] <- a_n_vals[n]
params_a_low["a_p_low"] <- a_p_vals[n]

inits_a_low <- init_def2
inits_a_low["R_n_low"] <- a_n_vals[n] * 10
inits_a_low["R_n_p"] <- a_n_vals[n] * 10
inits_a_low["R_p_low"] <- a_p_vals[n] * 10
inits_a_low["R_p_n"] <- a_p_vals[n] * 10

plant_a_low <- virus2_model_sim(params_a_low, "RPV", V0_b = V0, V0_c = V0, 
                                plant_time = 11, res_time = 12, 
                                inv_time = 100-11-12,
                                inits = inits_a_low) %>%
  virus2_model_format(params_a_low)
# can ignore warning message about unknown levels in 'f'
# occurs if one of the nutrients is never limiting
# function still works

plant_fig_fun(plant_a_low, params_a_low, -2e-3, -5e-4)
# Q is only drawn down so fast, giving viruses some nutrients


#### later invasion ####

# use default parameters
param_fun(params_def2, "r_b", 0.0961, "r_c", 0.221, "RPV",
          plant_time = 11, res_time = 100-11, inv_time = 100) %>%
  filter(variable2 %in% c("PAV_conc", "PAV_pop"))
# abundances are comparable to very low nutrients
# viruses can always invade when Q is at Qmin?


#### minimum nutrient concentrations ####

# values for low nutrient supply
q_nb_vals <- 10^(-3:6)
q_nc_vals <- 10^(-3:6)

# data frame
q_n_in <- tibble(param_foc1 = "q_nb",
                 param_val1 = q_nb_vals) %>%
  expand_grid(tibble(param_foc2 = "q_nc",
                     param_val2 = q_nc_vals))

# PAV invades RPV
pav_inv_q_n <- q_n_in %>%
  mutate(first_virus = "RPV") %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2, first_virus) param_fun(params_def2, param_foc1, param_val1, param_foc2, param_val2, first_virus))) %>%
  unnest(cols = c(sim_out))

# figure
pav_inv_q_n %>%
  filter(variable2 == "PAV_conc") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = virus_abund)) +
  geom_point(size = 10) +
  facet_wrap(~ nutrient) +
  scale_color_viridis_c() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "PAV min N", y = "RPV min N")

pav_inv_q_n %>%
  filter(variable2 == "RPV_conc") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = virus_abund)) +
  geom_point(size = 10) +
  facet_wrap(~ nutrient) +
  scale_color_viridis_c() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "PAV min N", y = "RPV min N")
# when q_ values are too high, viruses can't invade
# this is about the virus-plant interaction, not virus-virus


#### minimum nutrient concentration/nutrient content of virus ####

#### start here ####
# below values make simulation crash and output data can't be processed
# adjust values to be less extreme

# values for low nutrient supply
z_nc_vals <- 10^seq(-18, 0, length.out = 10)
q_nc_vals2 <- 10^(-11:-2)

# data frame
qz_nc_in <- tibble(param_foc1 = "z_nc",
                   param_val1 = z_nc_vals) %>%
  expand_grid(tibble(param_foc2 = "q_nc",
                     param_val2 = q_nc_vals2))

# PAV invades RPV
pav_inv_qz_nc <- qz_nc_in %>%
  mutate(first_virus = "RPV") %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2, first_virus) param_fun(params_def2, param_foc1, param_val1, param_foc2, param_val2, first_virus))) %>%
  unnest(cols = c(sim_out))

# figure
pav_inv_qz_nc %>%
  filter(variable2 == "PAV_conc") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = virus_abund)) +
  geom_point(size = 10) +
  facet_wrap(~ nutrient) +
  scale_color_viridis_c() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "PAV min N", y = "RPV min N")



#### older code ####

#### initial simulation settings ####

n <- 10

param_df <- tibble(param_foc = rep(c("q_nb", "q_pb", "q_nc", "q_pc", 
                                     "z_nb", "z_pb", "z_nc", "z_pc", 
                                     "c_b", "c_c", 
                                     "r_b", "r_c"), each = n), 
                   param_val = c(rep(seq(1e-6, 1e-2, length.out = n), 4),
                                 rep(seq(1e-20, 1e-17, length.out = n), 4),
                                 rep(seq(1e-4, 1e-2, length.out = n), 2),
                                 rep(seq(1e-3, 1, length.out = n), 2)))

#### PAV initial simulations ####

# alone
V0_b <- V0_init
V0_c <- 0
pav_only_params <- param_df %>%
  mutate(first_virus = "RPV") %>%
  mutate(pav_conc = pmap(., function(param_foc, param_val, first_virus) param_fun(params_def, param_foc, param_val, first_virus, max_days))) %>%
  unnest(cols = c(pav_conc)) %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P"))

# invades
V0_b <- V0_c <- V0_init
pav_inv_params <- param_df %>%
  mutate(first_virus = "RPV") %>%
  mutate(pav_conc = pmap(., function(param_foc, param_val, first_virus) param_fun(params_def, param_foc, param_val, first_virus, max_days))) %>%
  unnest(cols = c(pav_conc)) %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P"))

# figures
ggplot(pav_only_params, 
       aes(x = param_val, y = virus_conc, color = nutrient, linetype = nutrient)) +
  geom_line() +
  facet_rep_wrap(~ param_foc, scales = "free_x") +
  fig_theme

ggplot(pav_inv_params, 
       aes(x = param_val, y = virus_conc, color = nutrient, linetype = nutrient)) +
  geom_line() +
  facet_rep_wrap(~ param_foc, scales = "free_x") +
  fig_theme
# when RPV has a lower q for the non-limiting nutrient, PAV concentration is much lower (q_nc)


#### RPV initial simulations ####

# alone
V0_c <- V0_init
V0_b <- 0
rpv_only_params <- param_df %>%
  mutate(first_virus = "PAV") %>%
  mutate(pav_conc = pmap(., function(param_foc, param_val, first_virus) param_fun(params_def, param_foc, param_val, first_virus, max_days))) %>%
  unnest(cols = c(pav_conc)) %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P"))

# invades
V0_b <- V0_c <- V0_init
rpv_inv_params <- param_df %>%
  mutate(first_virus = "PAV") %>%
  mutate(pav_conc = pmap(., function(param_foc, param_val, first_virus) param_fun(params_def, param_foc, param_val, first_virus, max_days))) %>%
  unnest(cols = c(pav_conc)) %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P"))

# figures
ggplot(rpv_only_params, 
       aes(x = param_val, y = virus_conc, color = nutrient, linetype = nutrient)) +
  geom_line() +
  facet_rep_wrap(~ param_foc, scales = "free_x") +
  fig_theme

ggplot(rpv_inv_params, 
       aes(x = param_val, y = virus_conc, color = nutrient, linetype = nutrient)) +
  geom_line() +
  facet_rep_wrap(~ param_foc, scales = "free_x") +
  fig_theme
# don't see the same q effect here, but range may be too limited


#### q simulation settings ####

max_days <- 100

n <- 7

param_q_pav <- tibble(param_foc = rep(c("q_nc", "q_pc"), each = n),
                      param_val = c(1.1e-9, 1.1e-8, 1.1e-7, 1.1e-6, 1.1e-5, 1.1e-4, 1.1e-3,
                                    7.4e-17, 7.4e-15, 7.4e-13, 7.4e-11, 7.4e-9, 7.4e-7, 7.4e-5))

param_q_rpv <- tibble(param_foc = rep(c("q_nb", "q_pb"), each = n), 
                      param_val = c(1.1e-15, 1.1e-13, 1.1e-11, 1.1e-9, 1.1e-7, 1.1e-5, 1.1e-3,
                                    7.4e-17, 7.4e-15, 7.4e-13, 7.4e-11, 7.4e-9, 7.4e-7, 7.4e-5))


#### PAV q simulations ####

V0_b <- V0_c <- V0_init
pav_q_params <- param_q_pav %>%
  mutate(first_virus = "RPV") %>%
  mutate(pav_conc = pmap(., function(param_foc, param_val, first_virus) param_fun(params_def, param_foc, param_val, first_virus, max_days))) %>%
  unnest(cols = c(pav_conc)) %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P"))

ggplot(pav_q_params, 
       aes(x = param_val, y = virus_conc, color = nutrient, linetype = nutrient)) +
  geom_line() +
  facet_rep_wrap(~ param_foc, scales = "free_x") +
  scale_x_log10() +
  fig_theme

# why is PAV sensitive to q_nc?
params_qnc <- params_def

params_qnc["q_nc"] <- 1.1e-6

pav_inv_q_nc <- sim_param_fun(params_qnc, "RPV", 11, 12, max_days-11-12) %>%
  sim_dat_form()

ggplot(pav_inv_q_nc, aes(x = time, y = value, color = nutrient, linetype = nutrient)) +
  geom_line() +
  facet_rep_wrap(~ variable, scales = "free_y") +
  fig_theme
# RPV explodes and kills the plant, PAV was still able to invade
  

#### RPV q simulations ####

V0_b <- V0_c <- V0_init
rpv_q_params <- param_q_rpv %>%
  mutate(first_virus = "PAV") %>%
  mutate(pav_conc = pmap(., function(param_foc, param_val, first_virus) param_fun(params_def, param_foc, param_val, first_virus, max_days))) %>%
  unnest(cols = c(pav_conc)) %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P"))

ggplot(rpv_q_params, 
       aes(x = param_val, y = virus_conc, color = nutrient, linetype = nutrient)) +
  geom_line() +
  facet_rep_wrap(~ param_foc, scales = "free_x") +
  scale_x_log10() +
  fig_theme


#### r simulation settings ####

max_days <- 100

n <- 4

param_r_pav <- tibble(param_foc = rep("r_c", n),
                      param_val = c(1e-2, 1e-1, 1, 10))

param_r_rpv <- tibble(param_foc = rep("r_b", n),
                      param_val = c(1e-2, 1e-1, 1, 10))


#### PAV r simulations ####

V0_b <- V0_c <- V0_init
pav_r_params <- param_r_pav %>%
  mutate(first_virus = "RPV") %>%
  mutate(pav_conc = pmap(., function(param_foc, param_val, first_virus) param_fun(params_def, param_foc, param_val, first_virus, max_days))) %>%
  unnest(cols = c(pav_conc)) %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P"))

ggplot(pav_r_params, 
       aes(x = param_val, y = virus_conc, color = nutrient, linetype = nutrient)) +
  geom_line() +
  facet_rep_wrap(~ param_foc, scales = "free_x") +
  scale_x_log10() +
  fig_theme

# time series
params_rc <- params_def

params_rc["r_c"] <- 10

pav_inv_rc <- sim_param_fun(params_rc, "RPV", 11, 12, max_days-11-12) %>%
  sim_dat_form()

ggplot(pav_inv_rc, aes(x = time, y = value, color = nutrient, linetype = nutrient)) +
  geom_line() +
  facet_rep_wrap(~ variable, scales = "free_y") +
  fig_theme
# plant is getting smaller, but not dead
# PAV can't invade (decreasing almost immediately)
# RPV values are probably unrealistic

ggplot(filter(pav_inv_rc, time <= 25), 
       aes(x = time, y = value, color = nutrient, linetype = nutrient)) +
  geom_line() +
  facet_rep_wrap(~ variable, scales = "free_y") +
  fig_theme

# look at internal nutrients without viruses
V0_b <- V0_c <- 0
no_inv_rc <- sim_param_fun(params_rc, "RPV", 11, 12, max_days-11-12) %>%
  sim_dat_form()
ggplot(filter(no_inv_rc, time <= 25), aes(x = time, y = value, color = nutrient, linetype = nutrient)) +
  geom_line() +
  facet_rep_wrap(~ variable, scales = "free_y") +
  fig_theme
# the fast virus drags down internal nutrients
# probably similar mechanism to the q simulation above with RPV


#### RPV r simulations ####

V0_b <- V0_c <- V0_init
rpv_r_params <- param_r_rpv %>%
  mutate(first_virus = "PAV") %>%
  mutate(pav_conc = pmap(., function(param_foc, param_val, first_virus) param_fun(params_def, param_foc, param_val, first_virus, max_days))) %>%
  unnest(cols = c(pav_conc)) %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P"))

ggplot(rpv_r_params, 
       aes(x = param_val, y = virus_conc, color = nutrient, linetype = nutrient)) +
  geom_line() +
  facet_rep_wrap(~ param_foc, scales = "free_x") +
  scale_x_log10() +
  fig_theme

# time series
params_rb <- params_def

params_rb["r_b"] <- 10

rpv_inv_rb <- sim_param_fun(params_rb, "PAV", 11, 12, max_days-11-12) %>%
  sim_dat_form()

ggplot(rpv_inv_rb, aes(x = time, y = value, color = nutrient, linetype = nutrient)) +
  geom_line() +
  facet_rep_wrap(~ variable, scales = "free_y") +
  fig_theme
# same as PAV r simulation


#### resource supply rate settings ####

max_days <- 100

n <- 6

param_a_lo <- tibble(param_foc = rep(c("a_n_lo", "a_p_lo"), each = n),
                     param_val = c(1.1e-11, 1.1e-10, 1.1e-9, 1.1e-8, 1.1e-7, 1.1e-6,
                                   1.6e-12, 1.6e-11, 1.6e-10, 1.6e-9, 1.6e-8, 1.6e-7))


#### PAV a simulations ####

V0_b <- V0_c <- V0_init
pav_a_params <- param_a_lo %>%
  mutate(first_virus = "RPV") %>%
  mutate(pav_conc = pmap(., function(param_foc, param_val, first_virus) param_fun(params_def, param_foc, param_val, first_virus, max_days))) %>%
  unnest(cols = c(pav_conc)) %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P"))

ggplot(pav_a_params, 
       aes(x = param_val, y = virus_conc, color = nutrient, linetype = nutrient)) +
  geom_line() +
  facet_rep_grid(nutrient ~ param_foc, scales = "free") +
  scale_x_log10() +
  fig_theme

# time series
params_alo <- params_def

params_alo["a_n_lo"] <- 1.1e-13

pav_inv_alo <- sim_param_fun(params_alo, "RPV", 12, 11, max_days-11-12) %>%
  sim_dat_form()

ggplot(pav_inv_alo, aes(x = time, y = value, color = nutrient, linetype = nutrient)) +
  geom_line() +
  facet_rep_wrap(~ variable, scales = "free_y") +
  fig_theme

# compare to regular supply rates
pav_inv <- sim_param_fun(params_def, "RPV", 12, 11, max_days-11-12) %>%
  sim_dat_form()

ggplot(pav_inv, aes(x = time, y = value, color = nutrient, linetype = nutrient)) +
  geom_line() +
  facet_rep_wrap(~ variable, scales = "free_y") +
  fig_theme

#### Start here: double check simulations in ode_simulations ####
# the above figure doesn't match the manuscript figure and it should
# error in one of the models?
# can initial resource amount be changed to be a function of low/high?


#### RPV a simulations ####

V0_b <- V0_c <- V0_init
rpv_a_params <- param_a_lo %>%
  mutate(first_virus = "PAV") %>%
  mutate(pav_conc = pmap(., function(param_foc, param_val, first_virus) param_fun(params_def, param_foc, param_val, first_virus, max_days))) %>%
  unnest(cols = c(pav_conc)) %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P"))

ggplot(rpv_a_params, 
       aes(x = param_val, y = virus_conc, color = nutrient, linetype = nutrient)) +
  geom_line() +
  facet_rep_grid(nutrient ~ param_foc, scales = "free") +
  scale_x_log10() +
  fig_theme
