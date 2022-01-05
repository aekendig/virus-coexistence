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
    inits_in["R_n_low"] <- 5.6e-5 * 7 + param_val1 * 3
    inits_in["R_n_p"] <- 5.6e-5 * 7 + param_val1 * 3
  }
  if(param_foc2 == "a_p_lo"){
    inits_in["R_p_low"] <- 8.2e-6 * 7 + param_val2 * 3
    inits_in["R_p_n"] <- 8.2e-6 * 7 + param_val2 * 3
  }
  
  # run model
  mod_out <- virus2_model_sim(params = params_in, first_virus = first_virus, 
                              V0_b = V0_b, V0_c = V0_c,
                              plant_time = plant_time, res_time = res_time, inv_time = inv_time, 
                              inits = inits_in)
  
  # edit output
  mod_out2 <- mod_out %>%
    filter(time == max(time)) %>%
    pivot_longer(cols = everything(),
                 names_to = "variable",
                 values_to = "value") %>%
    mutate(nutrient = case_when(str_ends(variable, "low") == T ~ "low",
                                str_ends(variable, "np") == T ~ "+N+P",
                                str_ends(variable, "n") == T ~ "+N",
                                str_ends(variable, "p") == T ~ "+P") %>%
             fct_relevel("low", "+N", "+P"),
           variable = case_when(str_starts(variable, "V_b") == T ~ "V_b",
                                str_starts(variable, "V_c") == T ~ "V_c",
                                TRUE ~ variable),
           virus_conc = log10(value + 1))
  
  # output: invading virus titer
  if(first_virus == "PAV"){
    output <- mod_out2 %>%
      filter(variable == "V_c") %>%
      select(nutrient, virus_conc)
  }else{
    output <- mod_out2 %>%
      filter(variable == "V_b") %>%
      select(nutrient, virus_conc)
  }
  
  # return output
  return(output)

}

# test function
virus2_model_sim(params_def2, "RPV", V0_b = V0, V0_c = V0,
                 plant_time = 11, res_time = 12, 
                 inv_time = 100-11-12) %>%
  virus2_model_format(params_def2) %>%
  filter(time == max(time) & variable2 == "PAV") %>%
  mutate(virus_conc = log10(value + 1))
param_fun(params_def2, "r_b", 0.0961, "r_c", 0.221, "RPV")
param_fun(params_def2, "a_n_lo", 1.1e-10, "r_c", 0.221, "RPV")


#### resource supply rates ####

# values for low 
a_n_vals <- 1.1*10^(-15:-6)
a_p_vals <- 1.6*10^(-16:-7)

# data frame
a_in <- tibble(param_foc1 = "a_n_lo",
               param_val1 = a_n_vals) %>%
  expand_grid(tibble(param_foc2 = "a_p_lo",
                     param_val2 = a_p_vals))

# PAV invades RPV
pav_inv_a <- a_in %>%
  mutate(first_virus = "RPV") %>%
  mutate(pav_conc = pmap(., function(param_foc1, param_val1, param_foc2, param_val2, first_virus) param_fun(params_def2, param_foc1, param_val1, param_foc2, param_val2, first_virus))) %>%
  unnest(cols = c(pav_conc)) %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P"))

# figure
pav_inv_a %>%
  filter(nutrient == "low") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = virus_conc)) +
  geom_point(size = 10) +
  scale_x_log10() +
  scale_y_log10()

# RPV invades PAV
rpv_inv_a <- a_in %>%
  mutate(first_virus = "PAV") %>%
  mutate(rpv_conc = pmap(., function(param_foc1, param_val1, param_foc2, param_val2, first_virus) param_fun(params_def2, param_foc1, param_val1, param_foc2, param_val2, first_virus))) %>%
  unnest(cols = c(rpv_conc)) %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P"))

# figure
rpv_inv_a %>%
  filter(nutrient == "low") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = virus_conc)) +
  geom_point(size = 10) +
  scale_x_log10() +
  scale_y_log10()

# minimum titer
min(pav_inv_a$virus_conc)
min(rpv_inv_a$virus_conc)

#### start here ####
# both viruses plateau at a minimum titer as resources are reduced
# the plant could be dead in a lot of these cases
# total titer is probably more relevant than concentration
# this was one of Eric's comments


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
