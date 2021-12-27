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

# run models for all 4 nutrients
sim_param_fun <- function(params, first_virus, plant_time, res_time, inv_time){
  
  # simulate plant growth
  plant_mod <- ode(c(R_n_low = R0_n_lo, R_p_low = R0_p_lo, Q_n_low = Q0_n, Q_p_low = Q0_p, H_low = H0, 
                     V_b_low = 0, V_c_low = 0,
                     R_n_n = R0_n_hi, R_p_n = R0_p_lo, Q_n_n = Q0_n, Q_p_n = Q0_p, H_n = H0, 
                     V_b_n = 0, V_c_n = 0,
                     R_n_p = R0_n_lo, R_p_p = R0_p_hi, Q_n_p = Q0_n, Q_p_p = Q0_p, H_p = H0, 
                     V_b_p = 0, V_c_p = 0,
                     R_n_np = R0_n_hi, R_p_np = R0_p_hi, Q_n_np = Q0_n, Q_p_np = Q0_p, H_np = H0, 
                     V_b_np = 0, V_c_np = 0),
                   times = seq(0, plant_time, length.out = 100), plant_virus_param_model, params) %>%
    as_tibble() %>%
    mutate(across(everything(), as.double)) 
  
  plant_init <- plant_mod %>%
    filter(time == plant_time)
  
  y0_first_virus <- c(R_n_low = pull(plant_init, R_n_low), 
                      R_p_low = pull(plant_init, R_p_low), 
                      Q_n_low = pull(plant_init, Q_n_low), 
                      Q_p_low = pull(plant_init, Q_p_low), 
                      H_low = pull(plant_init, H_low), 
                      V_b_low = if_else(first_virus == "PAV", V0_b, 0),
                      V_c_low = if_else(first_virus == "RPV", V0_c, 0), 
                      R_n_n = pull(plant_init, R_n_n), 
                      R_p_n = pull(plant_init, R_p_n), 
                      Q_n_n = pull(plant_init, Q_n_n), 
                      Q_p_n = pull(plant_init, Q_p_n), 
                      H_n = pull(plant_init, H_n), 
                      V_b_n = if_else(first_virus == "PAV", V0_b, 0),
                      V_c_n = if_else(first_virus == "RPV", V0_c, 0), 
                      R_n_p = pull(plant_init, R_n_p), 
                      R_p_p = pull(plant_init, R_p_p), 
                      Q_n_p = pull(plant_init, Q_n_p), 
                      Q_p_p = pull(plant_init, Q_p_p), 
                      H_p = pull(plant_init, H_p), 
                      V_b_p = if_else(first_virus == "PAV", V0_b, 0),
                      V_c_p = if_else(first_virus == "RPV", V0_c, 0), 
                      R_n_np = pull(plant_init, R_n_np), 
                      R_p_np = pull(plant_init, R_p_np), 
                      Q_n_np = pull(plant_init, Q_n_np), 
                      Q_p_np = pull(plant_init, Q_p_np), 
                      H_np = pull(plant_init, H_np), 
                      V_b_np = if_else(first_virus == "PAV", V0_b, 0),
                      V_c_np = if_else(first_virus == "RPV", V0_c, 0))

  # first virus model
  first_virus_mod <- ode(y0_first_virus, seq(0, res_time, length.out = 100),
                         plant_virus_param_model, params) %>%
    as_tibble() %>%
    mutate(across(everything(), as.double)) 
  
  first_virus_init <- first_virus_mod %>%
    filter(time == res_time)
  
  # virus initial conditions
  y0_second_virus <- c(R_n_low = pull(first_virus_init, R_n_low), 
                      R_p_low = pull(first_virus_init, R_p_low), 
                      Q_n_low = pull(first_virus_init, Q_n_low), 
                      Q_p_low = pull(first_virus_init, Q_p_low), 
                      H_low = pull(first_virus_init, H_low), 
                      V_b_low = if_else(first_virus == "PAV", 
                                        pull(first_virus_init, V_b_low), V0_b),
                      V_c_low = if_else(first_virus == "RPV", 
                                        pull(first_virus_init, V_c_low), V0_c), 
                      R_n_n = pull(first_virus_init, R_n_n), 
                      R_p_n = pull(first_virus_init, R_p_n), 
                      Q_n_n = pull(first_virus_init, Q_n_n), 
                      Q_p_n = pull(first_virus_init, Q_p_n), 
                      H_n = pull(first_virus_init, H_n), 
                      V_b_n = if_else(first_virus == "PAV", 
                                      pull(first_virus_init, V_b_n), V0_b),
                      V_c_n = if_else(first_virus == "RPV", 
                                      pull(first_virus_init, V_c_n), V0_c),
                      R_n_p = pull(first_virus_init, R_n_p), 
                      R_p_p = pull(first_virus_init, R_p_p), 
                      Q_n_p = pull(first_virus_init, Q_n_p), 
                      Q_p_p = pull(first_virus_init, Q_p_p), 
                      H_p = pull(first_virus_init, H_p), 
                      V_b_p = if_else(first_virus == "PAV", 
                                      pull(first_virus_init, V_b_p), V0_b),
                      V_c_p = if_else(first_virus == "RPV", 
                                      pull(first_virus_init, V_c_p), V0_c),
                      R_n_np = pull(first_virus_init, R_n_np), 
                      R_p_np = pull(first_virus_init, R_p_np), 
                      Q_n_np = pull(first_virus_init, Q_n_np), 
                      Q_p_np = pull(first_virus_init, Q_p_np), 
                      H_np = pull(first_virus_init, H_np), 
                      V_b_np = if_else(first_virus == "PAV", 
                                       pull(first_virus_init, V_b_np), V0_b),
                      V_c_np = if_else(first_virus == "RPV", 
                                       pull(first_virus_init, V_c_np), V0_c))
  
  # first virus model
  second_virus_mod <- ode(y0_second_virus, seq(0, inv_time, length.out = 100),
                          plant_virus_param_model, params) %>%
    as_tibble() %>%
    mutate(across(everything(), as.double))
  
  # combine simulations
  comb_mod <- plant_mod %>%
    filter(time != max(time)) %>%
    full_join(first_virus_mod %>%
                filter(time != max(time)) %>%
                mutate(time = time + plant_time)) %>%
    full_join(second_virus_mod %>%
                mutate(time = time + plant_time + res_time))
  
  return(comb_mod)
  
}

# change parameter value
param_fun <- function(params_in, param_foc, param_val, first_virus, max_days){
  
  # update parameter value
  params_in[param_foc] <- param_val
  
  # run model
  mod_out <- sim_param_fun(params = params_in, first_virus = first_virus, 
                           plant_time = 11, res_time = 12, inv_time = max_days-11-12)
  
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

# format output data from simulation
sim_dat_form <- function(dat_in){
  
  dat_out <- dat_in %>%
    pivot_longer(cols = -c(time),
                 names_to = "variable",
                 values_to = "value") %>%
    mutate(nutrient = case_when(str_ends(variable, "low") == T ~ "low",
                                str_ends(variable, "np") == T ~ "+N+P",
                                str_ends(variable, "n") == T ~ "+N",
                                str_ends(variable, "p") == T ~ "+P") %>%
             fct_relevel("low", "+N", "+P"),
           variable = case_when(str_starts(variable, "V_b") == T ~ "V_b",
                                str_starts(variable, "V_c") == T ~ "V_c",
                                str_starts(variable, "R_n") == T ~ "R_n",
                                str_starts(variable, "R_p") == T ~ "R_p",
                                str_starts(variable, "Q_n") == T ~ "Q_n",
                                str_starts(variable, "Q_p") == T ~ "Q_p",
                                str_starts(variable, "H") == T ~ "H"))
  
  return(dat_out)
  
}

# default parameters
params_def <- c(a_n_lo = a_n_lo,
                a_n_hi = a_n_hi,
                a_p_lo = a_p_lo,
                a_p_hi = a_p_hi,
                u_n = u_n,
                u_p = u_p,
                k_n = k_n,
                k_p = k_p,
                Qmin_n = Qmin_n,
                Qmin_p = Qmin_p,
                q_nb = q_nb,
                q_nc = q_nc,
                q_pb = q_pb,
                q_pc = q_pc,
                z_nb = z_nb,
                z_pb = z_pb,
                z_nc = z_nc,
                z_pc = z_pc,
                m = m,
                g = g,
                c_b = c_b,
                c_c = c_c,
                r_b = r_b,
                r_c = r_c)

# test function
V0_b <- V0_c <- V0_init
param_fun(params_def, "r_b", 0, "RPV", 1000)
param_fun(params_def, "r_c", 0.08, "PAV", 100)
sim_param_fun(params_def, "RPV", 11, 12, 100-11-12) %>%
  filter(time %in% c(11, (11+12)))


#### initial simulation settings ####

max_days <- 100

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
