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


#### model function ####

# run model for all 4 nutrients
sim_param_fun <- function(params, first_virus, plant_time, res_time, inv_time){
  
  # simulate plant growth
  plant_init <- ode(c(R_n_low = R0_n_lo, R_p_low = R0_p_lo, Q_n_low = Q0_n, Q_p_low = Q0_p, H_low = H0, 
                      V_b_low = 0, V_c_low = 0,
                      R_n_n = R0_n_hi, R_p_n = R0_p_lo, Q_n_n = Q0_n, Q_p_n = Q0_p, H_n = H0, 
                      V_b_n = 0, V_c_n = 0,
                      R_n_p = R0_n_lo, R_p_p = R0_p_hi, Q_n_p = Q0_n, Q_p_p = Q0_p, H_p = H0, 
                      V_b_p = 0, V_c_p = 0,
                      R_n_np = R0_n_hi, R_p_np = R0_p_hi, Q_n_np = Q0_n, Q_p_np = Q0_p, H_np = H0, 
                      V_b_np = 0, V_c_np = 0),
                    times = seq(0, plant_time, length.out = 100), plant_virus_param_model, params) %>%
    as_tibble() %>%
    mutate(across(everything(), as.double)) %>%
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
  first_virus_init <- ode(y0_first_virus, seq(0, res_time, length.out = 100),
                         plant_virus_param_model, params) %>%
    as_tibble() %>%
    mutate(across(everything(), as.double)) %>%
    filter(time == res_time)
  
  # virus initial conditions
  y0_second_virus <- c(R_n_low = pull(plant_init, R_n_low), 
                      R_p_low = pull(plant_init, R_p_low), 
                      Q_n_low = pull(plant_init, Q_n_low), 
                      Q_p_low = pull(plant_init, Q_p_low), 
                      H_low = pull(plant_init, H_low), 
                      V_b_low = if_else(first_virus == "PAV", pull(first_virus_init, V_b_low), V0_b),
                      V_c_low = if_else(first_virus == "RPV", pull(first_virus_init, V_c_low), V0_c), 
                      R_n_n = pull(plant_init, R_n_n), 
                      R_p_n = pull(plant_init, R_p_n), 
                      Q_n_n = pull(plant_init, Q_n_n), 
                      Q_p_n = pull(plant_init, Q_p_n), 
                      H_n = pull(plant_init, H_n), 
                      V_b_n = if_else(first_virus == "PAV", pull(first_virus_init, V_b_n), V0_b),
                      V_c_n = if_else(first_virus == "RPV", pull(first_virus_init, V_c_n), V0_c),
                      R_n_p = pull(plant_init, R_n_p), 
                      R_p_p = pull(plant_init, R_p_p), 
                      Q_n_p = pull(plant_init, Q_n_p), 
                      Q_p_p = pull(plant_init, Q_p_p), 
                      H_p = pull(plant_init, H_p), 
                      V_b_p = if_else(first_virus == "PAV", pull(first_virus_init, V_b_p), V0_b),
                      V_c_p = if_else(first_virus == "RPV", pull(first_virus_init, V_c_p), V0_c),
                      R_n_np = pull(plant_init, R_n_np), 
                      R_p_np = pull(plant_init, R_p_np), 
                      Q_n_np = pull(plant_init, Q_n_np), 
                      Q_p_np = pull(plant_init, Q_p_np), 
                      H_np = pull(plant_init, H_np), 
                      V_b_np = if_else(first_virus == "PAV", pull(first_virus_init, V_b_np), V0_b),
                      V_c_np = if_else(first_virus == "RPV", pull(first_virus_init, V_c_np), V0_c))
  
  # first virus model
  second_virus_mod <- ode(y0_second_virus, seq(0, inv_time, length.out = 100),
                          plant_virus_param_model, params) %>%
    as_tibble() %>%
    mutate(across(everything(), as.double)) %>%
    filter(time == inv_time) %>%
    pivot_longer(cols = everything(),
                 names_to = "variable",
                 values_to = "value") %>%
    mutate(nutrient = case_when(str_ends(variable, "low") == T ~ "low",
                                str_ends(variable, "np") == T ~ "+N+P",
                                str_ends(variable, "n") == T ~ "+N",
                                str_ends(variable, "p") == T ~ "+P") %>%
             fct_relevel("low", "+N", "+P"),
           variable = case_when(str_starts(V_b) == T ~ "V_b",
                                str_starts(V_c) == T ~ "V_c",
                                TRUE ~ variable),
           virus_conc = log10(value + 1))
  
  # output: invading virus titer
  if(first_virus == "PAV"){
    output <- second_virus_mod %>%
      filter(variable == "V_c") %>%
      select(nutrient, virus_conc)
  }else{
    output <- second_virus_mod %>%
      filter(variable == "V_b") %>%
      select(nutrient, virus_conc)
  }
  
  return(output)
  
}

#### start here ####
# writing function so that parameters are input and global parameter values aren't reset

param_fun <- function(params, param_foc, param_val, first_virus){
  
  # update parameter value
  params[param_foc] <- param_val
  
  # extract parameter values
  q_nb <<- as.numeric(params["q_nb"])
  q_pb <<- as.numeric(params["q_pb"])
  q_nc <<- as.numeric(params["q_nc"])
  q_pc <<- as.numeric(params["q_pc"])
  z_nb <<- as.numeric(params["z_nb"])
  z_pb <<- as.numeric(params["z_pb"])
  z_nc <<- as.numeric(params["z_nc"])
  z_pc <<- as.numeric(params["z_pc"])
  c_b <<- as.numeric(params["c_b"])
  c_c <<- as.numeric(params["c_c"])
  r_b <<- as.numeric(params["r_b"])
  r_c <<- as.numeric(params["r_c"])
  
  # run model
  mod_low <- sim_fun("low", "low", "low", first_virus = first_virus, 
                 plant_time = 12, res_time = 11, inv_time = max_days-11-12)
  mod_n <- sim_fun("high", "low", "+N", first_virus = first_virus, 
                 plant_time = 12, res_time = 11, inv_time = max_days-11-12)
  mod_p <- sim_fun("low", "high", "+P", first_virus = first_virus, 
                 plant_time = 12, res_time = 11, inv_time = max_days-11-12)
  mod_np <- sim_fun("high", "high", "+N+P", first_virus = first_virus, 
                 plant_time = 12, res_time = 11, inv_time = max_days-11-12)
  
  # combine models
  mod <- sim_dat_fun(mod_low, mod_n, mod_p, mod_np)
  
  # output: invading virus titer
  if(first_virus == "PAV"){
    output <- mod %>%
      select(nutrient, RPV_log10)
  }else{
    output <- mod %>%
      select(nutrient, PAV_log10)
  }
  
  # return output
  return(output)

}


# default parameters (only run once)
params_def <- list("q_nb" = q_nb,
                   "q_pb" = q_pb,
                   "q_nc" = q_nc,
                   "q_pc" = q_pc,
                   "z_nb" = z_nb,
                   "z_pb" = z_pb,
                   "z_nc" = z_nc,
                   "z_pc" = z_pc,
                   "c_b" = c_b,
                   "c_c" = c_c,
                   "r_b" = r_b,
                   "r_c" = r_c)

# test function
pav_inv <- param_fun(params_def, "r_b", 0, "RPV")
rpv_inv <- param_fun(params_def, "r_c", 0.08, "PAV")


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
  mutate(pav_conc = pmap(., function(param_foc, param_val, first_virus) param_fun(params_def, param_foc, param_val, first_virus))) %>%
  unnest(cols = c(pav_conc)) %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P"))

# invades
V0_b <- V0_c <- V0_init
pav_inv_params <- param_df %>%
  mutate(first_virus = "RPV") %>%
  mutate(pav_conc = pmap(., function(param_foc, param_val, first_virus) param_fun(params_def, param_foc, param_val, first_virus))) %>%
  unnest(cols = c(pav_conc)) %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P"))

# figures
ggplot(pav_only_params, 
       aes(x = param_val, y = PAV_log10, color = nutrient, linetype = nutrient)) +
  geom_line() +
  facet_rep_wrap(~ param_foc, scales = "free_x") +
  fig_theme

ggplot(pav_inv_params, 
       aes(x = param_val, y = PAV_log10, color = nutrient, linetype = nutrient)) +
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
  mutate(pav_conc = pmap(., function(param_foc, param_val, first_virus) param_fun(params_def, param_foc, param_val, first_virus))) %>%
  unnest(cols = c(pav_conc)) %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P"))

# invades
V0_b <- V0_c <- V0_init
rpv_inv_params <- param_df %>%
  mutate(first_virus = "PAV") %>%
  mutate(pav_conc = pmap(., function(param_foc, param_val, first_virus) param_fun(params_def, param_foc, param_val, first_virus))) %>%
  unnest(cols = c(pav_conc)) %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P"))

# figures
ggplot(rpv_only_params, 
       aes(x = param_val, y = RPV_log10, color = nutrient, linetype = nutrient)) +
  geom_line() +
  facet_rep_wrap(~ param_foc, scales = "free_x") +
  fig_theme

ggplot(rpv_inv_params, 
       aes(x = param_val, y = RPV_log10, color = nutrient, linetype = nutrient)) +
  geom_line() +
  facet_rep_wrap(~ param_foc, scales = "free_x") +
  fig_theme
# don't see the same q effect here, but range may be too limited


#### q simulation settings ####

max_days <- 100

n <- 7

param_df_pav <- tibble(param_foc = rep(c("q_nc", "q_pc"), each = n),
                       param_val = c(1.1e-9, 1.1e-8, 1.1e-7, 1.1e-6, 1.1e-5, 1.1e-4, 1.1e-3,
                                     7.4e-17, 7.4e-15, 7.4e-13, 7.4e-11, 7.4e-9, 7.4e-7, 7.4e-5))

param_df_rpv <- tibble(param_foc = rep(c("q_nb", "q_pb"), each = n), 
                       param_val = c(1.1e-15, 1.1e-13, 1.1e-11, 1.1e-9, 1.1e-7, 1.1e-5, 1.1e-3,
                                     7.4e-17, 7.4e-15, 7.4e-13, 7.4e-11, 7.4e-9, 7.4e-7, 7.4e-5))


#### PAV q simulations ####

V0_b <- V0_c <- V0_init
q_nb <- 1.1e-3
q_nc <- 7.4e-5
pav_q_params <- param_df_pav %>%
  mutate(first_virus = "RPV") %>%
  mutate(pav_conc = pmap(., function(param_foc, param_val, first_virus) param_fun(params_def, param_foc, param_val, first_virus))) %>%
  unnest(cols = c(pav_conc)) %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P"))

ggplot(pav_q_params, 
       aes(x = param_val, y = PAV_log10, color = nutrient, linetype = nutrient)) +
  geom_line() +
  facet_rep_wrap(~ param_foc, scales = "free_x") +
  scale_x_log10() +
  fig_theme

# explore P simulation
pav_q_p <- sim_fun("low", "high", "+P", first_virus = "RPV", 
       plant_time = 12, res_time = 11, inv_time = max_days-11-12)


#### RPV q simulations ####

V0_b <- V0_c <- V0_init
rpv_q_params <- param_df_rpv %>%
  mutate(first_virus = "PAV") %>%
  mutate(pav_conc = pmap(., function(param_foc, param_val, first_virus) param_fun(params_def, param_foc, param_val, first_virus))) %>%
  unnest(cols = c(pav_conc)) %>%
  mutate(nutrient = fct_relevel(nutrient, "low", "+N", "+P"))

ggplot(rpv_q_params, 
       aes(x = param_val, y = RPV_log10, color = nutrient, linetype = nutrient)) +
  geom_line() +
  facet_rep_wrap(~ param_foc, scales = "free_x") +
  scale_x_log10() +
  fig_theme



