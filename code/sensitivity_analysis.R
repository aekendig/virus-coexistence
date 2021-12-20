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

# simulation data formatting
# revised from ode_simulations
sim_dat_fun <- function(low_sim, n_sim, p_sim, np_sim){
  
  dat_out <- low_sim %>%
    full_join(n_sim) %>%
    full_join(p_sim) %>%
    full_join(np_sim) %>%
    filter(time == max_days) %>%
    mutate(PAV_log10 = log10(V_b + 1),
           RPV_log10 = log10(V_c + 1)) %>%
    arrange(nutrient)
  
  return(dat_out)
  
}

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
                                 7.4e-11, 7.4e-10, 7.4e-9, 7.4e-8, 7.4e-7, 7.4e-6, 7.4e-5))

param_df_rpv <- tibble(param_foc = rep(c("q_nb", "q_pb"), each = n), 
                       param_val = c(1.1e-9, 1.1e-8, 1.1e-7, 1.1e-6, 1.1e-5, 1.1e-4, 1.1e-3,
                                     7.4e-11, 7.4e-10, 7.4e-9, 7.4e-8, 7.4e-7, 7.4e-6, 7.4e-5))


#### PAV q simulations ####

V0_b <- V0_c <- V0_init
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

#### start here ####
# why does only one q out of the four have an effect?