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
V0_b <- V0_c <- V0_init
max_days <- 100


#### model function ####

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
  mod <- sim_fun("high", "high", "N+P", first_virus = first_virus, 
                 plant_time = 12, res_time = 11, inv_time = max_days-11-12)
  
  # output: invading virus titer
  if(first_virus == "PAV"){
    output <- mod %>%
      filter(time == max_days) %>%
      pull(V_c)
  }else{
    output <- mod %>%
      filter(time == max_days) %>%
      pull(V_b)
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


#### parameter ranges ####

n <- 10

pav_params <- tibble(param_foc = rep(c("q_nb", "q_pb", "q_nc", "q_pc", 
                                       "z_nb", "z_pb", "z_nc", "z_pc", 
                                       "c_b", "c_c", 
                                       "r_b", "r_c"), each = n), 
                     param_val = c(rep(seq(1e-6, 1e-2, length.out = n), 4),
                                   rep(seq(1e-20, 1e-17, length.out = n), 4),
                                   rep(seq(1e-4, 1e-2, length.out = n), 2),
                                   rep(seq(1e-3, 1, length.out = n), 2)),
                     first_virus = "RPV") %>%
  mutate(pav_conc = pmap(., function(param_foc, param_val, first_virus) param_fun(params_def, param_foc, param_val, first_virus))) %>%
  unnest(cols = c(pav_conc))

ggplot(pav_params, aes(x = param_val, y = log10(pav_conc + 1))) +
  geom_line() +
  facet_rep_wrap(~ param_foc, scales = "free_x") +
  fig_theme
