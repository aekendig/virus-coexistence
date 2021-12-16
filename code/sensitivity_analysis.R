## Goal: sensitivity analysis of virus parameters


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(deSolve)
library(cowplot)
library(lemon)
library(sensitivity)

# load model and settings
source("code/model_settings.R")
V0_b <- V0_c <- V0_init


#### model function ####

param_fun <- function(params, first_virus){
  
  # extract parameter values
  q_nb <- params["q_nb"]
  q_pb <- params["q_pb"]
  q_nc <- params["q_nc"]
  q_pc <- params["q_pc"]
  z_nb <- params["z_nb"]
  z_pb <- params["z_pb"]
  z_nc <- params["z_nc"]
  z_pc <- params["z_pc"]
  c_b <- params["c_b"]
  c_c <- params["c_c"]
  r_b <- params["r_b"]
  r_c <- params["r_c"]
  
  # run model
  mod <- sim_fun("high", "high", "N+P", first_virus = first_virus, 
                 plant_time = 12, res_time = 11, inv_time = 100-11-12)
  
  # output: invading virus titer
  if(first_virus == "PAV"){
    output <- mod %>%
      filter(time == 100) %>%
      pull(V_c)
  }else{
    output <- mod %>%
      filter(time == 100) %>%
      pull(V_b)
  }
  
  # return output
  return(output)

}


# test function
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

pav_inv <- param_fun(params_def, "RPV")
rpv_inv <- param_fun(params_def, "PAV")

#### start with long-term invasion simulations in ode_simulations ####
# need stable pop size


#### monotonic relationships ####


#### PRCC ####