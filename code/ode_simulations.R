## Goal: simulate mutual invasions


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(deSolve)


#### parameters ####

# supply rates
a_n_lo <- 1.1e-6
a_n_hi <- 5.6e-5
a_p_lo <- 1.6e-7
a_p_hi <- 8.2e-6

# other parameters
u_n <- 2.9e-3
u_p <- 1.4e-3
k_n <- 4.9e-5
k_p <- 3.0e-5
Qmin_n <- 1.1e-3
Qmin_p <- 7.4e-5
q_nb <- q_nc <- 1.1e-3
q_pb <- q_pc <- 7.4e-5
z_nb <- 1.6e-18
z_pb <- 2.6e-19
z_nc <- 1.7e-18
z_pc <- 2.6e-19
m <- 0.04
g <- 0.28
c_b <- 0.5
c_c <- 0.8
r_b <- 0.82
r_c <- 1.4

# initial values
R0_n_lo <- a_n_lo*2
R0_p_lo <- a_p_lo*2
R0_n_hi <- a_n_hi*2
R0_p_hi <- a_p_hi*2
Q0_n <- Qmin_n
Q0_p <- Qmin_p
B0 <- 1e-3
V0_b <- V0_c <- 100000


#### model ####

plant_virus_model = function (t, yy, parms) { 
  
  # supply rates
  a_n = parms[1];
  a_p = parms[2];
  
  # set initial values
  R_n = yy[1];
  R_p = yy[2];
  Q_n = yy[3];
  Q_p = yy[4];
  B = yy[5];
  V_b = yy[6];
  V_c = yy[7];
  
  # model
  dR_n = a_n - (u_n * R_n * B) / (R_n + k_n);
  dR_p = a_p - (u_p * R_p * B) / (R_p + k_p);
  dQ_n = (u_n * R_n) / (R_n + k_n) - min((1 - Qmin_n / Q_n), (1 - Qmin_p / Q_p)) * g * Q_n - min((1 - q_nb / Q_n), (1 - q_pb / Q_p)) * z_nb * r_b * V_b - min((1 - q_nc / Q_n), (1 - q_pc / Q_p)) * z_nc * r_c * V_c
  dQ_p = (u_p * R_p) / (R_p + k_p) - min((1 - Qmin_n / Q_n), (1 - Qmin_p / Q_p)) * g * Q_p - min((1 - q_nb / Q_n), (1 - q_pb / Q_p)) * z_pb * r_b * V_b - min((1 - q_nc / Q_n), (1 - q_pc / Q_p)) * z_pc * r_c * V_c
  dB = min((1 - Qmin_n / Q_n), (1 - Qmin_p / Q_p)) * g * B - m * B
  dV_b = min((1 - q_nb / Q_n), (1 - q_pb / Q_p)) * r_b * V_b - c_b * V_b
  dV_c = min((1 - q_nc / Q_n), (1 - q_pc / Q_p)) * r_c * V_c - c_c * V_c
  
  return(list(c(dR_n, dR_p, dQ_n, dQ_p, dB, dV_b, dV_c)))
}


#### simulation wrapper ####

sim_fun <- function(N_level, P_level, nut_trt, first_virus){
  
  # initial resource values
  if(N_level == "low"){
    R0_n <- R0_n_lo
    a_n <- a_n_lo
  }else{
    R0_n <- R0_n_hi
    a_n <- a_n_hi
  }
  
  if(P_level == "low"){
    R0_p <- R0_p_lo
    a_p <- a_p_lo
  }else{
    R0_p <- R0_p_hi
    a_p <- a_p_hi
  }
  
  # plant model
  plant_mod <- ode(c(R_n = R0_n, R_p = R0_p, Q_n = Q0_n, Q_p = Q0_p, B = B0, V_b = 0, V_c = 0),
                   seq(0, 11, length.out = 100), plant_virus_model, c(a_n = a_n, a_p = a_p)) %>%
    as_tibble() %>%
    mutate(across(everything(), as.double))
    
  # final time
  plant_init <- plant_mod %>%
    filter(time == 11)
  
  # virus initial conditions
  y0_first_virus <- c(R_n = pull(plant_init, R_n), 
                      R_p = pull(plant_init, R_p), 
                      Q_n = pull(plant_init, Q_n), 
                      Q_p = pull(plant_init, Q_p), 
                      B = pull(plant_init, B), 
                      V_b = if_else(first_virus == "PAV", V0_b, 0),
                      V_c = if_else(first_virus == "RPV", V0_c, 0))
  
  # first virus model
  first_virus_mod <- ode(y0_first_virus, seq(0, 12, length.out = 100),
                         plant_virus_model, c(a_n = a_n, a_p = a_p)) %>%
    as_tibble() %>%
    mutate(across(everything(), as.double))
  
  # final time
  first_virus_init <- first_virus_mod %>%
    filter(time == 12)
  
  # virus initial conditions
  y0_second_virus <- c(R_n = pull(first_virus_init, R_n), 
                      R_p = pull(first_virus_init, R_p), 
                      Q_n = pull(first_virus_init, Q_n), 
                      Q_p = pull(first_virus_init, Q_p), 
                      B = pull(first_virus_init, B), 
                      V_b = if_else(first_virus == "PAV", pull(first_virus_init, V_b), V0_b),
                      V_c = if_else(first_virus == "RPV", pull(first_virus_init, V_c), V0_c))
  
  # first virus model
  second_virus_mod <- ode(y0_second_virus, seq(0, 19, length.out = 100),
                         plant_virus_model, c(a_n = a_n, a_p = a_p)) %>%
    as_tibble() %>%
    mutate(across(everything(), as.double))
  
  # combine models
  mod_out <- plant_mod %>%
    full_join(first_virus_mod %>%
                mutate(time = time + 11)) %>%
    full_join(second_virus_mod %>%
                mutate(time = time + 11 + 12)) %>%
    mutate(nutrient = nut_trt,
           resident = first_virus,
           invader = if_else(resident == "PAV", "RPV", "PAV"))
  
}


#### simulations ####

low_pav_sim <- sim_fun("low", "low", "low", "PAV")
n_pav_sim <- sim_fun("high", "low", "N", "PAV")
p_pav_sim <- sim_fun("low", "high", "P", "PAV")
np_pav_sim <- sim_fun("high", "high", "N+P", "PAV")

low_rpv_sim <- sim_fun("low", "low", "low", "RPV")
n_rpv_sim <- sim_fun("high", "low", "N", "RPV")
p_rpv_sim <- sim_fun("low", "high", "P", "RPV")
np_rpv_sim <- sim_fun("high", "high", "N+P", "RPV")


#### virus figures ####

virus_dat <- low_pav_sim %>%
  full_join(n_pav_sim) %>%
  full_join(p_pav_sim) %>%
  full_join(np_pav_sim) %>%
  full_join(low_rpv_sim) %>%
  full_join(n_rpv_sim) %>%
  full_join(p_rpv_sim) %>%
  full_join(np_rpv_sim) %>%
  mutate(V_b_rel = V_b / max(V_b),
         V_c_rel = V_c / max(V_c)) %>%
  select(time, nutrient, invader, resident, V_b_rel, V_c_rel) %>%
  pivot_longer(cols = c(V_b_rel, V_c_rel),
               names_to = "virus",
               values_to = "concentration") %>%
  mutate(virus = fct_recode(virus, "BYDV-PAV" = "V_b_rel",
                            "CYDV-RPV" = "V_c_rel"))

# facet by virus
ggplot(virus_dat, aes(time, concentration, linetype = virus, color = nutrient)) +
  geom_line() +
  facet_wrap(~invader)

# facet by nutrient
ggplot(virus_dat, aes(time, concentration, linetype = virus, color = nutrient)) +
  geom_line() +
  facet_wrap(nutrient~invader, scales = "free")
# RPV can growth with N or N+P
# PAV can only growth with N+P