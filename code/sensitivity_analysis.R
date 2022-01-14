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
    mutate(VbH_low = V_b_low * H_low,
           VbH_n = V_b_n * H_n,
           VbH_p = V_b_p * H_p,
           VbH_np = V_b_np * H_np,
           VcH_low = V_c_low * H_low,
           VcH_n = V_c_n * H_n,
           VcH_p = V_c_p * H_p,
           VcH_np = V_c_np * H_np) %>%
    pivot_longer(cols = -time,
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
           value = if_else(str_starts(variable, "V") == T, 
                                 log10(value + 1), value),
           abund_type = if_else(str_starts(variable, "V") == T, 
                                str_sub(variable2, 5, 7), NA_character_))
  
  # save time series
  print(ggplot(mod_out2, aes(x = time, y = value, 
                             color = nutrient, linetype = nutrient)) +
          geom_line() +
          facet_wrap(~ variable2, scales = "free_y") +
          scale_color_viridis_d() +
          fig_theme +
          labs(title = paste0(param_foc1, " = ", param_val1, "; ", param_foc2, " = ", param_val2)))
  
  # last time point
  mod_out3 <- mod_out2 %>%
    filter(time == max(time))
  
  # return output
  return(mod_out3)

}

# test function
virus2_model_sim(params_def2, "RPV", V0_b = V0, V0_c = V0,
                 plant_time = 11, res_time = 12, 
                 inv_time = 100-11-12) %>%
  virus2_model_format(params_def2) %>%
  filter(time == max(time) & variable2 %in% c("PAV_conc", "PAV_pop")) %>%
  mutate(virus_conc = log10(value + 1))

pdf("output/temp_sensitivity_analysis_fig.pdf")
param_fun(params_def2, "r_b", 0.0961, "r_c", 0.221, "RPV") %>%
  filter(variable2 %in% c("PAV_conc", "PAV_pop"))
dev.off()

pdf("output/temp_sensitivity_analysis_fig.pdf")
param_fun(params_def2, "a_n_lo", 1.1e-10, "r_c", 0.221, "RPV") %>%
  filter(variable2 %in% c("PAV_conc", "PAV_pop"))
dev.off()


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

# values for z and q
z_nc_vals <- 10^seq(-18, 0, length.out = 10)
q_nc_vals2 <- 10^(-4:5) # simulations crash for low q's

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
  filter(variable2 == "RPV_pop") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = virus_abund)) +
  geom_point(size = 10) +
  facet_wrap(~ nutrient) +
  scale_color_viridis_c() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "RPV N conc", y = "RPV min N")
# for lowest min N, increasing N conc decreases RPV pop
# for higher min N, increasing N conc allows RPV to invade

pav_inv_qz_nc %>%
  filter(variable2 == "PAV_pop") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = virus_abund)) +
  geom_point(size = 10) +
  facet_wrap(~ nutrient) +
  scale_color_viridis_c() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "RPV N conc", y = "RPV min N")
# similar patterns for PAV

# time series to clarify why
params_qz_nc <- params_def2
params_qz_nc["z_nc"] <- z_nc_vals[1]
params_qz_nc["q_nc"] <- q_nc_vals2[1]

sim_qz_nc <- virus2_model_sim(params_qz_nc, "RPV", V0_b = V0, V0_c = V0, 
                              plant_time = 11, res_time = 12, 
                              inv_time = 100-11-12) %>%
  virus2_model_format(params_qz_nc)

sim_qz_nc %>%
  filter(variable2 == "RPV_pop") %>%
  ggplot(aes(time, log10(value + 1), 
             color = nutrient, linetype = nutrient, size = lim_nut_H)) +
  geom_line() +
  scale_color_viridis_d() +
  scale_size_manual(values = c(1.2, 0.6)) +
  fig_theme

plant_fig_fun(sim_qz_nc, params_qz_nc, -2e-3, -5e-4)
# increasing z by itself: makes plant N crash > plant biomass declines >
# RPV population size quickly plateaus at smaller size
# increase q: RPV conc doesn't change, pop size is smaller than default >
# plant N stops at virus's q > plant size is bigger than default
# this is an interaction between q and z -- doesn't happen with increased q alone
# only occurs once z is set to 1e-6, this is probably when virus drawdown of Q
# starts to have a large impact. Why does it equilibrate at q_nc though?


#### q time series ####

# change q value
params_q_nc_lo <- params_def2
params_q_nc_lo["q_nc"] <- 1.1e-5
params_q_nc_hi <- params_def2
params_q_nc_hi["q_nc"] <- 1.1e-1

# run simulations
sim_q_nc_lo <- virus2_model_sim(params_q_nc_lo, "RPV", V0_b = 0, V0_c = V0, 
                                plant_time = 11, res_time = 12, 
                                inv_time = 100-11-12) %>%
  virus2_model_format(params_q_nc_lo)

sim_q_nc_df <- virus2_model_sim(params_def2, "RPV", V0_b = 0, V0_c = V0, 
                                plant_time = 11, res_time = 12, 
                                inv_time = 100-11-12) %>%
  virus2_model_format(params_def2)

sim_q_nc_hi <- virus2_model_sim(params_q_nc_hi, "RPV", V0_b = 0, V0_c = V0, 
                                plant_time = 11, res_time = 12, 
                                inv_time = 100-11-12) %>%
  virus2_model_format(params_q_nc_hi)

# combine simulations
sim_q_nc <- sim_q_nc_lo %>%
  mutate(q_nc = "low (1.1e-5)") %>%
  full_join(sim_q_nc_df %>%
              mutate(q_nc = "med (1.1e-3)")) %>%
  full_join(sim_q_nc_hi %>%
              mutate(q_nc = "high (1.1e-1)"))

sim_q_nc %>%
  filter(variable2 == "RPV_pop") %>%
  ggplot(aes(time, log10(value + 1), 
             color = q_nc, size = lim_nut_H)) +
  geom_line() +
  scale_color_viridis_d(option = "cividis") +
  scale_size_manual(values = c(1.2, 0.6)) +
  facet_wrap(~ nutrient) +
  fig_theme


#### q-z interaction effect on Q ####

# values for z and q
z_nc_vals2 <- 10^(-18:-5)
q_nc_vals3 <- 10^(-5:-1)

# data frame
qz_nc_in2 <- tibble(param_foc1 = "z_nc",
                   param_val1 = z_nc_vals2) %>%
  expand_grid(tibble(param_foc2 = "q_nc",
                     param_val2 = q_nc_vals3))

# RPV by itself
pdf("output/sensitivity_analysis_qz_nc.pdf")
rpv_res_qz_nc <- qz_nc_in2 %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2) param_fun(params_def2, param_foc1, param_val1, param_foc2, param_val2, first_virus = "RPV", V0_b = 0))) %>%
  unnest(cols = c(sim_out))
dev.off()
# when z_n is low, RPV does better with low q_n and cannot infect with high q_n
# when z_n is high, RPV kills the plant with low q_n by pulling Q_n down very low, RPV pop then goes to zero

# edit dataframe
rpv_res_qz_nc2 <- rpv_res_qz_nc %>%
  mutate(z_nc_10 = log10(param_val1),
         z_nc_f = as.factor(z_nc_10),
         Q_10 = if_else(variable2 %in% "Q_n", log10(value + 1e-6), NA_real_),
         plant_alive = if_else(variable2 == "H" & value > 0, 1, 0),
         RPV_present = if_else(variable2 == "RPV_pop" & value > 0, 1, 0),
         PAV_present = if_else(variable2 == "PAV_pop" & value > 0, 1, 0)) %>%
  group_by(param_val1, param_val2, nutrient) %>%
  mutate(plant_alive = sum(plant_alive, na.rm = T),
         RPV_present = sum(RPV_present, na.rm = T),
         PAV_present = sum(PAV_present, na.rm = T)) %>%
  ungroup()

# missing values
rpv_res_qz_nc2 %>%
  filter(variable2 == "Q_n" & is.na(value))
# z = 1e-6, q = 1e-5

rpv_res_qz_nc2 %>%
  filter(variable2 == "Q_n" & is.na(Q_10))
# same as above plus z = 1e-5, q = 1e-5 for low (negative Q)

# figure
rpv_res_qz_nc2 %>%
  filter(variable2 == "Q_n" & plant_alive == 1) %>%
  ggplot(aes(x = param_val2, y = value, color = z_nc_10, group = z_nc_f)) +
  geom_line() +
  facet_wrap(~ nutrient, scales = "free") +
  scale_color_viridis_c(option = "magma", name = "z_nc") +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "RPV min N", y = "Plant N conc.") +
  fig_theme
# consistent with expected trend

# why do some decline at higher RPV min N?
# looks like z_nc = 1e-6 or 1e-7, q_nc = 1e-1
# plant-virus competition

# present/absent diagram
# may be better to do this with longer time series
rpv_res_qz_nc3 <- rpv_res_qz_nc2 %>%
  select(param_foc1, param_val1, param_foc2, param_val2, nutrient, plant_alive, RPV_present) %>%
  unique() %>%
  mutate(final_present = case_when(plant_alive == 1 & RPV_present == 1 ~ "plant and virus",
                                   plant_alive == 1 & RPV_present == 0 ~ "plant only",
                                   plant_alive == 0 & RPV_present == 1 ~ "virus only",
                                   plant_alive == 0 & RPV_present == 0 ~ "neither"))

ggplot(rpv_res_qz_nc3, aes(x = param_val1, y = param_val2, color = final_present)) +
  geom_point(size = 10) +
  facet_wrap(~ nutrient) +
  scale_color_viridis_d(option = "plasma") +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "RPV N content", y = "RPV min N") +
  fig_theme 
# when q_nc = Q_n,min (1e-3), plant and virus always coexist
# when q_nc decreases, virus can only persist if it's N content isn't too high (kills plant)
# when q_nc increases, virus can only persist if it's N content is high (controls plant's N content)
# all are "neither" for a simulation that crashed (no values after a certain time point)
# plant and virus present in low bottom right is from a weird simulation (probably an error)

#### start here ####
# build off above to include virus competition