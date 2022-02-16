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

# load model and settings
source("code/model_settings.R")


#### model functions ####

# change parameter value
param_fun <- function(params_in, param_foc1, param_val1, param_foc2, param_val2,
                      first_virus, output_type = "last",
                      V0_b = V0, V0_c = V0, inits_in = init_def2,
                      plant_time = 11, res_time = 12, inv_time = 100-11-12){
  
  # update parameter value
  params_in[param_foc1] <- param_val1
  params_in[param_foc2] <- param_val2
  
  # update resource supply
  # if(param_foc1 == "a_n_lo"){
  #   inits_in["R_n_low"] <- param_val1 * 10
  #   inits_in["R_n_p"] <- param_val1 * 10
  # }
  # if(param_foc2 == "a_p_lo"){
  #   inits_in["R_p_low"] <- param_val2 * 10
  #   inits_in["R_p_n"] <- param_val2 * 10
  # }
  
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
           value2 = if_else(str_starts(variable, "V") == T, 
                                 log10(value + 1), value),
           abund_type = if_else(str_starts(variable, "V") == T, 
                                str_sub(variable2, 5, 7), NA_character_))
  
  # save time series
  print(ggplot(mod_out2, aes(x = time, y = value2, 
                             color = nutrient, linetype = nutrient)) +
          geom_line() +
          facet_wrap(~ variable2, scales = "free_y") +
          scale_color_viridis_d() +
          fig_theme +
          labs(title = paste0(param_foc1, " = ", param_val1, "; ", param_foc2, " = ", param_val2)))
  
  # outputs
  if(output_type == "last"){
    
    # extract last time point
    mod_out3 <- mod_out2 %>%
      filter(time == max(time))
    
  } else {
    
    # calculate virus growth rates
    mod_out3 <- virus2_growth_rate(mod_out2, first_virus)
    
  }
  
  
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
param_fun(params_def2, "q_nb", 1.1e-3, "q_pb", 7.4e-5, "RPV") %>%
  filter(variable2 %in% c("PAV_conc", "PAV_pop"))
dev.off()
# output should match output of above

pdf("output/temp_sensitivity_analysis_fig.pdf")
param_fun(params_def2, "a_n_lo", 1.1e-10, "q_pb", 7.4e-5, "RPV") %>%
  filter(variable2 %in% c("PAV_conc", "PAV_pop"))
dev.off()
# output should match output of above unless N supply is low


#### resource supply rates ####

# values for low nutrient supply
a_n_vals <- c(1.1e-13, 1.1e-10, 1.1e-6, 5.6e-5)
a_p_vals <- c(1.6e-14, 1.6e-11, 1.6e-7, 8.2e-6)

# data frame
a_in <- tibble(param_foc1 = "a_n_lo",
               param_val1 = a_n_vals) %>%
  expand_grid(tibble(param_foc2 = "a_p_lo",
                     param_val2 = a_p_vals))

# PAV invades RPV
pdf("output/sensitivity_analysis_pav_inv_a.pdf")
pav_inv_a <- a_in %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2) param_fun(params_def2, param_foc1, param_val1, param_foc2, param_val2, first_virus = "RPV", output_type = "last"))) %>%
  unnest(cols = c(sim_out))
dev.off()

# figure
pav_inv_a %>%
  filter(variable2 %in% c("PAV_conc", "RPV_conc") & nutrient == "low") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = value2)) +
  geom_point(size = 10) +
  facet_wrap(~ variable2) +
  scale_color_viridis_c(name = "Concentration", direction = -1) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "N supply", y = "P supply")

pav_inv_a %>%
  filter(variable2 %in% c("PAV_pop", "RPV_pop") & nutrient == "low") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = value2)) +
  geom_point(size = 10) +
  facet_wrap(~ variable2) +
  scale_color_viridis_c(name = "Population size", direction = -1) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "N supply", y = "P supply")
# Q is only drawn down so fast, giving viruses some nutrients
# regardless of resource supply rate

# invasion growth rate
pav_inv_a2 <- a_in %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2) param_fun(params_def2, param_foc1, param_val1, param_foc2, param_val2, first_virus = "RPV", output_type = "first"))) %>%
  unnest(cols = c(sim_out))

pav_inv_a2 %>%
  filter(variable2 %in% c("PAV_conc", "RPV_conc") & nutrient == "low") %>%
  mutate(growth_sign = ifelse(growth <= 0, "<= 0", "> 0")) %>%
  ggplot(aes(x = param_val1, y = param_val2, color = growth, shape = growth_sign)) +
  geom_point(size = 10) +
  facet_wrap(~ variable2) +
  scale_color_viridis_c(name = "Growth rate", direction = -1) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "N supply", y = "P supply")


#### later invasion ####

# use default parameters
param_fun(params_def2, "r_b", 0.288, "r_c", 0.451, "RPV", output_type = "first",
          plant_time = 11, res_time = 100-11, inv_time = 100) %>%
  filter(variable2 %in% c("PAV_conc", "PAV_pop"))
# virus can't invade later, but decline in concentration is tiny


#### minimum N concentrations (q_n) ####

# values for low nutrient supply
q_nb_vals <- 10^seq(-3, 6, length.out = 5)
q_nc_vals <- 10^seq(-3, 6, length.out = 5)

# data frame
q_n_in <- tibble(param_foc1 = "q_nb",
                 param_val1 = q_nb_vals) %>%
  expand_grid(tibble(param_foc2 = "q_nc",
                     param_val2 = q_nc_vals))

# PAV invades RPV
pdf("output/sensitivity_analysis_pav_inv_q_n.pdf")
pav_inv_q_n <- q_n_in %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2) param_fun(params_def2, param_foc1, param_val1, param_foc2, param_val2, first_virus = "RPV", output_type = "last"))) %>%
  unnest(cols = c(sim_out))
dev.off()

# figure
pav_inv_q_n %>%
  filter(variable2 == "PAV_conc") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = value2)) +
  geom_point(size = 10) +
  facet_wrap(~ nutrient) +
  scale_color_viridis_c(name = "PAV conc", direction = -1) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "PAV min N", y = "RPV min N")
# increasing the virus's q above the plants makes it grow slower, no matter how much larger

pav_inv_q_n %>%
  filter(variable2 == "RPV_conc") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = value2)) +
  geom_point(size = 10) +
  facet_wrap(~ nutrient) +
  scale_color_viridis_c(name = "RPV conc", direction = -1) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "PAV min N", y = "RPV min N")


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