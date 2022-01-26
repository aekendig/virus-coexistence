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
pdf("output/sensitivity_analysis_pav_inv_a.pdf")
pav_inv_a <- a_in %>%
  mutate(first_virus = "RPV") %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2, first_virus) param_fun(params_def2, param_foc1, param_val1, param_foc2, param_val2, first_virus))) %>%
  unnest(cols = c(sim_out))
dev.off()

# figure
pav_inv_a %>%
  filter(variable2 == "PAV_conc") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = value)) +
  geom_point(size = 10) +
  facet_wrap(~ nutrient) +
  scale_color_viridis_c() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "N supply", y = "P supply")

pav_inv_a %>%
  filter(variable2 == "PAV_pop") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = value)) +
  geom_point(size = 10) +
  facet_wrap(~ nutrient) +
  scale_color_viridis_c() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "N supply", y = "P supply")

pav_inv_a %>%
  filter(variable2 == "RPV_pop") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = value)) +
  geom_point(size = 10) +
  facet_wrap(~ nutrient) +
  scale_color_viridis_c() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "N supply", y = "P supply")

# why are the viruses always present?
# looking at plant dynamics pdf:
# Q is only drawn down so fast, giving viruses some nutrients


#### later invasion ####

# use default parameters
param_fun(params_def2, "r_b", 0.0961, "r_c", 0.221, "RPV",
          plant_time = 11, res_time = 100-11, inv_time = 100) %>%
  filter(variable2 %in% c("PAV_conc", "PAV_pop"))
# abundances are comparable to very low nutrients
# viruses can always invade when Q is at Qmin because of their q values


#### minimum nutrient concentrations (q) ####

# values for low nutrient supply
q_nb_vals <- 10^(-3:6)
q_nc_vals <- 10^(-3:6)

# data frame
q_n_in <- tibble(param_foc1 = "q_nb",
                 param_val1 = q_nb_vals) %>%
  expand_grid(tibble(param_foc2 = "q_nc",
                     param_val2 = q_nc_vals))

# PAV invades RPV
pdf("output/sensitivity_analysis_pav_inv_q_n.pdf")
pav_inv_q_n <- q_n_in %>%
  mutate(first_virus = "RPV") %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2, first_virus) param_fun(params_def2, param_foc1, param_val1, param_foc2, param_val2, first_virus))) %>%
  unnest(cols = c(sim_out))
dev.off()
# increasing q makes virus go extinct
# the time to extinction decreases with larger q
# this is about the virus-plant interaction, not virus-virus

# figure
pav_inv_q_n %>%
  filter(variable2 == "PAV_conc") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = value)) +
  geom_point(size = 10) +
  facet_wrap(~ nutrient) +
  scale_color_viridis_c() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "PAV min N", y = "RPV min N")

pav_inv_q_n %>%
  filter(variable2 == "RPV_conc") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = value)) +
  geom_point(size = 10) +
  facet_wrap(~ nutrient) +
  scale_color_viridis_c() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "PAV min N", y = "RPV min N")


#### minimum nutrient concentration (q_n)/nutrient content of virus (z_n) ####

# values for z and q
z_nc_vals <- 10^seq(-18, 0, length.out = 10)
q_nc_vals2 <- 10^(-4:5) # simulations crash for low q's

# data frame
qz_nc_in <- tibble(param_foc1 = "q_nc",
                   param_val1 = q_nc_vals2) %>%
  expand_grid(tibble(param_foc2 = "z_nc",
                   param_val2 = z_nc_vals))

# PAV invades RPV
pdf("output/sensitivity_analysis_pav_inv_qz_nc.pdf")
pav_inv_qz_nc <- qz_nc_in %>%
  mutate(first_virus = "RPV") %>%
  mutate(sim_out = pmap(., function(param_foc1, param_val1, param_foc2, param_val2, first_virus) param_fun(params_def2, param_foc1, param_val1, param_foc2, param_val2, first_virus))) %>%
  unnest(cols = c(sim_out))
dev.off()
# lowest q (RPV min N) 1e-4:
# increasing z (RPV N conc) up to 1e-4 gives RPV more control of Q_n
# Q_n decreases -> host and PAV decrease
# this occurs because 1e-4 is below the host's Q_min (1.1e-3)
# when z >= 1e-4 -> Q_n < 0 (unrealistic/simulation error) -> P becomes the limiting nutrient
# q = 1e-3:
# I think Q_n is pulled to 1e-3 with higher z, decreasing plant growth (because 1e-3 < 1.1e-3), but not killing the plant (in this time frame)
# q = 1e-2:
# RPV can't persist until it's z is high enough to influence Q
# PAV and host can persist because their q's are lower

# figure
pav_inv_qz_nc %>%
  filter(variable2 == "RPV_pop") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = value)) +
  geom_point(size = 10) +
  facet_wrap(~ nutrient) +
  scale_color_viridis_c() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "RPV min N (q)", y = "RPV N conc (z)")

pav_inv_qz_nc %>%
  filter(variable2 == "PAV_pop") %>%
  ggplot(aes(x = param_val1, y = param_val2, color = value)) +
  geom_point(size = 10) +
  facet_wrap(~ nutrient) +
  scale_color_viridis_c() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "RPV min N (q)", y = "RPV N conc (z)")
# PAV can always persist unless plant is killed


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
  geom_hline(yintercept = as.numeric(params_def1["Qmin_n"]),
             color = "black", linetype = "dashed") +
  geom_line() +
  facet_wrap(~ nutrient) +
  scale_x_log10() +
  scale_color_viridis_d(name = "RPV N conc (z)") +
  labs(x = "PAV N conc (z)", y = "Relative virus concentration") +
  fig_theme

#### start here ####
# look at saved PDF and above figure
# make figure for plant biomass
# when/why do plants benefit?
# does this make sense for the model?
# revise last simulation to prevent extreme plant biomass values


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