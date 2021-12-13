## Goal: simulate mutual invasions


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(deSolve)
library(cowplot)


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
m <- 0.007
g <- 0.144
c_b <- 0.002
c_c <- 0.005
r_b <- 0.09
r_c <- 0.27

# initial values
H0 <- 1e-2
Q0_const <- 10
Q0_n <- Qmin_n * Q0_const
Q0_p <- Qmin_p * Q0_const
R0_const <- 10
R0_n_lo <- a_n_hi * R0_const
R0_n_hi <- a_n_hi * R0_const
R0_p_lo <- a_p_hi * R0_const
R0_p_hi <- a_p_hi * R0_const
V0_init <- 5e5

# figure settings
fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing.x = unit(0,"line"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.background = element_blank(),
        legend.key.width = unit(7, "mm"),
        legend.key.heigh = unit(5, "mm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 8, hjust = 0, face = "bold"),
        strip.placement = "outside",
        plot.title = element_text(size = 8, vjust = 0, face = "bold"))

# figure labels
Qn_lab <- tibble(time = 0, Q_n = Qmin_n, label = "Q['min,N']")
Qp_lab <- tibble(time = 0, Q_p = Qmin_p, label = "Q['min,P']")


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
  H = yy[5];
  V_b = yy[6];
  V_c = yy[7];
  
  # model
  dR_n = a_n - (u_n * R_n * H) / (R_n + k_n);
  dR_p = a_p - (u_p * R_p * H) / (R_p + k_p);
  dQ_n = (u_n * R_n) / (R_n + k_n) - min((1 - Qmin_n / Q_n), (1 - Qmin_p / Q_p)) * g * Q_n - min((1 - q_nb / Q_n), (1 - q_pb / Q_p)) * z_nb * r_b * V_b - min((1 - q_nc / Q_n), (1 - q_pc / Q_p)) * z_nc * r_c * V_c
  dQ_p = (u_p * R_p) / (R_p + k_p) - min((1 - Qmin_n / Q_n), (1 - Qmin_p / Q_p)) * g * Q_p - min((1 - q_nb / Q_n), (1 - q_pb / Q_p)) * z_pb * r_b * V_b - min((1 - q_nc / Q_n), (1 - q_pc / Q_p)) * z_pc * r_c * V_c
  dH = min((1 - Qmin_n / Q_n), (1 - Qmin_p / Q_p)) * g * H - m * H
  dV_b = min((1 - q_nb / Q_n), (1 - q_pb / Q_p)) * r_b * V_b - c_b * V_b
  dV_c = min((1 - q_nc / Q_n), (1 - q_pc / Q_p)) * r_c * V_c - c_c * V_c
  
  return(list(c(dR_n, dR_p, dQ_n, dQ_p, dH, dV_b, dV_c)))
}


#### simulation wrapper ####

sim_fun <- function(N_level, P_level, nut_trt, first_virus, plant_time, res_time, inv_time){
  
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
  plant_mod <- ode(c(R_n = R0_n, R_p = R0_p, Q_n = Q0_n, Q_p = Q0_p, H = H0, V_b = 0, V_c = 0),
                   seq(0, plant_time, length.out = 100), plant_virus_model, c(a_n = a_n, a_p = a_p)) %>%
    as_tibble() %>%
    mutate(across(everything(), as.double))
    
  # final time
  plant_init <- plant_mod %>%
    filter(time == plant_time)
  
  # virus initial conditions
  y0_first_virus <- c(R_n = pull(plant_init, R_n), 
                      R_p = pull(plant_init, R_p), 
                      Q_n = pull(plant_init, Q_n), 
                      Q_p = pull(plant_init, Q_p), 
                      H = pull(plant_init, H), 
                      V_b = if_else(first_virus == "PAV", V0_b, 0),
                      V_c = if_else(first_virus == "RPV", V0_c, 0))
  
  # first virus model
  first_virus_mod <- ode(y0_first_virus, seq(0, res_time, length.out = 100),
                         plant_virus_model, c(a_n = a_n, a_p = a_p)) %>%
    as_tibble() %>%
    mutate(across(everything(), as.double))
  
  # final time
  first_virus_init <- first_virus_mod %>%
    filter(time == res_time)
  
  # virus initial conditions
  y0_second_virus <- c(R_n = pull(first_virus_init, R_n),
                      R_p = pull(first_virus_init, R_p),
                      Q_n = pull(first_virus_init, Q_n),
                      Q_p = pull(first_virus_init, Q_p),
                      H = pull(first_virus_init, H),
                      V_b = if_else(first_virus == "PAV", pull(first_virus_init, V_b), V0_b),
                      V_c = if_else(first_virus == "RPV", pull(first_virus_init, V_c), V0_c))
  
  # first virus model
  second_virus_mod <- ode(y0_second_virus, seq(0, inv_time, length.out = 100),
                         plant_virus_model, c(a_n = a_n, a_p = a_p)) %>%
    as_tibble() %>%
    mutate(across(everything(), as.double))
  
  # combine models
  mod_out <- plant_mod %>%
    full_join(first_virus_mod %>%
                mutate(time = time + plant_time)) %>%
    full_join(second_virus_mod %>%
                mutate(time = time + plant_time + res_time)) %>%
    mutate(nutrient = nut_trt,
           resident = first_virus,
           invader = if_else(resident == "PAV", "RPV", "PAV"))
  
}


#### figure functions ####

# simulation data formatting
sim_dat_fun <- function(low_sim, n_sim, p_sim, np_sim){
  
  dat_out <- low_sim %>%
    full_join(n_sim) %>%
    full_join(p_sim) %>%
    full_join(np_sim) %>%
    mutate(Qlim_n = Qmin_n / Q_n,
           Qlim_p = Qmin_p / Q_p) %>%
    rowwise() %>%
    mutate(Qlim = max(Qlim_n, Qlim_p)) %>%
    ungroup() %>%
    mutate(nutrient = fct_recode(nutrient, "+N" = "N",
                                 "+P" = "P",
                                 "+N+P" = "N+P") %>%
             fct_relevel("low", "+N", "+P"),
           lim_nut_H = case_when(Qlim == Qlim_n ~ "Q[N]",
                                 Qlim == Qlim_p ~ "Q[P]"),
           PAV_log10 = log10(V_b + 1),
           RPV_log10 = log10(V_c + 1)) %>%
    rename(PAV = V_b, RPV = V_c)
  
  return(dat_out)
  
}

# plant_figure function
plant_fig_fun <- function(low_sim, n_sim, p_sim, np_sim, q_adj_n, q_adj_p){
  
  # combine
  plant_dat <- sim_dat_fun(low_sim, n_sim, p_sim, np_sim)
  
  # panels
  plant_H_fig <- ggplot(plant_dat, aes(time, H, linetype = nutrient, color = nutrient)) +
    geom_line(aes(size = lim_nut_H)) +
    scale_color_viridis_d(end = 0.7, name = "Nutrient\nsupply") +
    scale_linetype(name = "Nutrient\nsupply") +
    scale_size_manual(values = c(1.2, 0.6), guide= "none") +
    labs(x = "Time (days)", y = "Plant biomass (g)", title = "(A)") +
    fig_theme +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 1.2)))
  
  plant_lim_fig <- ggplot(plant_dat, aes(time, Qlim, linetype = nutrient, color = nutrient)) +
    geom_line(aes(size = lim_nut_H)) +
    scale_color_viridis_d(end = 0.7, guide = "none") +
    scale_linetype(guide = "none") +
    scale_size_manual(values = c(1.2, 0.6), name = "Limiting\nnutrient",
                      labels = c("N", "P")) +
    labs(x = "Time (days)", title = "(B)", 
         y = expression(paste("Limiting nutrient ratio (", Q[min], "/Q)", sep = ""))) +
    fig_theme +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank())
  
  plant_Rn_fig <- ggplot(plant_dat, aes(time, R_n, linetype = nutrient, color = nutrient)) +
    geom_line(aes(size = lim_nut_H)) +
    scale_color_viridis_d(end = 0.7) +
    scale_size_manual(values = c(1.2, 0.6)) +
    labs(x = "Time (days)", y = "Environment N (g)", title = "(C)") +
    fig_theme +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_blank())
  
  plant_Rp_fig <- ggplot(plant_dat, aes(time, R_p, linetype = nutrient, color = nutrient)) +
    geom_line(aes(size = lim_nut_H)) +
    scale_color_viridis_d(end = 0.7) +
    scale_size_manual(values = c(1.2, 0.6)) +
    labs(x = "", y = "Environment P (g)", title = "(D)") +
    fig_theme +
    theme(legend.position = "none")
  
  plant_Qn_fig <- ggplot(plant_dat, aes(time, Q_n)) +
    geom_hline(yintercept = Qmin_n, color = "black") +
    geom_text(data = Qn_lab, aes(label = label), parse = T, color = "black", fontface = "italic",
              size = 3, hjust = 0, vjust = 0, nudge_y = q_adj_n) +
    geom_line(aes(linetype = nutrient, color = nutrient, size = lim_nut_H)) +
    scale_color_viridis_d(end = 0.7) +
    scale_size_manual(values = c(1.2, 0.6)) +
    labs(x = "Time (days)", y = "Plant N concentration", title = "(E)") +
    fig_theme +
    theme(legend.position = "none")
  
  plant_Qp_fig <- ggplot(plant_dat, aes(time, Q_p)) +
    geom_hline(yintercept = Qmin_p, color = "black") +
    geom_text(data = Qp_lab, aes(label = label), parse = T, color = "black", fontface = "italic", 
              size = 3, hjust = 0, vjust = 0, nudge_y = q_adj_p) +
    geom_line(aes(linetype = nutrient, color = nutrient, size = lim_nut_H)) +
    scale_color_viridis_d(end = 0.7) +
    scale_size_manual(values = c(1.2, 0.6)) +
    labs(x = "", y = "Plant P concentration", title = "(F)") +
    fig_theme +
    theme(legend.position = "none")
  
  # extract legends
  nut_leg <- get_legend(plant_H_fig)
  lim_leg <- get_legend(plant_lim_fig)
  
  # combine figures
  plant_top_fig <- plot_grid(plant_H_fig + theme(legend.position = "none"), 
                             plant_lim_fig + theme(legend.position = "none"), 
                             plant_Rn_fig,
                             nut_leg,
                             nrow = 1,
                             rel_widths = c(1, 1, 1, 0.35))
  plant_bot_fig <- plot_grid(plant_Rp_fig, plant_Qn_fig, plant_Qp_fig, 
                             lim_leg,
                             nrow = 1,
                             rel_widths = c(1, 1, 1, 0.35))
  
  # output
  return(plot_grid(plant_top_fig, plant_bot_fig,
                   nrow = 2, rel_heights = c(0.85, 1)))
  
}


#### long-term plant ####

# time intervals
plant_days <- 250
res_days <- 12
inv_days <- 1000-plant_days-res_days

# viruses
V0_b <- V0_c <- 0

# simulations
low_plant_long_sim <- sim_fun("low", "low", "low", "PAV", plant_days, res_days, inv_days)
n_plant_long_sim <- sim_fun("high", "low", "N", "PAV", plant_days, res_days, inv_days)
p_plant_long_sim <- sim_fun("low", "high", "P", "PAV", plant_days, res_days, inv_days)
np_plant_long_sim <- sim_fun("high", "high", "N+P", "PAV", plant_days, res_days, inv_days)

# figure
pdf("output/long_term_plant_simulation_figure.pdf", width = 6.5, height = 4.5)
plant_fig_fun(low_plant_long_sim, n_plant_long_sim, p_plant_long_sim, np_plant_long_sim, -3e-3, -1e-3)
dev.off()


#### short-term plant ####

# time intervals
plant_days <- 11
res_days <- 12
inv_days <- 19

# viruses
V0_b <- V0_c <- 0

# simulations
low_plant_short_sim <- sim_fun("low", "low", "low", "PAV", plant_days, res_days, inv_days)
n_plant_short_sim <- sim_fun("high", "low", "N", "PAV", plant_days, res_days, inv_days)
p_plant_short_sim <- sim_fun("low", "high", "P", "PAV", plant_days, res_days, inv_days)
np_plant_short_sim <- sim_fun("high", "high", "N+P", "PAV", plant_days, res_days, inv_days)

# figure
pdf("output/short_term_plant_simulation_figure.pdf", width = 6.5, height = 4.5)
plant_fig_fun(low_plant_short_sim, n_plant_short_sim, p_plant_short_sim, np_plant_short_sim, -2e-3, -5e-4)
dev.off()


#### single virus simulations ####

# time intervals
plant_days <- 11
res_days <- 12
inv_days <- 19

# PAV with plant
V0_b <- V0_init
V0_c <- 0

low_pav_first_sim <- sim_fun("low", "low", "low", "PAV", plant_days, res_days, inv_days)
n_pav_first_sim <- sim_fun("high", "low", "N", "PAV", plant_days, res_days, inv_days)
p_pav_first_sim <- sim_fun("low", "high", "P", "PAV", plant_days, res_days, inv_days)
np_pav_first_sim <- sim_fun("high", "high", "N+P", "PAV", plant_days, res_days, inv_days)

low_pav_second_sim <- sim_fun("low", "low", "low", "RPV", plant_days, res_days, inv_days)
n_pav_second_sim <- sim_fun("high", "low", "N", "RPV", plant_days, res_days, inv_days)
p_pav_second_sim <- sim_fun("low", "high", "P", "RPV", plant_days, res_days, inv_days)
np_pav_second_sim <- sim_fun("high", "high", "N+P", "RPV", plant_days, res_days, inv_days)

# RPV with plant
V0_b <- 0
V0_c <- V0_init

low_rpv_first_sim <- sim_fun("low", "low", "low", "RPV", plant_days, res_days, inv_days)
n_rpv_first_sim <- sim_fun("high", "low", "N", "RPV", plant_days, res_days, inv_days)
p_rpv_first_sim <- sim_fun("low", "high", "P", "RPV", plant_days, res_days, inv_days)
np_rpv_first_sim <- sim_fun("high", "high", "N+P", "RPV", plant_days, res_days, inv_days)

low_rpv_second_sim <- sim_fun("low", "low", "low", "PAV", plant_days, res_days, inv_days)
n_rpv_second_sim <- sim_fun("high", "low", "N", "PAV", plant_days, res_days, inv_days)
p_rpv_second_sim <- sim_fun("low", "high", "P", "PAV", plant_days, res_days, inv_days)
np_rpv_second_sim <- sim_fun("high", "high", "N+P", "PAV", plant_days, res_days, inv_days)


#### virus invasion simulations ####

# time intervals
plant_days <- 11
res_days <- 12
inv_days <- 19

# viruses
V0_b <- V0_c <- V0_init

# simulations
low_pav_inv_sim <- sim_fun("low", "low", "low", "RPV", plant_days, res_days, inv_days)
n_pav_inv_sim <- sim_fun("high", "low", "N", "RPV", plant_days, res_days, inv_days)
p_pav_inv_sim <- sim_fun("low", "high", "P", "RPV", plant_days, res_days, inv_days)
np_pav_inv_sim <- sim_fun("high", "high", "N+P", "RPV", plant_days, res_days, inv_days)

low_rpv_inv_sim <- sim_fun("low", "low", "low", "PAV", plant_days, res_days, inv_days)
n_rpv_inv_sim <- sim_fun("high", "low", "N", "PAV", plant_days, res_days, inv_days)
p_rpv_inv_sim <- sim_fun("low", "high", "P", "PAV", plant_days, res_days, inv_days)
np_rpv_inv_sim <- sim_fun("high", "high", "N+P", "PAV", plant_days, res_days, inv_days)


#### edit simulation data ####

pav_first_sim <- sim_dat_fun(low_pav_first_sim, n_pav_first_sim, 
                             p_pav_first_sim, np_pav_first_sim) %>%
  filter(time > plant_days) %>%
  mutate(virus_conc = PAV_log10,
         virus = resident,
         scenario = "(D) Resident alone")

pav_second_sim <- sim_dat_fun(low_pav_second_sim, n_pav_second_sim, 
                              p_pav_second_sim, np_pav_second_sim) %>%
  filter(time > (plant_days + res_days)) %>%
  mutate(virus_conc = PAV_log10,
         virus = invader,
         scenario = "(B) Invade alone")

rpv_first_sim <- sim_dat_fun(low_rpv_first_sim, n_rpv_first_sim, 
                             p_rpv_first_sim, np_rpv_first_sim) %>%
  filter(time > plant_days) %>%
  mutate(virus_conc = RPV_log10,
         virus = resident,
         scenario = "(D) Resident alone")

rpv_second_sim <- sim_dat_fun(low_rpv_second_sim, n_rpv_second_sim, 
                              p_rpv_second_sim, np_rpv_second_sim) %>%
  filter(time > (plant_days + res_days)) %>%
  mutate(virus_conc = RPV_log10,
         virus = invader,
         scenario = "(B) Invade alone")

pav_res_sim <- sim_dat_fun(low_rpv_inv_sim, n_rpv_inv_sim, 
                           p_rpv_inv_sim, np_rpv_inv_sim) %>%
  filter(time > plant_days) %>%
  mutate(virus_conc = PAV_log10,
         virus = resident,
         scenario = "(C) Resident virus")

rpv_inv_sim <- pav_res_sim %>%
  filter(time > (plant_days + res_days)) %>%
  mutate(virus_conc = RPV_log10,
         virus = invader,
         scenario = "(A) Invade resident virus")

rpv_res_sim <- sim_dat_fun(low_pav_inv_sim, n_pav_inv_sim, 
                           p_pav_inv_sim, np_pav_inv_sim) %>%
  filter(time > plant_days) %>%
  mutate(virus_conc = RPV_log10,
         virus = resident,
         scenario = "(C) Resident virus")

pav_inv_sim <- rpv_res_sim %>%
  filter(time > (plant_days + res_days)) %>%
  mutate(virus_conc = PAV_log10,
         virus = invader,
         scenario = "(A) Invade resident virus")

# combine
virus_sim <- pav_inv_sim %>%
  full_join(rpv_inv_sim) %>%
  full_join(pav_second_sim) %>%
  full_join(rpv_second_sim) %>%
  full_join(pav_res_sim) %>%
  full_join(rpv_res_sim) %>%
  full_join(pav_first_sim) %>%
  full_join(rpv_first_sim)


#### virus figure ####

# facet labels
fac_lab <- c(PAV = "BYDV-PAV~titer~(log[10])",
             RPV = "BYDV-PAV~titer~(log[10])")

#### start here ####
# add in axes
# move y labels to left
# only x strip text is bold
# parse labels

ggplot(virus_sim, 
       aes(x = time, y = virus_conc, color = nutrient, linetype = nutrient)) +
  geom_line(aes(size = lim_nut_H)) +
  facet_grid(rows = vars(virus),
             cols = vars(scenario),
             labeller = labeller(virus = fac_lab)) +
  scale_color_viridis_d(end = 0.7, name = "Nutrient\nsupply") +
  scale_linetype(name = "Nutrient\nsupply") +
  scale_size_manual(values = c(1.2, 0.6), name = "Limiting\nnutrient",
                    labels = c("N", "P")) +
  labs(x = "Time (days)") +
  fig_theme +
  theme(axis.title.y = element_blank())



#### single virus figures ####

# PAV only figures


pav_first_fig <- ggplot(pav_first_sim, 
                        aes(x = time, y = PAV_log10, color = nutrient, linetype = nutrient)) +
  geom_line(aes(size = lim_nut_H)) +
  scale_color_viridis_d(end = 0.7, name = "Nutrient supply") +
  scale_linetype(name = "Nutrient supply") +
  scale_size_manual(values = c(1.2, 0.6), guide= "none") +
  labs(x = "Time (days)", title = "(A) Single virus", 
       y = expression(paste("BYDV-PAV titer (", log[10], ")", sep = ""))) +
  fig_theme 



pav_second_fig <- ggplot(pav_second_sim, 
                         aes(x = time, y = PAV_log10, color = nutrient, linetype = nutrient)) +
  geom_line(aes(size = lim_nut_H)) +
  scale_color_viridis_d(end = 0.7, name = "Nutrient supply") +
  scale_linetype(name = "Nutrient supply") +
  scale_size_manual(values = c(1.2, 0.6), guide= "none") +
  labs(title = "(B) Invade alone", 
       y = expression(paste("BYDV-PAV titer (", log[10], ")", sep = ""))) +
  fig_theme +
  theme(legend.position = "none",
        axis.title = element_blank())

# RPV only figures


rpv_first_fig <- ggplot(rpv_first_sim, 
                        aes(x = time, y = RPV_log10, color = nutrient, linetype = nutrient)) +
  geom_line(aes(size = lim_nut_H)) +
  scale_color_viridis_d(end = 0.7, name = "Nutrient supply") +
  scale_linetype(name = "Nutrient supply") +
  scale_size_manual(values = c(1.2, 0.6), guide= "none") +
  labs(x = "Time (days)", title = "(B)", 
       y = expression(paste("CYDV-RPV titer (", log[10], ")", sep = ""))) +
  fig_theme 



rpv_second_fig <- ggplot(rpv_second_sim, 
                         aes(x = time, y = RPV_log10, color = nutrient, linetype = nutrient)) +
  geom_line(aes(size = lim_nut_H)) +
  scale_color_viridis_d(end = 0.7, name = "Nutrient supply") +
  scale_linetype(name = "Nutrient supply") +
  scale_size_manual(values = c(1.2, 0.6), guide= "none") +
  labs(x = "Time (days)", 
       y = expression(paste("CYDV-RPV titer (", log[10], ")", sep = ""))) +
  fig_theme 





#### virus invasion figures ####

# edit data


# pav invades
pav_inv_fig <- ggplot(pav_inv_sim %>% filter(time > (plant_days + res_days)), 
       aes(x = time, y = PAV_log10, color = nutrient, linetype = nutrient)) +
  geom_line(aes(size = lim_nut_H)) +
  scale_color_viridis_d(end = 0.7, name = "Nutrient supply") +
  scale_linetype(name = "Nutrient supply") +
  scale_size_manual(values = c(1.2, 0.6), guide= "none") +
  labs(title = "(A) Invade other virus", 
       y = expression(paste("BYDV-PAV titer (", log[10], ")", sep = ""))) +
  fig_theme +
  theme(axis.title.x = element_blank())

# rpv invades
rpv_inv_fig <- ggplot(rpv_inv_sim %>% filter(time > (plant_days + res_days)), 
       aes(x = time, y = RPV_log10, color = nutrient, linetype = nutrient)) +
  geom_line(aes(size = lim_nut_H)) +
  scale_color_viridis_d(end = 0.7, name = "Nutrient supply") +
  scale_linetype(name = "Nutrient supply") +
  scale_size_manual(values = c(1.2, 0.6), guide= "none") +
  labs(x = "Time (days)",
       y = expression(paste("CYDV-RPV titer (", log[10], ")", sep = ""))) +
  fig_theme +
  theme(legend.position = "none")

# PAV resident
pav_res_fig <- ggplot(rpv_inv_sim, 
       aes(x = time, y = PAV_log10, color = nutrient, linetype = nutrient)) +
  geom_line(aes(size = lim_nut_H)) +
  scale_color_viridis_d(end = 0.7, name = "Nutrient supply") +
  scale_linetype(name = "Nutrient supply") +
  scale_size_manual(values = c(1.2, 0.6), guide= "none") +
  labs(x = "Time (days)", title = "(A) Single virus", 
       y = expression(paste("BYDV-PAV titer (", log[10], ")", sep = ""))) +
  fig_theme 

repv_res_fig <- ggplot(pav_inv_sim, 
       aes(x = time, y = RPV_log10, color = nutrient, linetype = nutrient)) +
  geom_line(aes(size = lim_nut_H)) +
  scale_color_viridis_d(end = 0.7, name = "Nutrient supply") +
  scale_linetype(name = "Nutrient supply") +
  scale_size_manual(values = c(1.2, 0.6), guide= "none") +
  labs(x = "Time (days)", title = "(A) Single virus", 
       y = expression(paste("CYDV-RPV titer (", log[10], ")", sep = ""))) +
  fig_theme 


# combine
tiff("output/invasion_simulation_figure.tiff", width = 4.5, height = 4.5, units = "in", res = 300)
plot_grid(res_fig, inv_fig,
          nrow = 2,
          rel_heights = c(0.95, 1))
dev.off()


#### plant figure with PAV data ####

# edit plant dat (change name, plant_dat used above)
plant_dat <- comb_dat %>%
  mutate(Qlim_n = Qmin_n / Q_n,
         Qlim_p = Qmin_p / Q_p,
         qlim_n = q_nb / Q_n, # used PAV q's, but they are equal to RPV's
         qlim_p = q_pb / Q_p) %>%
  rowwise() %>%
  mutate(Qlim = max(Qlim_n, Qlim_p),
         qlim = max(qlim_n, qlim_p)) %>%
  ungroup() %>%
  mutate(nutrient = fct_recode(nutrient, "+N" = "N",
                               "+P" = "P",
                               "+N+P" = "N+P") %>%
           fct_relevel("low", "+N", "+P"),
         lim_nut_H = case_when(Qlim == Qlim_n ~ "Q[N]",
                               Qlim == Qlim_p ~ "Q[P]"),
         lim_nut_V = case_when(qlim == qlim_n ~ "Q[N]",
                               qlim == qlim_p ~ "Q[P]"))
  
pav_invader_dat <- plant_dat %>%
  filter(invader == "PAV")

rpv_invader_dat <- plant_dat %>%
  filter(invader == "RPV")


# plant figures
# show one invader in main text and other in supplement because they look about the same
pav_Rn_fig <- ggplot(pav_invader_dat, aes(time, R_n, linetype = nutrient, color = nutrient)) +
  geom_line(aes(size = lim_nut_V)) +
  scale_color_viridis_d(end = 0.7) +
  scale_size_manual(values = c(1.2, 0.6)) +
  labs(x = "Time (days)", y = "Environment N (g)", title = "(C)") +
  fig_theme +
  theme(legend.position = "none")

pav_Rp_fig <- ggplot(pav_invader_dat, aes(time, R_p, linetype = nutrient, color = nutrient)) +
  geom_line(aes(size = lim_nut_V)) +
  scale_color_viridis_d(end = 0.7) +
  scale_size_manual(values = c(1.2, 0.6)) +
  labs(x = "Time (days)", y = "Environment P (g)", title = "(D)") +
  fig_theme +
  theme(legend.position = "none")

pav_H_fig <- ggplot(pav_invader_dat, aes(time, H, linetype = nutrient, color = nutrient)) +
  geom_line(aes(size = lim_nut_V)) +
  scale_color_viridis_d(end = 0.7, name = "Nutrient") +
  scale_linetype(name = "Nutrient") +
  scale_size_manual(values = c(1.2, 0.6), guide = "none") +
  labs(x = "Time (days)", y = "Plant biomass (g)", title = "(B)") +
  fig_theme +
  theme(legend.position = c(0.23, 0.68)) +
  guides(color = guide_legend(override.aes = list(size = 1.2)))

pav_Qn_fig <- ggplot(pav_invader_dat, aes(time, Q_n)) +
  geom_hline(yintercept = Qmin_n, color = "black") +
  geom_text(data = Qn_lab, aes(label = label), parse = T, color = "black", fontface = "italic",
            size = 3, hjust = 0, vjust = 0, nudge_y = 1e-4) +
  geom_line(aes(linetype = nutrient, color = nutrient, size = lim_nut_V)) +
  scale_color_viridis_d(end = 0.7) +
  scale_size_manual(values = c(1.2, 0.6)) +
  labs(x = "Time (days)", y = "Plant N concentration", title = "(E)") +
  fig_theme +
  theme(legend.position = "none")

pav_Qp_fig <- ggplot(pav_invader_dat, aes(time, Q_p)) +
  geom_hline(yintercept = Qmin_p, color = "black") +
  geom_text(data = Qp_lab, aes(label = label), parse = T, color = "black", fontface = "italic", 
            size = 3, hjust = 0, vjust = 0, nudge_y = 4e-5) +
  geom_line(aes(linetype = nutrient, color = nutrient, size = lim_nut_V)) +
  scale_color_viridis_d(end = 0.7) +
  scale_size_manual(values = c(1.2, 0.6)) +
  labs(x = "Time (days)", y = "Plant P concentration", title = "(F)") +
  fig_theme +
  theme(legend.position = "none")

pav_lim_fig <- ggplot(pav_invader_dat, aes(time, qlim, linetype = nutrient, color = nutrient)) +
  geom_line(aes(size = lim_nut_V)) +
  scale_color_viridis_d(end = 0.7, guide = "none") +
  scale_linetype(guide = "none") +
  scale_size_manual(values = c(1.2, 0.6), name = "Limiting nutrient",
                    labels = c("N", "P")) +
  labs(x = "Time (days)", title = "(A)", 
       y = expression(paste("Limiting nutrient ratio (", Q[min], "/Q)", sep = ""))) +
  fig_theme +
  theme(legend.position = c(0.3, 0.8))

ggplot(pav_invader_dat, aes(time, Qlim, linetype = nutrient, color = nutrient)) +
  geom_line(aes(size = lim_nut_H)) +
  scale_color_viridis_d(end = 0.7) +
  scale_size_manual(values = c(1.2, 0.6)) +
  labs(x = "Time (days)", title = "(F)", 
       y = expression(paste("Limiting nutrient ratio (", Q[min], "/Q)", sep = ""))) +
  fig_theme +
  theme(legend.position = "none")
# same as the virus one because the qmin values were set equal to Qmin

# combine figures
pav_top_fig <- plot_grid(pav_lim_fig, pav_H_fig, pav_Rn_fig,  
                         nrow = 1)
pav_bot_fig <- plot_grid(pav_Rp_fig, pav_Qn_fig, pav_Qp_fig, 
                         nrow = 1)

tiff("output/pav_invasion_plant_simulation_figure.tiff", width = 6.5, height = 4.5, units = "in", res = 300)
plot_grid(pav_top_fig, pav_bot_fig,
          nrow = 2)
dev.off()


#### plant figure with RPV data ####

# recreate figures
rpv_Rn_fig <- pav_Rn_fig %+% rpv_invader_dat
rpv_Rp_fig <- pav_Rp_fig %+% rpv_invader_dat
rpv_H_fig <- pav_H_fig %+% rpv_invader_dat
rpv_Qn_fig <- pav_Qn_fig %+% rpv_invader_dat
rpv_Qp_fig <- pav_Qp_fig %+% rpv_invader_dat
rpv_lim_fig <- pav_lim_fig %+% rpv_invader_dat

# combine figures
rpv_top_fig <- plot_grid(rpv_lim_fig, rpv_H_fig, rpv_Rn_fig,  
                         nrow = 1)
rpv_bot_fig <- plot_grid(rpv_Rp_fig, rpv_Qn_fig, rpv_Qp_fig, 
                         nrow = 1)

tiff("output/rpv_invasion_plant_simulation_figure.tiff", width = 6.5, height = 4.5, units = "in", res = 300)
plot_grid(rpv_top_fig, rpv_bot_fig,
          nrow = 2)
dev.off()