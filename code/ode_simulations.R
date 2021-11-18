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
m <- 0.05
g <- 0.29
c_b <- 0.02
c_c <- 0.1
r_b <- 0.28
r_c <- 0.63

# initial values
H0 <- 1e-3
Q0_const <- 10
Q0_n <- Qmin_n * Q0_const
Q0_p <- Qmin_p * Q0_const
R0_const <- 3
R0_n_lo <- a_n_hi * R0_const
R0_n_hi <- a_n_hi * R0_const
R0_p_lo <- a_p_hi * R0_const
R0_p_hi <- a_p_hi * R0_const

# figure settings
fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing.x = unit(0,"line"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 10),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.background = element_blank(),
        legend.key.width = unit(7, "mm"),
        legend.key.heigh = unit(5, "mm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 11, hjust = 0),
        strip.placement = "outside",
        plot.title = element_text(size = 11, vjust = 0))

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
  plant_mod <- ode(c(R_n = R0_n, R_p = R0_p, Q_n = Q0_n, Q_p = Q0_p, H = H0, V_b = 0, V_c = 0),
                   seq(0, plant_days, length.out = 100), plant_virus_model, c(a_n = a_n, a_p = a_p)) %>%
    as_tibble() %>%
    mutate(across(everything(), as.double))
    
  # final time
  plant_init <- plant_mod %>%
    filter(time == plant_days)
  
  # virus initial conditions
  y0_first_virus <- c(R_n = pull(plant_init, R_n), 
                      R_p = pull(plant_init, R_p), 
                      Q_n = pull(plant_init, Q_n), 
                      Q_p = pull(plant_init, Q_p), 
                      H = pull(plant_init, H), 
                      V_b = if_else(first_virus == "PAV", V0_b, 0),
                      V_c = if_else(first_virus == "RPV", V0_c, 0))
  
  # first virus model
  first_virus_mod <- ode(y0_first_virus, seq(0, res_days, length.out = 100),
                         plant_virus_model, c(a_n = a_n, a_p = a_p)) %>%
    as_tibble() %>%
    mutate(across(everything(), as.double))
  
  # final time
  first_virus_init <- first_virus_mod %>%
    filter(time == res_days)
  
  # virus initial conditions
  y0_second_virus <- c(R_n = pull(first_virus_init, R_n),
                      R_p = pull(first_virus_init, R_p),
                      Q_n = pull(first_virus_init, Q_n),
                      Q_p = pull(first_virus_init, Q_p),
                      H = pull(first_virus_init, H),
                      V_b = if_else(first_virus == "PAV", pull(first_virus_init, V_b), V0_b),
                      V_c = if_else(first_virus == "RPV", pull(first_virus_init, V_c), V0_c))
  
  # first virus model
  second_virus_mod <- ode(y0_second_virus, seq(0, inv_days, length.out = 100),
                         plant_virus_model, c(a_n = a_n, a_p = a_p)) %>%
    as_tibble() %>%
    mutate(across(everything(), as.double))
  
  # combine models
  mod_out <- plant_mod %>%
    full_join(first_virus_mod %>%
                mutate(time = time + plant_days)) %>%
    full_join(second_virus_mod %>%
                mutate(time = time + plant_days + res_days)) %>%
    mutate(nutrient = nut_trt,
           resident = first_virus,
           invader = if_else(resident == "PAV", "RPV", "PAV"))
  
}


#### long-term plant ####

# time intervals
plant_days <- 250
res_days <- 12
inv_days <- 500-plant_days-res_days
# ran with 1000 days and system equilibrates by 500

# viruses
V0_b <- V0_c <- 0

# simulations
low_plant_sim <- sim_fun("low", "low", "low", "PAV")
n_plant_sim <- sim_fun("high", "low", "N", "PAV")
p_plant_sim <- sim_fun("low", "high", "P", "PAV")
np_plant_sim <- sim_fun("high", "high", "N+P", "PAV")


#### long-term plant figure ####

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
plant_fig_fun <- function(low_sim, n_sim, p_sim, np_sim, 
                          nut_leg, lim_leg){
  
  # combine
  plant_dat <- sim_dat_fun(low_sim, n_sim, p_sim, np_sim)
  
  # panels
  plant_Rn_fig <- ggplot(plant_dat, aes(time, R_n, linetype = nutrient, color = nutrient)) +
    geom_line(aes(size = lim_nut_H)) +
    scale_color_viridis_d(end = 0.7) +
    scale_size_manual(values = c(1.2, 0.6)) +
    labs(x = "Time (days)", y = "Environment N (g)", title = "(C)") +
    fig_theme +
    theme(legend.position = "none")
  
  plant_Rp_fig <- ggplot(plant_dat, aes(time, R_p, linetype = nutrient, color = nutrient)) +
    geom_line(aes(size = lim_nut_H)) +
    scale_color_viridis_d(end = 0.7) +
    scale_size_manual(values = c(1.2, 0.6)) +
    labs(x = "Time (days)", y = "Environment P (g)", title = "(D)") +
    fig_theme +
    theme(legend.position = "none")
  
  plant_H_fig <- ggplot(plant_dat, aes(time, H, linetype = nutrient, color = nutrient)) +
    geom_line(aes(size = lim_nut_H)) +
    scale_color_viridis_d(end = 0.7, guide = "none") +
    scale_linetype(guide = "none") +
    scale_size_manual(values = c(1.2, 0.6), name = "Limiting nutrient",
                      labels = c("N", "P")) +
    labs(x = "Time (days)", y = "Plant biomass (g)", title = "(B)") +
    fig_theme +
    theme(legend.position = nut_leg)
  
  plant_Qn_fig <- ggplot(plant_dat, aes(time, Q_n)) +
    geom_hline(yintercept = Qmin_n, color = "black") +
    geom_text(data = Qn_lab, aes(label = label), parse = T, color = "black", fontface = "italic",
              size = 3, hjust = 0, vjust = 0, nudge_y = -3e-3) +
    geom_line(aes(linetype = nutrient, color = nutrient, size = lim_nut_H)) +
    scale_color_viridis_d(end = 0.7) +
    scale_size_manual(values = c(1.2, 0.6)) +
    labs(x = "Time (days)", y = "Plant N concentration", title = "(E)") +
    fig_theme +
    theme(legend.position = "none")
  
  plant_Qp_fig <- ggplot(plant_dat, aes(time, Q_p)) +
    geom_hline(yintercept = Qmin_p, color = "black") +
    geom_text(data = Qp_lab, aes(label = label), parse = T, color = "black", fontface = "italic", 
              size = 3, hjust = 0, vjust = 0, nudge_y = -1e-3) +
    geom_line(aes(linetype = nutrient, color = nutrient, size = lim_nut_H)) +
    scale_color_viridis_d(end = 0.7) +
    scale_size_manual(values = c(1.2, 0.6)) +
    labs(x = "Time (days)", y = "Plant P concentration", title = "(F)") +
    fig_theme +
    theme(legend.position = "none")
  
  plant_lim_fig <- ggplot(plant_dat, aes(time, Qlim, linetype = nutrient, color = nutrient)) +
    geom_line(aes(size = lim_nut_H)) +
    scale_color_viridis_d(end = 0.7, name = "Nutrient supply") +
    scale_linetype(name = "Nutrient supply") +
    scale_size_manual(values = c(1.2, 0.6), guide= "none") +
    labs(x = "Time (days)", title = "(A)", 
         y = expression(paste("Limiting nutrient ratio (", Q[min], "/Q)", sep = ""))) +
    fig_theme +
    theme(legend.position = lim_leg) +
    guides(color = guide_legend(override.aes = list(size = 1.2)))
  
  # combine figures
  plant_top_fig <- plot_grid(plant_lim_fig, plant_H_fig, plant_Rn_fig,  
                             nrow = 1)
  plant_bot_fig <- plot_grid(plant_Rp_fig, plant_Qn_fig, plant_Qp_fig, 
                             nrow = 1)
  
  # output
  return(plot_grid(plant_top_fig, plant_bot_fig,
                   nrow = 2))
  
}

# apply function
pdf("output/long_term_plant_simulation_figure.pdf", width = 6.5, height = 4.5)
plant_fig_fun(low_plant_sim, n_plant_sim, p_plant_sim, np_plant_sim,
              nut_leg = c(0.5, 0.5), lim_leg = c(0.5, 0.4))
dev.off()


#### single virus ####

# time intervals
plant_days <- 250
res_days <- 12
inv_days <- 5000-plant_days-res_days

# PAV with plant
V0_b <- 100000
V0_c <- 0

low_pav_only_sim <- sim_fun("low", "low", "low", "PAV")
n_pav_only_sim <- sim_fun("high", "low", "N", "PAV")
p_pav_only_sim <- sim_fun("low", "high", "P", "PAV")
np_pav_only_sim <- sim_fun("high", "high", "N+P", "PAV")

# RPV with plant
V0_b <- 0
V0_c <- 100000

low_rpv_only_sim <- sim_fun("low", "low", "low", "RPV")
n_rpv_only_sim <- sim_fun("high", "low", "N", "RPV")
p_rpv_only_sim <- sim_fun("low", "high", "P", "RPV")
np_rpv_only_sim <- sim_fun("high", "high", "N+P", "RPV")


#### single virus figures ####

# plant figures
pdf("output/long_term_pav_only_plant_simulation_figure.pdf", width = 6.5, height = 4.5)
plant_fig_fun(low_pav_only_sim, n_pav_only_sim, p_pav_only_sim, np_pav_only_sim,
              nut_leg = c(0.5, 0.5), lim_leg = c(0.5, 0.4))
dev.off()

pdf("output/long_term_rpv_only_plant_simulation_figure.pdf", width = 6.5, height = 4.5)
plant_fig_fun(low_rpv_only_sim, n_rpv_only_sim, p_rpv_only_sim, np_rpv_only_sim,
              nut_leg = c(0.5, 0.5), lim_leg = c(0.5, 0.4))
dev.off()

# PAV only figure
pav_only_sim <- sim_dat_fun(low_pav_only_sim, n_pav_only_sim, 
                            p_pav_only_sim, np_pav_only_sim)

ggplot(pav_only_sim, aes(x = time, y = PAV_log10, color = nutrient, linetype = nutrient)) +
  geom_line(aes(size = lim_nut_H)) +
  scale_color_viridis_d(end = 0.7, name = "Nutrient supply") +
  scale_linetype(name = "Nutrient supply") +
  scale_size_manual(values = c(1.2, 0.6), guide= "none") +
  labs(x = "Time (days)", title = "(A) Single virus", 
       y = expression(paste("BYDV-PAV titer (", log[10], ")", sep = ""))) +
  fig_theme 

# RPV only figure
rpv_only_sim <- sim_dat_fun(low_rpv_only_sim, n_rpv_only_sim, 
                            p_rpv_only_sim, np_rpv_only_sim)

ggplot(rpv_only_sim, aes(x = time, y = RPV_log10, color = nutrient, linetype = nutrient)) +
  geom_line(aes(size = lim_nut_H)) +
  scale_color_viridis_d(end = 0.7, name = "Nutrient supply") +
  scale_linetype(name = "Nutrient supply") +
  scale_size_manual(values = c(1.2, 0.6), guide= "none") +
  labs(x = "Time (days)", title = "(B)", 
       y = expression(paste("CYDV-RPV titer (", log[10], ")", sep = ""))) +
  fig_theme 


#### start here ####

# PAV can't establish when host has only been growing for 11 days (maybe okay?)
# RPV titer increases immediately in +N simulation before it's introduced (coding error?)
# Both viruses reach unreasonable titers (re-fit parameters)
# Viruses kill the plants (okay)
# Ran simulations for 500 days and the plant output exactly matched plant alone


#### virus invasion ####

# time intervals
plant_days <- 11
res_days <- 12
inv_days <- 19

# viruses
V0_b <- V0_c <- 100000

# simulations
low_pav_sim <- sim_fun("low", "low", "low", "PAV")
n_pav_sim <- sim_fun("high", "low", "N", "PAV")
p_pav_sim <- sim_fun("low", "high", "P", "PAV")
np_pav_sim <- sim_fun("high", "high", "N+P", "PAV")

low_rpv_sim <- sim_fun("low", "low", "low", "RPV")
n_rpv_sim <- sim_fun("high", "low", "N", "RPV")
p_rpv_sim <- sim_fun("low", "high", "P", "RPV")
np_rpv_sim <- sim_fun("high", "high", "N+P", "RPV")


#### virus figure ####



# combine
comb_dat <- low_pav_sim %>%
  full_join(n_pav_sim) %>%
  full_join(p_pav_sim) %>%
  full_join(np_pav_sim) %>%
  full_join(low_rpv_sim) %>%
  full_join(n_rpv_sim) %>%
  full_join(p_rpv_sim) %>%
  full_join(np_rpv_sim)

# edit virus data
virus_dat <- comb_dat %>%
  mutate(PAV_rel = V_b / max(V_b),
         RPV_rel = V_c / max(V_c)) %>%
  rename(PAV_conc = V_b, RPV_conc = V_c) %>%
  select(time, nutrient, invader, resident, PAV_conc, RPV_conc, PAV_rel, RPV_rel) %>%
  pivot_longer(cols = c(PAV_conc, RPV_conc, PAV_rel, RPV_rel),
               names_to = c("virus", ".value"),
               names_pattern = "(.+)_(.+)") %>%
  mutate(virus = fct_recode(virus, "BYDV-PAV" = "PAV",
                            "CYDV-RPV" = "RPV"),
         nutrient = fct_recode(nutrient, "+N" = "N",
                               "+P" = "P",
                               "+N+P" = "N+P") %>%
           fct_relevel("low", "+N", "+P"),
         log_conc = log10(conc + 1),
         dpiR = time - plant_days)

res_dat <- virus_dat %>%
  filter((virus == "BYDV-PAV" & invader == "RPV") | (virus == "CYDV-RPV" & invader == "PAV")) %>%
  mutate(scenario = case_when(invader == "PAV" ~ "(A) CYDV-RPV resident",
                              invader == "RPV" ~ "(B) BYDV-PAV resident")) %>%
  filter(dpiR > 0)

inv_dat <- virus_dat %>%
  filter((virus == "BYDV-PAV" & invader == "PAV") | (virus == "CYDV-RPV" & invader == "RPV")) %>%
  mutate(scenario = case_when(invader == "PAV" ~ "(C) BYDV-PAV invades",
                              invader == "RPV" ~ "(D) CYDV-RPV invades")) %>%
  filter(dpiR > 0)

# virus figures
res_fig <- ggplot(res_dat, aes(dpiR, log_conc, linetype = nutrient, color = nutrient)) +
  geom_line(size = 1.2) +
  facet_wrap(~scenario, scales = "free") +
  scale_color_viridis_d(end = 0.7, name = "Nutrient") +
  scale_linetype(name = "Nutrient") +
  labs(x = "Time (days)", y = expression(paste("Resident virus titer (", log[10], ")", sep = ""))) +
  scale_y_continuous(limits = c(0, 8.2)) +
  fig_theme +
  theme(legend.position = c(0.1, 0.3),
        axis.title.x = element_blank())

inv_fig <- ggplot(inv_dat, aes(dpiR, log_conc, linetype = nutrient, color = nutrient)) +
  geom_line(size = 1.2) +
  facet_wrap(~scenario, scales = "free") +
  scale_color_viridis_d(end = 0.7) +
  labs(x = "Time (days)", y = expression(paste("Invading virus titer (", log[10], ")", sep = ""))) +
  scale_y_continuous(limits = c(0, 8.2)) +
  fig_theme +
  theme(legend.position = "none")

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