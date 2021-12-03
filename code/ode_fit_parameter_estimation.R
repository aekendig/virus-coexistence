## Goal: estimate plant and virus growth rates


#### set up ####

# import data
# data repository for Kendig et al. 2020:
# https://doi.org/10.6073/pasta/00a35cbd4a9b2a007433c3d2be0d1742
# saved as a directory in this project
# created sub-directories for code and data

# clears environment, loads tidyverse, processes qPCR data
source("./edi.411.2/code/qPCR_raw_data_processing.R") 
# edited the above code so that data were imported from correct folder

# clear all except dataset
rm(list = setdiff(ls(), c("dat")))

# load packages
library(FME)
library(manipulate)
library(minpack.lm)
library(cowplot)
library(lemon)

# import data
sdat <- read_csv("./edi.411.2/data/sample_exp_molc_data.csv")

# figure settings
fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing.x = unit(0,"line"),
        axis.text.y = element_text(size = 7, color = "black"),
        axis.text.x = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 9),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.position = "none",
        legend.key.size = unit(5, "mm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 9),
        strip.placement = "outside",
        plot.title = element_text(size = 9, vjust = 0))


#### edit data ####
# code source: concentration_analysis.R in edi folder

# days post inoculation
dpi <- tibble(
  time = 1:8,
  dpi = c(5, 8, 12, 16, 19, 22, 26, 29))

# remove samples:
# poor standard curve efficiency
# quantities below standard curve, but greater than 1e3 (not sure if these should be zeros or not; standards removed if contamination had higher concentration)
# multiple qPCR tests of the same sample and the sample wasn't detected in one or had the higher variance in detected in multiple
# specific cases: low volume, known contamination, mis-labelling
# make values below the standard curve 0 (will be removed in this analysis)
dat2 <- dat %>%
  filter(remove == 0 & material == "shoot") %>%
  mutate(quant_adj = case_when(target == "PAV" & quant_adj < PAVmin ~ 0,
                               target == "RPV" & quant_adj < RPVmin ~ 0,
                               is.na(quant_adj) ~ 0,
                               TRUE ~ quant_adj)) %>%
  full_join(dpi)

# remove values above standard curve
dat3 <- dat2 %>%
  filter((target == "RPV" & quant_adj <= RPVmax) | (target == "PAV" & quant_adj <= PAVmax))

# average technical replicates
dat4 <- dat3 %>%
  group_by(target, dpi, time, inoc, high_N, high_P, nutrient, round, replicate, sample, shoot_mass_g, root_mass_g, leaf_area_mm2, leaves, mass_ext_mg, PAVmin, PAVint, PAVslope, RPVmin, RPVint, RPVslope, q_group) %>%
  summarise(tech_cycle = mean(cycle, na.rm = T)) %>%
  mutate(quant = case_when(
    target == "PAV" ~ 10 ^ ((tech_cycle - PAVint) / PAVslope),
    target == "RPV" ~ 10 ^ ((tech_cycle - RPVint) / RPVslope)),
    quant_adj = case_when(target == "PAV" & quant < PAVmin ~ 0,
                          target == "RPV" & quant < RPVmin ~ 0,
                          is.na(quant) ~ 0,
                          TRUE ~ quant),
    quant_ul = case_when(q_group %in% c("01","02") ~ quant_adj / 2.5,
                         TRUE ~ quant_adj / 7),
    quant_zero = case_when(quant_adj == 0 ~ 1,
                           TRUE ~ 0)) %>%
  ungroup()

# check for same sample in multiple qPCR groups
dups <- dat4 %>%
  group_by(target, sample) %>%
  mutate(dup = duplicated(sample)) %>%
  filter(dup == T) %>%
  select(sample, target) %>%
  ungroup()

dups %>%
  left_join(dat4) %>%
  select(sample, target, quant_ul) %>%
  data.frame() 
# all are zero's, none are healthy

# select rows from sdat
sdat2 <- sdat %>%
  filter(material == "shoot") %>%
  select(round, time, inoc, nutrient, replicate, RTPCR_PAV, RTPCR_RPV)

# accidental RPV inoculations
(accr <- dat4 %>%
    left_join(sdat2) %>%
    filter(inoc %in% c("PAV", "healthy") & ((target == "RPV" & quant_zero == 0) | RTPCR_RPV == 1)) %>%
    select(sample, round, time, inoc, nutrient, replicate, quant_zero, RTPCR_RPV))

# accidental PAV inoculations  
(accp <- dat4 %>%
    left_join(sdat2) %>%
    filter(inoc %in% c("RPV", "healthy") & ((target == "PAV" & quant_zero == 0) | RTPCR_PAV == 1)) %>%
    select(sample, round, time, inoc, nutrient, replicate, quant_zero, RTPCR_PAV)) 

# accidentally inoculated samples
accs <- unique(c(accr$sample, accp$sample))

# make wide by qPCR target
# remove accidental infections
# remove undetected infections
dat5 <- dat4 %>%
  select(target:RPVslope, quant_ul) %>%
  pivot_wider(names_from = target,
              values_from = quant_ul,
              names_glue = "quant_{target}") %>% # quant_target = quantity / ul
  filter(!(sample %in% accs) & 
           !(inoc %in% c("PAV", "coinfection") & quant_PAV == 0) & 
           !(inoc %in% c("RPV", "coinfection") & quant_RPV == 0)) %>%
  mutate(conc_PAV = 1000 * quant_PAV * 50 / mass_ext_mg, # convert to quantity / g
         conc_RPV = 1000 * quant_RPV * 50 / mass_ext_mg,
         co = ifelse(inoc == "coinfection", 1, 0),
         dpp = case_when(round == 1 ~ dpi + 10,
                         TRUE ~ dpi + 11),
         full_mass_g = shoot_mass_g + root_mass_g)

# check for duplicates
dat5[duplicated(dat5$sample) == T, ]
# all were removed


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
q_n <- 1.1e-3 # same for both viruses
q_p <- 7.4e-5 # same for both viruses
z_nb <- 1.6e-18
z_pb <- 2.6e-19
z_nc <- 1.7e-18
z_pc <- 2.6e-19
 
# initial values
H0 <- 1e-2 # changed this 11/24/21 from 1e-3
Q0_const <- 10
Q0_n <- Qmin_n * Q0_const
Q0_p <- Qmin_p * Q0_const
R0_const <- 10 # used to be 3, which was based on nutrients supplied, but nutrients are also in the seed
R0_n_lo <- a_n_hi * R0_const
R0_n_hi <- a_n_hi * R0_const
R0_p_lo <- a_p_hi * R0_const
R0_p_hi <- a_p_hi * R0_const
V0_init <- 100000

# time
plant_days <- 11


#### plant model ####

plant_model = function (t, yy, parms) { 
  
  # supply rates
  m = parms[1];
  g = ifelse(length(parms) == 2, parms[2], g);
  
  # set initial values
  R_n_low = yy[1];
  R_p_low = yy[2];
  Q_n_low = yy[3];
  Q_p_low = yy[4];
  H_low = yy[5];
  
  R_n_n = yy[6];
  R_p_n = yy[7];
  Q_n_n = yy[8];
  Q_p_n = yy[9];
  H_n = yy[10];
  
  R_n_p = yy[11];
  R_p_p = yy[12];
  Q_n_p = yy[13];
  Q_p_p = yy[14];
  H_p = yy[15];
  
  R_n_np = yy[16];
  R_p_np = yy[17];
  Q_n_np = yy[18];
  Q_p_np = yy[19];
  H_np = yy[20];
  
  # model
  dR_n_low = a_n_lo - (u_n * R_n_low * H_low) / (R_n_low + k_n);
  dR_p_low = a_p_lo - (u_p * R_p_low * H_low) / (R_p_low + k_p);
  dQ_n_low = (u_n * R_n_low) / (R_n_low + k_n) - min((1 - Qmin_n / Q_n_low), (1 - Qmin_p / Q_p_low)) * g * Q_n_low
  dQ_p_low = (u_p * R_p_low) / (R_p_low + k_p) - min((1 - Qmin_n / Q_n_low), (1 - Qmin_p / Q_p_low)) * g * Q_p_low
  dH_low = min((1 - Qmin_n / Q_n_low), (1 - Qmin_p / Q_p_low)) * g * H_low - m * H_low
  
  dR_n_n = a_n_hi - (u_n * R_n_n * H_n) / (R_n_n + k_n);
  dR_p_n = a_p_lo - (u_p * R_p_n * H_n) / (R_p_n + k_p);
  dQ_n_n = (u_n * R_n_n) / (R_n_n + k_n) - min((1 - Qmin_n / Q_n_n), (1 - Qmin_p / Q_p_n)) * g * Q_n_n
  dQ_p_n = (u_p * R_p_n) / (R_p_n + k_p) - min((1 - Qmin_n / Q_n_n), (1 - Qmin_p / Q_p_n)) * g * Q_p_n
  dH_n = min((1 - Qmin_n / Q_n_n), (1 - Qmin_p / Q_p_n)) * g * H_n - m * H_n
  
  dR_n_p = a_n_lo - (u_n * R_n_p * H_p) / (R_n_p + k_n);
  dR_p_p = a_p_hi - (u_p * R_p_p * H_p) / (R_p_p + k_p);
  dQ_n_p = (u_n * R_n_p) / (R_n_p + k_n) - min((1 - Qmin_n / Q_n_p), (1 - Qmin_p / Q_p_p)) * g * Q_n_p
  dQ_p_p = (u_p * R_p_p) / (R_p_p + k_p) - min((1 - Qmin_n / Q_n_p), (1 - Qmin_p / Q_p_p)) * g * Q_p_p
  dH_p = min((1 - Qmin_n / Q_n_p), (1 - Qmin_p / Q_p_p)) * g * H_p - m * H_p
  
  dR_n_np = a_n_hi - (u_n * R_n_np * H_np) / (R_n_np + k_n);
  dR_p_np = a_p_hi - (u_p * R_p_np * H_np) / (R_p_np + k_p);
  dQ_n_np = (u_n * R_n_np) / (R_n_np + k_n) - min((1 - Qmin_n / Q_n_np), (1 - Qmin_p / Q_p_np)) * g * Q_n_np
  dQ_p_np = (u_p * R_p_np) / (R_p_np + k_p) - min((1 - Qmin_n / Q_n_np), (1 - Qmin_p / Q_p_np)) * g * Q_p_np
  dH_np = min((1 - Qmin_n / Q_n_np), (1 - Qmin_p / Q_p_np)) * g * H_np - m * H_np
  
  return(list(c(dR_n_low, dR_p_low, dQ_n_low, dQ_p_low, dH_low,
                dR_n_n, dR_p_n, dQ_n_n, dQ_p_n, dH_n,
                dR_n_p, dR_p_p, dQ_n_p, dQ_p_p, dH_p,
                dR_n_np, dR_p_np, dQ_n_np, dQ_p_np, dH_np)))
}


#### visualize plant parameters ####

# data
mock <- dat5 %>%
  filter(inoc == "healthy") %>%
  mutate(variable = case_when(nutrient == "low" ~ "H_low",
                              nutrient == "N" ~ "H_n",
                              nutrient == "P" ~ "H_p",
                              nutrient == "N+P" ~ "H_np"),
         variable2 = "H",
         value = full_mass_g,
         time = dpp,
         nutrient = fct_recode(nutrient, "+N" = "N",
                               "+P" = "P",
                               "+N+P" = "N+P") %>%
           fct_relevel("low", "+N", "+P"))

# wrapper function
plant_wrapper <- function(m, g){
  
  out <- ode(c(R_n_low = R0_n_lo, R_p_low = R0_p_lo, Q_n_low = Q0_n, Q_p_low = Q0_p, H_low = H0,
               R_n_n = R0_n_hi, R_p_n = R0_p_lo, Q_n_n = Q0_n, Q_p_n = Q0_p, H_n = H0,
               R_n_p = R0_n_lo, R_p_p = R0_p_hi, Q_n_p = Q0_n, Q_p_p = Q0_p, H_p = H0,
               R_n_np = R0_n_hi, R_p_np = R0_p_hi, Q_n_np = Q0_n, Q_p_np = Q0_p, H_np = H0),
             times, plant_model, c(m = m, g = g)) %>%
    as_tibble() %>%
    mutate(across(everything(), as.double)) %>%
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
                                 str_starts(variable, "H") == T ~ "H"))
  
  # use to visualize all variables
  # ggplot(out, aes(x = time, y = value, color = nutrient)) +
  #   geom_line() +
  #   stat_summary(data = mock, geom = "errorbar", width = 0, fun.data = "mean_se") +
  #   stat_summary(data = mock, geom = "point", size = 2, fun = "mean") +
  #   facet_wrap(~ variable2, scales = "free")
  
  # mock_mean <- mock %>%
  #   group_by(nutrient) %>%
  #   summarize(value = mean(value)) %>%
  #   ungroup()
  # 
  # # visualize only H
  # ggplot(filter(out, variable2 == "H"), aes(x = time, y = value)) +
  #   geom_line() +
  #   geom_hline(data = mock_mean, aes(yintercept = value), linetype = "dashed") +
  #   geom_point(data = mock) +
  #   facet_wrap(~ nutrient, scales = "free")

  # figure for supplement
  ggplot(filter(out, variable2 == "H"), 
         aes(x = time, y = value, linetype = nutrient, color = nutrient)) +
    geom_line() +
    geom_point(data = mock) +
    facet_rep_wrap(~ nutrient, nrow = 1) +
    scale_color_viridis_d(end = 0.7, name = "Nutrient supply") +
    scale_linetype(name = "Nutrient supply") +
    labs(x = "Time (days)", y = "Plant biomass (g)", title = "(A)") +
    fig_theme
}

# initiate slider for ggplot
manipulate(plot(1:5, cex=size), size = slider(0.5,10,step=0.5))

# time 
times <- seq(0, max(dat5$dpp)*2, length.out = 100)

manipulate(plant_wrapper(m, g), m = slider(0, 0.1), g = slider(0, 1))
# g ~ 0.144
# m ~ 0.004

# set g so that we only estimate one parameter
g <- 0.144


#### compare plant model to observations ####

# data
mock_fit <- mock %>%
  rename("name" = "variable") %>%
  select(name, time, value) %>%
  data.frame()

# time
times <- seq(0, max(dat5$dpp), length.out = 100)

# cost function
plant_cost <- function(input_plant){ 
  m = input_plant[1];
  out = ode(y = c(R_n_low = R0_n_lo, R_p_low = R0_p_lo, Q_n_low = Q0_n, Q_p_low = Q0_p, H_low = H0,
                  R_n_n = R0_n_hi, R_p_n = R0_p_lo, Q_n_n = Q0_n, Q_p_n = Q0_p, H_n = H0,
                  R_n_p = R0_n_lo, R_p_p = R0_p_hi, Q_n_p = Q0_n, Q_p_p = Q0_p, H_p = H0,
                  R_n_np = R0_n_hi, R_p_np = R0_p_hi, Q_n_np = Q0_n, Q_p_np = Q0_p, H_np = H0),
            times = times, func = plant_model, parms = c(m = m))
  return(modCost(model = out[ , c("time", "H_low", "H_n", "H_p", "H_np")], obs = mock_fit, y = "value"))   
}


#### estimate plant parameters ####

#initial guess
input_plant <- c(0.004) 

# fit model
fit_plant <- modFit(plant_cost, input_plant, lower = c(0))
summary(fit_plant)
deviance(fit_plant)
fit_plant$ssr # sum of squared residuals
fit_plant$ms # mean squared residuals


#### visualize plant model fit ####

# save value
m <- fit_plant$par[1]; # 0.007

# time 
times <- seq(0, max(dat5$dpp), length.out = 100)

# figure
plant_fig <- plant_wrapper(m, g)


#### fit virus parameters ####


#### virus model ####

virus_model = function (t, yy, parms) { 
  
  # unknown parameters and supply rates
  c = parms[1];
  r = ifelse(length(parms) == 2, parms[2], r);
  
  # set initial values
  R_n_low = yy[1];
  R_p_low = yy[2];
  Q_n_low = yy[3];
  Q_p_low = yy[4];
  H_low = yy[5];
  V_low = yy[6];
  
  R_n_n = yy[7];
  R_p_n = yy[8];
  Q_n_n = yy[9];
  Q_p_n = yy[10];
  H_n = yy[11];
  V_n = yy[12];
  
  R_n_p = yy[13];
  R_p_p = yy[14];
  Q_n_p = yy[15];
  Q_p_p = yy[16];
  H_p = yy[17]
  V_p = yy[18];
  
  R_n_np = yy[19];
  R_p_np = yy[20];
  Q_n_np = yy[21];
  Q_p_np = yy[22];
  H_np = yy[23];
  V_np = yy[24];
  
  # model
  dR_n_low = a_n_lo - (u_n * R_n_low * H_low) / (R_n_low + k_n);
  dR_p_low = a_p_lo - (u_p * R_p_low * H_low) / (R_p_low + k_p);
  dQ_n_low = (u_n * R_n_low) / (R_n_low + k_n) - min((1 - Qmin_n / Q_n_low), (1 - Qmin_p / Q_p_low)) * g * Q_n_low - min((1 - q_n / Q_n_low), (1 - q_p / Q_p_low)) * z_n * r * V_low
  dQ_p_low = (u_p * R_p_low) / (R_p_low + k_p) - min((1 - Qmin_n / Q_n_low), (1 - Qmin_p / Q_p_low)) * g * Q_p_low - min((1 - q_n / Q_n_low), (1 - q_p / Q_p_low)) * z_p * r * V_low
  dH_low = min((1 - Qmin_n / Q_n_low), (1 - Qmin_p / Q_p_low)) * g * H_low - m * H_low
  dV_low = min((1 - q_n / Q_n_low), (1 - q_p / Q_p_low)) * r * V_low - c * V_low
  
  dR_n_n = a_n_hi - (u_n * R_n_n * H_n) / (R_n_n + k_n);
  dR_p_n = a_p_lo - (u_p * R_p_n * H_n) / (R_p_n + k_p);
  dQ_n_n = (u_n * R_n_n) / (R_n_n + k_n) - min((1 - Qmin_n / Q_n_n), (1 - Qmin_p / Q_p_n)) * g * Q_n_n - min((1 - q_n / Q_n_n), (1 - q_p / Q_p_n)) * z_n * r * V_n
  dQ_p_n = (u_p * R_p_n) / (R_p_n + k_p) - min((1 - Qmin_n / Q_n_n), (1 - Qmin_p / Q_p_n)) * g * Q_p_n - min((1 - q_n / Q_n_n), (1 - q_p / Q_p_n)) * z_p * r * V_n
  dH_n = min((1 - Qmin_n / Q_n_n), (1 - Qmin_p / Q_p_n)) * g * H_n - m * H_n
  dV_n = min((1 - q_n / Q_n_n), (1 - q_p / Q_p_n)) * r * V_n - c * V_n
  
  dR_n_p = a_n_lo - (u_n * R_n_p * H_p) / (R_n_p + k_n);
  dR_p_p = a_p_hi - (u_p * R_p_p * H_p) / (R_p_p + k_p);
  dQ_n_p = (u_n * R_n_p) / (R_n_p + k_n) - min((1 - Qmin_n / Q_n_p), (1 - Qmin_p / Q_p_p)) * g * Q_n_p - min((1 - q_n / Q_n_p), (1 - q_p / Q_p_p)) * z_n * r * V_p
  dQ_p_p = (u_p * R_p_p) / (R_p_p + k_p) - min((1 - Qmin_n / Q_n_p), (1 - Qmin_p / Q_p_p)) * g * Q_p_p - min((1 - q_n / Q_n_p), (1 - q_p / Q_p_p)) * z_p * r * V_p
  dH_p = min((1 - Qmin_n / Q_n_p), (1 - Qmin_p / Q_p_p)) * g * H_p - m * H_p
  dV_p = min((1 - q_n / Q_n_p), (1 - q_p / Q_p_p)) * r * V_p - c * V_p
  
  dR_n_np = a_n_hi - (u_n * R_n_np * H_np) / (R_n_np + k_n);
  dR_p_np = a_p_hi - (u_p * R_p_np * H_np) / (R_p_np + k_p);
  dQ_n_np = (u_n * R_n_np) / (R_n_np + k_n) - min((1 - Qmin_n / Q_n_np), (1 - Qmin_p / Q_p_np)) * g * Q_n_np - min((1 - q_n / Q_n_np), (1 - q_p / Q_p_np)) * z_n * r * V_np
  dQ_p_np = (u_p * R_p_np) / (R_p_np + k_p) - min((1 - Qmin_n / Q_n_np), (1 - Qmin_p / Q_p_np)) * g * Q_p_np - min((1 - q_n / Q_n_np), (1 - q_p / Q_p_np)) * z_p * r * V_np
  dH_np = min((1 - Qmin_n / Q_n_np), (1 - Qmin_p / Q_p_np)) * g * H_np - m * H_np
  dV_np = min((1 - q_n / Q_n_np), (1 - q_p / Q_p_np)) * r * V_np - c * V_np
  
  return(list(c(dR_n_low, dR_p_low, dQ_n_low, dQ_p_low, dH_low, dV_low,
                dR_n_n, dR_p_n, dQ_n_n, dQ_p_n, dH_n, dV_n,
                dR_n_p, dR_p_p, dQ_n_p, dQ_p_p, dH_p, dV_p,
                dR_n_np, dR_p_np, dQ_n_np, dQ_p_np, dH_np, dV_np)))
}


#### visualize virus parameters ####

# data
pav <- dat5 %>%
  filter(inoc == "PAV") %>%
  mutate(variable = case_when(nutrient == "low" ~ "V_low",
                              nutrient == "N" ~ "V_n",
                              nutrient == "P" ~ "V_p",
                              nutrient == "N+P" ~ "V_np"),
         variable2 = "V",
         value = conc_PAV,
         time = dpi,
         nutrient = fct_recode(nutrient, "+N" = "N",
                               "+P" = "P",
                               "+N+P" = "N+P") %>%
           fct_relevel("low", "+N", "+P"))

rpv <- dat5 %>%
  filter(inoc == "RPV") %>%
  mutate(variable = case_when(nutrient == "low" ~ "V_low",
                              nutrient == "N" ~ "V_n",
                              nutrient == "P" ~ "V_p",
                              nutrient == "N+P" ~ "V_np"),
         variable2 = "V",
         value = conc_RPV,
         time = dpi,
         nutrient = fct_recode(nutrient, "+N" = "N",
                               "+P" = "P",
                               "+N+P" = "N+P") %>%
           fct_relevel("low", "+N", "+P"))

# min dpp
min(pav$dpp)

# min virions
min(pav$conc_PAV)
min(rpv$conc_RPV)

# initial plant values
virus_init_fun <- function(plant_time){
  
  # simulate plant growth
  plant_init <- ode(c(R_n_low = R0_n_lo, R_p_low = R0_p_lo, Q_n_low = Q0_n, Q_p_low = Q0_p, H_low = H0,
                      R_n_n = R0_n_hi, R_p_n = R0_p_lo, Q_n_n = Q0_n, Q_p_n = Q0_p, H_n = H0,
                      R_n_p = R0_n_lo, R_p_p = R0_p_hi, Q_n_p = Q0_n, Q_p_p = Q0_p, H_p = H0,
                      R_n_np = R0_n_hi, R_p_np = R0_p_hi, Q_n_np = Q0_n, Q_p_np = Q0_p, H_np = H0),
                    times = seq(0, plant_time, length.out = 100), plant_model, c(m, g)) %>%
    as_tibble() %>%
    mutate(across(everything(), as.double)) %>%
    filter(time == plant_time)
  
  y0_virus <- c(R_n_low = pull(plant_init, R_n_low), 
                R_p_low = pull(plant_init, R_p_low), 
                Q_n_low = pull(plant_init, Q_n_low), 
                Q_p_low = pull(plant_init, Q_p_low), 
                H_low = pull(plant_init, H_low), 
                V_low = V0,
                R_n_n = pull(plant_init, R_n_n), 
                R_p_n = pull(plant_init, R_p_n), 
                Q_n_n = pull(plant_init, Q_n_n), 
                Q_p_n = pull(plant_init, Q_p_n), 
                H_n = pull(plant_init, H_n), 
                V_n = V0,
                R_n_p = pull(plant_init, R_n_p), 
                R_p_p = pull(plant_init, R_p_p), 
                Q_n_p = pull(plant_init, Q_n_p), 
                Q_p_p = pull(plant_init, Q_p_p), 
                H_p = pull(plant_init, H_p), 
                V_p = V0,
                R_n_np = pull(plant_init, R_n_np), 
                R_p_np = pull(plant_init, R_p_np), 
                Q_n_np = pull(plant_init, Q_n_np), 
                Q_p_np = pull(plant_init, Q_p_np), 
                H_np = pull(plant_init, H_np), 
                V_np = V0)
  
  return(y0_virus)
}

# wrapper function
virus_wrapper <- function(c, r, species, plant_time){
  
  out <- ode(virus_init_fun(plant_time),
             times, virus_model, c(c = c, r = r)) %>%
    as_tibble() %>%
    mutate(across(everything(), as.double)) %>%
    pivot_longer(cols = -time,
                 names_to = "variable",
                 values_to = "value") %>%
    mutate(nutrient = case_when(str_ends(variable, "low") == T ~ "low",
                                str_ends(variable, "np") == T ~ "+N+P",
                                str_ends(variable, "n") == T ~ "+N",
                                str_ends(variable, "p") == T ~ "+P"),
           variable2 = case_when(str_starts(variable, "R_n") == T ~ "R_n",
                                 str_starts(variable, "R_p") == T ~ "R_p",
                                 str_starts(variable, "Q_n") == T ~ "Q_n",
                                 str_starts(variable, "Q_p") == T ~ "Q_p",
                                 str_starts(variable, "H") == T ~ "H",
                                 str_starts(variable, "V") == T ~ "V"))
  
  if(species == "PAV"){
    virus <- pav
  }else{
    virus <- rpv
  }
  
  # use to visualize all variables
  # ggplot(out, aes(x = time, y = value, color = nutrient)) +
  #   geom_line() +
  #   stat_summary(data = virus, geom = "errorbar", width = 0, fun.data = "mean_se") +
  #   stat_summary(data = virus, geom = "point", size = 2, fun = "mean") +
  #   facet_wrap(~ variable2, scales = "free")
  
  # virus_mean <- virus %>%
  #   group_by(nutrient) %>%
  #   summarize(value = mean(value)) %>%
  #   ungroup()
  # 
  # # average raw data
  # ggplot(filter(out, variable2 == "V"), aes(x = time, y = value)) +
  #   geom_line() +
  #   geom_hline(data = virus_mean, aes(yintercept = value), linetype = "dashed") +
  #   stat_summary(data = virus, geom = "errorbar", width = 0, fun.data = "mean_se") +
  #   stat_summary(data = virus, geom = "point", size = 2, fun = "mean") +
  #   facet_wrap(~ nutrient, scales = "free")
  
  # raw data points
  ggplot(filter(out, variable2 == "V"), aes(x = time, y = value)) +
    geom_line() +
    # geom_point(data = virus) +
    facet_wrap(~ nutrient, scales = "free")
}

# initiate slider for ggplot
manipulate(plot(1:5, cex=size), size = slider(0.5,10,step=0.5))

# time 
times <- seq(0, max(dat5$dpi), length.out = 100)

# PAV
z_n <- z_nb
z_p <- z_pb
V0 <- V0_init * 5
manipulate(virus_wrapper(c, r, species = "PAV", plant_time = plant_days), c = slider(0, 0.1), r = slider(0, 1))
# r ~ 0.092
# c ~ 0.0056

# RPV
z_n <- z_nc
z_p <- z_pc
V0 <- V0_init * 5
manipulate(virus_wrapper(c, r, species = "RPV", plant_time = plant_days), c = slider(0, 0.5), r = slider(0, 1))
# r ~ 0.27
# c ~ 0.005

# set c so that we only estimate one parameter
r_b <- 0.09
r_c <- r_b * 3


#### compare virus model to observations ####

# data
pav_fit <- pav %>%
  rename("name" = "variable") %>%
  select(name, time, value) %>%
  data.frame()

rpv_fit <- rpv %>%
  rename("name" = "variable") %>%
  select(name, time, value) %>%
  data.frame()

# cost functions
pav_cost <- function(input_virus){ 
  c = input_virus[1];
  out = ode(y = virus_init_fun(plant_days), 
            times = times_pav, func = virus_model, parms = c(c = c))
  return(modCost(model = out[ , c("time", "V_low", "V_n", "V_p", "V_np")], obs = pav_fit, y = "value"))   
}

rpv_cost <- function(input_virus){ 
  c = input_virus[1];
  out = ode(y = virus_init_fun(plant_days), 
            times = times_rpv, func = virus_model, parms = c(c = c))
  return(modCost(model = out[ , c("time", "V_low", "V_n", "V_p", "V_np")], obs = rpv_fit, y = "value"))   
}


#### estimate virus parameters ####

# times
times_pav <- seq(0, max(pav_fit$time), length.out = 100)
times_rpv <- seq(0, max(rpv_fit$time), length.out = 100)

#initial guess
input_pav <- c(0.0055)
input_rpv <- c(0.005)

# fit PAV model
z_n <- z_nb
z_p <- z_pb
r <- r_b
fit_pav <- modFit(pav_cost, input_pav, lower = c(0))
summary(fit_pav)
deviance(fit_pav)
fit_pav$ssr # sum of squared residuals
fit_pav$ms # mean squared residuals

# fit RPV model
z_n <- z_nc
z_p <- z_pc
r <- r_c
fit_rpv <- modFit(rpv_cost, input_rpv, lower = c(0))
summary(fit_rpv)
deviance(fit_rpv)
fit_rpv$ssr # sum of squared residuals
fit_rpv$ms # mean squared residuals


#### visualize virus model fit ####

# fit PAV model
c <- c_b <- fit_pav$par[1]
r <- r_b
z_n <- z_nb
z_p <- z_pb
times <- seq(0, max(dat5$dpi), length.out = 100)

# quick figure
virus_wrapper(c, r, species = "PAV", plant_time = (plant_days + 12)) +
  fig_theme

# pred_pav <- ode(y = virus_init_fun(plant_days),
#                 times = times_pav, func = virus_model, parms = c(c_b)) %>%
#   as_tibble() %>%
#   mutate(across(everything(), as.double))
# 
# # add time to data
# pav_fit2 <- pav_fit %>%
#   mutate(time = time + plant_days)
# 
# # visualize
# ggplot(pav_fit2, aes(time, value)) +
#   stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
#   stat_summary(geom = "point", fun = "mean") +
#   geom_line(data = pred_pav, aes(y = V)) +
#   labs(x = "Time (days)", y = expression(paste("BYDV-PAV conc. (", g^-1, ")", sep = "")), title = "(B) BYDV-PAV growth rate") +
#   fig_theme +
#   theme(plot.title = element_text(size = 9, vjust = 0, hjust = 0.5))

# fit RPV model
c <- c_c <- fit_rpv$par[1]
r <- r_c
z_n <- z_nc
z_p <- z_pc

# quick figure
virus_wrapper(c, r, species = "RPV", plant_time = (plant_days + 12)) +
  fig_theme

# pred_rpv <- ode(y = virus_init_fun("high", "high"), 
#                 times = times_rpv, func = virus_model, parms = c(r_c)) %>%
#   as_tibble() %>%
#   mutate(across(everything(), as.double),
#          time = time + 11)
# 
# # add time to data
# rpv_NP2 <- rpv %>%
#   filter(nutrient == "N+P") %>%
#   rename("name" = "variable") %>%
#   select(name, time, value) %>%
#   data.frame() %>%
#   mutate(time = time + 11)
# 
# # visualize
# rpv_fig <- ggplot(rpv_NP2, aes(time, value)) +
#   stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
#   stat_summary(geom = "point", fun = "mean") +
#   geom_line(data = pred_rpv, aes(y = V)) +
#   labs(x = "Time (days)", y = expression(paste("CYDV-RPV conc. (", g^-1, ")", sep = "")), title = "(C) CYDV-RPV growth rate") +
#   fig_theme +
#   theme(plot.title = element_text(size = 9, vjust = 0, hjust = 0.5))


#### output ####

# figure
pdf("output/growth_rate_fit_figure.pdf", width = 6.5, height = 2)
plot_grid(plant_fig, pav_fig, rpv_fig,
          nrow = 1)
dev.off()

# parameters
write_csv(tibble(parameter = c("g", "m", "r_b", "c_b", "r_c", "c_c"),
                 value = c(g, m, r_b, c_b, r_c, c_c)),
          "output/estimated_growth_mortality_rates.csv")

