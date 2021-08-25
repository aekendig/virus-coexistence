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

# import data
sdat <- read_csv("./edi.411.2/data/sample_exp_molc_data.csv")


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
              names_glue = "quant_{target}") %>%
  filter(!(sample %in% accs) & 
           !(inoc %in% c("PAV", "coinfection") & quant_PAV == 0) & 
           !(inoc %in% c("RPV", "coinfection") & quant_RPV == 0)) %>%
  mutate(conc_PAV = quant_PAV * 50 / mass_ext_mg,
         conc_RPV = quant_RPV * 50 / mass_ext_mg,
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
R0_n <- a_n_lo * 7
R0_p <- a_p_lo * 7
Q0_n <- Qmin_n
Q0_p <- Qmin_p
B0 <- 1e-3

# Combine values for model
y0_plant <- c(R_n = R0_n, R_p = R0_p, Q_n = Q0_n, Q_p = Q0_p, B = B0)


#### plant model ####

plant_model = function (t, yy, parms) { 
  
  # supply rates
  g = parms[1];
  a_n = ifelse(length(parms) > 1, parms[2], a_n);
  a_p = ifelse(length(parms) > 2, parms[3], a_p);
  m = ifelse(length(parms) == 4, parms[4], m);

  # set initial values
  R_n = yy[1];
  R_p = yy[2];
  Q_n = yy[3];
  Q_p = yy[4];
  B = yy[5];
 
  # model
  dR_n = a_n - (u_n * R_n * B) / (R_n + k_n);
  dR_p = a_p - (u_p * R_p * B) / (R_p + k_p);
  dQ_n = (u_n * R_n) / (R_n + k_n) - min((1 - Qmin_n / Q_n), (1 - Qmin_p / Q_p)) * g * Q_n
  dQ_p = (u_p * R_p) / (R_p + k_p) - min((1 - Qmin_n / Q_n), (1 - Qmin_p / Q_p)) * g * Q_p
  dB = min((1 - Qmin_n / Q_n), (1 - Qmin_p / Q_p)) * g * B - m * B
  
  return(list(c(dR_n, dR_p, dQ_n, dQ_p, dB)))
}


#### visualize plant parameters ####

# data
mock <- dat5 %>%
  filter(inoc == "healthy") %>%
  mutate(variable = "B",
         value = full_mass_g,
         time = dpp)

# initiate slider for ggplot
manipulate(plot(1:5, cex=size), size = slider(0.5,10,step=0.5))

# time 
times <- seq(0, max(dat5$dpp), length.out = 100)
# times <- seq(0, 100)

# wrapper function
plant_wrapper <- function(g, m){
  
  out_low <- ode(y0_plant, times, plant_model, c(g = g, a_n = a_n_lo, a_p = a_p_lo, m = m)) %>%
    as_tibble() %>%
    mutate(nutrient = "low",
           across(!nutrient, as.double),
           nutrient = as.character(nutrient))
  
  out_N <- ode(y0_plant, times, plant_model, c(g = g, a_n = a_n_hi, a_p = a_p_lo, m = m)) %>%
    as_tibble() %>%
    mutate(nutrient = "N",
           across(!nutrient, as.double),
           nutrient = as.character(nutrient))
  
  out_P <- ode(y0_plant, times, plant_model, c(g = g, a_n = a_n_lo, a_p = a_p_hi, m = m)) %>%
    as_tibble() %>%
    mutate(nutrient = "P",
           across(!nutrient, as.double),
           nutrient = as.character(nutrient))
  
  out_NP <- ode(y0_plant, times, plant_model, c(g = g, a_n = a_n_hi, a_p = a_p_hi, m = m)) %>%
    as_tibble() %>%
    mutate(nutrient = "N+P",
           across(!nutrient, as.double),
           nutrient = as.character(nutrient))
  
  out <- out_low %>%
    full_join(out_N) %>%
    full_join(out_P) %>%
    full_join(out_NP) %>%
    pivot_longer(cols = c(R_n, R_p, Q_n, Q_p, B),
                 names_to = "variable",
                 values_to = "value")
  
  ggplot(out, aes(x = time, y = value, color = nutrient)) +
    geom_line() +
    stat_summary(data = mock, geom = "errorbar", width = 0, fun.data = "mean_se") +
    stat_summary(data = mock, geom = "point", size = 2, fun = "mean") +
    facet_wrap(~ variable, scales = "free")
}

manipulate(plant_wrapper(g, m), g = slider(0.001, 1), m = slider(0.001, 1))
# only fits the high N+P data
# m ~ 0.03
# g ~ 0.3
# predicted biomass of other treatments way lower and don't increase with parameter adjustments

# set m so that we only estimate one parameter
m <- 0.03


#### compare plant model to observations ####

# data
mock_NP <- mock %>%
  filter(nutrient == "N+P") %>%
  rename("name" = "variable") %>%
  select(name, time, value) %>%
  data.frame()

# nutrient supply rates
a_n <- a_n_hi
a_p <- a_p_hi

# # Rename columns 
# colnames(rdatCH)=c("name","time","val")
# 
# # Combine data
# rdat3=subset(rdatCH,name=="Hh"|name=="Ph")

plant_cost <- function(input_plant){ 
  g = input_plant[1];
  out = ode(y = y0_plant, times = seq(0, max(mock_NP$time), length.out = 100), func = plant_model, parms = c(g = g))
  return(modCost(model = out[ , c("time", "B")], obs = mock_NP, y = "value"))   
  # return(out[ , c("time", "B")]) #for troubleshooting
}


#### estimate plant parameters ####

#initial guess
input_plant <- c(0.3) 

# fit model
fit_plant <- modFit(plant_cost, input_plant, lower = c(0))
summary(fit_plant)
deviance(fit_plant)
fit_plant$ssr # sum of squared residuals
fit_plant$ms # mean squared residuals


#### visualize plant model fit ####

# save value
g <- fit_plant$par[1]; # 0.27

# fit model
pred_plant <- ode(y = y0_plant, times = seq(0, max(mock_NP$time), length.out = 100), func = plant_model, parms = c(g)) %>%
  as_tibble() %>%
  mutate(across(everything(), as.double))

# visualize
ggplot(mock_NP, aes(time, value)) +
  geom_point() +
  geom_line(data = pred_plant, aes(y = B))


#### virus model ####

virus_model = function (t, yy, parms) { 
  
  # supply rates
  r = parms[1];
  a_n = ifelse(length(parms) > 1, parms[2], a_n);
  a_p = ifelse(length(parms) > 2, parms[3], a_p);
  c = ifelse(length(parms) == 4, parms[4], c);
  
  # set initial values
  R_n = yy[1];
  R_p = yy[2];
  Q_n = yy[3];
  Q_p = yy[4];
  B = yy[5];
  V = yy[6];
  
  # model
  dR_n = a_n - (u_n * R_n * B) / (R_n + k_n);
  dR_p = a_p - (u_p * R_p * B) / (R_p + k_p);
  dQ_n = (u_n * R_n) / (R_n + k_n) - min((1 - Qmin_n / Q_n), (1 - Qmin_p / Q_p)) * g * Q_n - min((1 - q_n / Q_n), (1 - q_p / Q_p)) * z_n * r * V
  dQ_p = (u_p * R_p) / (R_p + k_p) - min((1 - Qmin_n / Q_n), (1 - Qmin_p / Q_p)) * g * Q_p - min((1 - q_n / Q_n), (1 - q_p / Q_p)) * z_p * r * V
  dB = min((1 - Qmin_n / Q_n), (1 - Qmin_p / Q_p)) * g * B - m * B
  dV = min((1 - q_n / Q_n), (1 - q_p / Q_p)) * r * V - c * V
  
  return(list(c(dR_n, dR_p, dQ_n, dQ_p, dB, dV)))
}


#### visualize virus parameters ####

# data
pav <- dat5 %>%
  filter(inoc == "PAV") %>%
  mutate(variable = "V",
         value = conc_PAV,
         time = dpi)

#### start here: predict plant size at 0 dpi and change starting condition vector ####
# all code below this is for plant fitting, not virus

B0 <-
V0 <- 1e-3
y0_virus <- c(R_n = R0_n, R_p = R0_p, Q_n = Q0_n, Q_p = Q0_p, B = B0, V = V0)

# initiate slider for ggplot
manipulate(plot(1:5, cex=size), size = slider(0.5,10,step=0.5))

# time 
# times <- seq(0, max(dat5$dpi), length.out = 100)
times <- seq(0, 100)

# wrapper function
plant_wrapper <- function(g, m){
  
  out_low <- ode(y0_plant, times, plant_model, c(g = g, a_n = a_n_lo, a_p = a_p_lo, m = m)) %>%
    as_tibble() %>%
    mutate(nutrient = "low",
           across(!nutrient, as.double),
           nutrient = as.character(nutrient))
  
  out_N <- ode(y0_plant, times, plant_model, c(g = g, a_n = a_n_hi, a_p = a_p_lo, m = m)) %>%
    as_tibble() %>%
    mutate(nutrient = "N",
           across(!nutrient, as.double),
           nutrient = as.character(nutrient))
  
  out_P <- ode(y0_plant, times, plant_model, c(g = g, a_n = a_n_lo, a_p = a_p_hi, m = m)) %>%
    as_tibble() %>%
    mutate(nutrient = "P",
           across(!nutrient, as.double),
           nutrient = as.character(nutrient))
  
  out_NP <- ode(y0_plant, times, plant_model, c(g = g, a_n = a_n_hi, a_p = a_p_hi, m = m)) %>%
    as_tibble() %>%
    mutate(nutrient = "N+P",
           across(!nutrient, as.double),
           nutrient = as.character(nutrient))
  
  out <- out_low %>%
    full_join(out_N) %>%
    full_join(out_P) %>%
    full_join(out_NP) %>%
    pivot_longer(cols = c(R_n, R_p, Q_n, Q_p, B),
                 names_to = "variable",
                 values_to = "value")
  
  ggplot(out, aes(x = time, y = value, color = nutrient)) +
    geom_line() +
    stat_summary(data = mock, geom = "errorbar", width = 0, fun.data = "mean_se") +
    stat_summary(data = mock, geom = "point", size = 2, fun = "mean") +
    facet_wrap(~ variable, scales = "free")
}

manipulate(plant_wrapper(g, m), g = slider(0.001, 1), m = slider(0.001, 1))
# only fits the high N+P data
# m ~ 0.03
# g ~ 0.7
# predicted biomass of other treatments way lower and don't increase with parameter adjustments

# set m so that we only estimate one parameter
m <- 0.03


#### compare plant model to observations ####

# data
mock_NP <- mock %>%
  filter(nutrient == "N+P") %>%
  rename("name" = "variable") %>%
  select(name, time, value) %>%
  data.frame()

# nutrient supply rates
a_n <- a_n_hi
a_p <- a_p_hi

# # Rename columns 
# colnames(rdatCH)=c("name","time","val")
# 
# # Combine data
# rdat3=subset(rdatCH,name=="Hh"|name=="Ph")

plant_cost <- function(input_plant){ 
  g = input_plant[1];
  out = ode(y = y0_plant, times = seq(0, max(mock_NP$time), length.out = 100), func = plant_model, parms = c(g = g))
  return(modCost(model = out[ , c("time", "B")], obs = mock_NP, y = "value"))   
  # return(out[ , c("time", "B")]) #for troubleshooting
}


#### estimate plant parameters ####

#initial guess
input_plant <- c(0.7) 

# fit model
fit_plant <- modFit(plant_cost, input_plant, lower = c(0))
summary(fit_plant)
deviance(fit_plant)
fit_plant$ssr # sum of squared residuals
fit_plant$ms # mean squared residuals


#### visualize plant model fit ####

# fit model
g <- fit_plant$par[1];
pred_plant <- ode(y = y0_plant, times = seq(0, max(mock_NP$time), length.out = 100), func = plant_model, parms = c(g)) %>%
  as_tibble() %>%
  mutate(across(everything(), as.double))

# visualize
ggplot(mock_NP, aes(time, value)) +
  geom_point() +
  geom_line(data = pred_plant, aes(y = B))