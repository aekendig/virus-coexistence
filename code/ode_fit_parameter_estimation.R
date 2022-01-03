## Goal: estimate plant and virus growth rates

#### start here ####
# ode_simulations and model_settings are all updated
# current virus parameters make PAV grow linearly without limit and RPV crashes by ~ 600 days


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


#### edit data ####
# code source: concentration_analysis.R in edi folder

# days post inoculation
dpi <- tibble(time = 1:8,
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


#### parameters and model functions ####

# load model parameters and figure settings
source("code/model_settings.R")

# time
plant_days <- 11


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

mock_mean <- mock %>%
  group_by(nutrient) %>%
  summarize(value = mean(value)) %>%
  ungroup()

# wrapper function
plant_wrapper <- function(m, g){
  
  # set to default
  params_in <- params_def1
  
  # update parameter value
  params_in["m"] <- m
  params_in["g"] <- g
  
  # fit model
  out <- ode(init_def1, times, plant_model, params_in) %>%
    plant_model_format()

  # # use to visualize all variables
  # ggplot(out, aes(x = time, y = value, color = nutrient)) +
  #   geom_line() +
  #   stat_summary(data = mock, geom = "errorbar", width = 0, fun.data = "mean_se") +
  #   stat_summary(data = mock, geom = "point", size = 2, fun = "mean") +
  #   facet_wrap(~ variable2, scales = "free")
  
  # visualize only H
  ggplot(filter(out, variable2 == "H"), aes(x = time, y = value)) +
    geom_line() +
    geom_hline(data = mock_mean, aes(yintercept = value), linetype = "dashed") +
    geom_point(data = mock) +
    facet_wrap(~ nutrient, scales = "free")
}

# initiate slider for ggplot
manipulate(plot(1:5, cex=size), size = slider(0.5,10,step=0.5))

# time 
times <- seq(0, max(dat5$dpp)*2, length.out = 100)

manipulate(plant_wrapper(m, g), m = slider(0, 0.1), g = slider(0, 1))
# g ~ 0.144
# m ~ 0.004

# set g so that we only estimate one parameter
params_def1["g"] <- 0.144


#### compare plant model to observations ####

# data
mock_fit <- mock %>%
  rename("name" = "variable") %>%
  select(name, time, value) %>%
  data.frame()

# time
times <- seq(0, max(dat5$dpp), length.out = 100)

# cost function
plant_cost <- function(m_init){ 
  
  # set to default
  params_in <- params_def1
  
  # update parameter value
  params_in["m"] <- m_init
  
  # fit model
  out = ode(y = init_def1, times = times, func = plant_model, parms = params_in)
  
  # compare to data
  return(modCost(model = out[ , c("time", "H_low", "H_n", "H_p", "H_np")], obs = mock_fit, y = "value"))   
}


#### estimate plant parameters ####

# fit model
fit_plant <- modFit(plant_cost, 0.004, lower = c(0))
summary(fit_plant)
deviance(fit_plant)
fit_plant$ssr # sum of squared residuals
fit_plant$ms # mean squared residuals


#### visualize plant model fit ####

# save value
params_def1["m"] <- fit_plant$par[1]

# time 
times <- seq(0, max(dat5$dpp), length.out = 100)

# data for plant figure
plant_pred_dat <- ode(init_def1, times, plant_model, params_def1) %>%
  plant_model_format() %>%
  filter(variable2 == "H")

# figure
plant_pred_dat %>%
  ggplot(aes(x = time, y = value)) +
  geom_line() +
  geom_point(data = mock) +
  facet_wrap(~ nutrient, scales = "free")


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

# initial values
init_virus1 <- virus1_model_init()

# wrapper function
virus_wrapper <- function(c, r, species){
  
  if(species == "PAV"){
    virus <- pav
  }else{
    virus <- rpv
  }
  
  # set to default
  params_in <- params_def1
  
  # update parameter value
  params_in["c"] <- c
  params_in["r"] <- r
  
  # fit model
  out <- ode(init_virus1, times, virus1_model, params_in) %>%
    virus1_model_format()
  
  # # visualize all variables
  # ggplot(out, aes(x = time, y = value, color = nutrient)) +
  #   geom_line() +
  #   stat_summary(data = virus, geom = "errorbar", width = 0, fun.data = "mean_se") +
  #   stat_summary(data = virus, geom = "point", size = 2, fun = "mean") +
  #   facet_wrap(~ variable2, scales = "free")
  
  # average raw data
  virus_mean <- virus %>%
    group_by(nutrient) %>%
    summarize(value = mean(value)) %>%
    ungroup()

  ggplot(filter(out, variable2 == "V"), aes(x = time, y = value)) +
    geom_line() +
    geom_hline(data = virus_mean, aes(yintercept = value), linetype = "dashed") +
    stat_summary(data = virus, geom = "errorbar", width = 0, fun.data = "mean_se") +
    stat_summary(data = virus, geom = "point", size = 2, fun = "mean") +
    facet_wrap(~ nutrient, scales = "free")
}

# initiate slider for ggplot
manipulate(plot(1:5, cex=size), size = slider(0.5,10,step=0.5))

# time 
times <- seq(0, max(dat5$dpi), length.out = 100)

# PAV
manipulate(virus_wrapper(c, r, species = "PAV"), c = slider(0, 0.1), r = slider(0, 1))
# r ~ 0.092
# c ~ 0.0056

# RPV
manipulate(virus_wrapper(c, r, species = "RPV"), c = slider(0, 0.5), r = slider(0, 1))
# r ~ 0.27
# c ~ 0.005

# set r so that we only estimate one parameter
params_pav <- params_def1
params_rpv <- params_def1

params_pav["r"] <- 0.09
params_rpv["r"] <- 0.27


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
pav_cost <- function(c_init){ 
  
  # set to default
  params_in <- params_pav
  
  # update parameter value
  params_in["c"] <- c_init
  
  # fit model
  out = ode(y = init_virus1, times = times_pav, func = virus1_model, parms = params_in)
  
  # compare to data
  return(modCost(model = out[ , c("time", "V_low", "V_n", "V_p", "V_np")], obs = pav_fit, y = "value"))   
}

rpv_cost <- function(c_init){ 
  
  # set to default
  params_in <- params_rpv
  
  # update parameter value
  params_in["c"] <- c_init
  
  # fit model
  out = ode(y = init_virus1, times = times_rpv, func = virus1_model, parms = params_in)
  
  # compare to data
  return(modCost(model = out[ , c("time", "V_low", "V_n", "V_p", "V_np")], obs = rpv_fit, y = "value"))   
}


#### estimate virus parameters ####

# times
times_pav <- seq(0, max(pav_fit$time), length.out = 100)
times_rpv <- seq(0, max(rpv_fit$time), length.out = 100)

# fit PAV model
fit_pav <- modFit(pav_cost, 0.0055, lower = c(0))
summary(fit_pav)
deviance(fit_pav)
fit_pav$ssr # sum of squared residuals
fit_pav$ms # mean squared residuals

# fit RPV model
fit_rpv <- modFit(rpv_cost, 0.005, lower = c(0))
summary(fit_rpv)
deviance(fit_rpv)
fit_rpv$ssr # sum of squared residuals
fit_rpv$ms # mean squared residuals


#### visualize virus model fit ####

# virus times
times <- seq(0, max(dat5$dpi), length.out = 100)

# add new values
params_pav["c"] <- fit_pav$par[1]
params_rpv["c"] <- fit_rpv$par[1]

# data for figure
pav_pred_dat <- ode(init_virus1, times, virus1_model, params_pav) %>%
  virus1_model_format() %>%
  filter(variable2 == "V")

rpv_pred_dat <- ode(init_virus1, times, virus1_model, params_rpv) %>%
  virus1_model_format() %>%
  filter(variable2 == "V")

# figure
pav_pred_dat %>%
  ggplot(aes(x = time, y = value)) +
  geom_line() +
  geom_point(data = pav) +
  facet_wrap(~ nutrient, scales = "free")

rpv_pred_dat %>%
  ggplot(aes(x = time, y = value)) +
  geom_line() +
  geom_point(data = pav) +
  facet_wrap(~ nutrient, scales = "free")


#### output ####

fig_pred_dat <- plant_pred_dat %>%
  mutate(fitType = "Plant~biomass~(g)") %>%
  full_join(pav_pred_dat %>%
              mutate(fitType = "BYDV-PAV~conc.~(g^-1)",
                     time = time + plant_days)) %>%
  full_join(rpv_pred_dat %>%
              mutate(fitType = "CYDV-RPV~conc.~(g^-1)",
                     time = time + plant_days)) %>%
  mutate(fitType = fct_relevel(fitType, "Plant~biomass~(g)"))

fig_raw_dat <- mock_fit %>%
  mutate(fitType = "Plant~biomass~(g)") %>%
  full_join(pav_fit %>%
              mutate(fitType = "BYDV-PAV~conc.~(g^-1)",
                     time = time + plant_days)) %>%
  full_join(rpv_fit %>%
              mutate(fitType = "CYDV-RPV~conc.~(g^-1)",
                     time = time + plant_days)) %>%
  mutate(nutrient = case_when(str_ends(name, "low") == T ~ "low",
                              str_ends(name, "np") == T ~ "+N+P",
                              str_ends(name, "n") == T ~ "+N",
                              str_ends(name, "p") == T ~ "+P") %>%
           fct_relevel("low", "+N", "+P"),
         fitType = fct_relevel(fitType, "Plant~biomass~(g)"))

pdf("output/growth_rate_fit_figure.pdf", width = 6.5, height = 5)
ggplot(fig_pred_dat, aes(x = time, y = value, color = nutrient)) +
  geom_line(size = 1.2) +
  geom_point(data = fig_raw_dat, alpha = 0.5) +
  facet_rep_grid(fitType ~ nutrient, labeller = label_parsed,
                 scales = "free_y", switch = "y") +
  scale_color_viridis_d(end = 0.7, name = "Nutrient supply") +
  labs(x = "Time (days)") +
  fig_theme +
  theme(axis.title.y = element_blank(),
        legend.position = "none")
dev.off()

# parameters
write_csv(tibble(parameter = c("g", "m", "r_b", "c_b", "r_c", "c_c"),
                 value = c(as.numeric(params_def1["g"]),
                           as.numeric(params_def1["m"]),
                           as.numeric(params_pav["r"]),
                           as.numeric(params_pav["c"]),
                           as.numeric(params_rpv["r"]),
                           as.numeric(params_rpv["c"]))),
          "output/estimated_growth_mortality_rates.csv")

