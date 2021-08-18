## Goal: figure comparing model estimates for informative and uninformative priors from concentration_analysis.R

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse) # version used: 1.2.1
library(cowplot) # version used: 0.9.4
library(brms)  # version used: 2.7.0
library(tidybayes)  # version used: 1.0.4

# load models
load("./output/concentration_analysis_informative_pav.rda")
load("./output/concentration_analysis_informative_rpv.rda")
load("./output/concentration_analysis_uninformative_pav.rda")
load("./output/concentration_analysis_uninformative_rpv.rda")

load("./output/transmission_pav_up_concentration_interaction_uninformative.rda")
load("./output/transmission_rpv_up_concentration_interaction_uninformative.rda")
load("./output/transmission_pav_up_concentration_interaction_informative.rda")
load("./output/transmission_rpv_up_concentration_interaction_informative.rda")

# color palette
col_pal = c("black", "darkgoldenrod2", "dodgerblue1", "palegreen4")

# text sizes
sm_txt = 6
lg_txt = 8


#### edit data ####

# posterior samples
postcru <- posterior_samples(m.lu.r)
postcpu <- posterior_samples(m.lu.p)
postcri <- posterior_samples(m.li.r)
postcpi <- posterior_samples(m.li.p)

posttru <- posterior_samples(mrciu)
posttpu <- posterior_samples(mpciu)
posttri <- posterior_samples(mrcii)
posttpi <- posterior_samples(mpcii)


# rename columns
colnames(postcru) <- colnames(postcpu) <- colnames(postcri) <- colnames(postcpi) <- c("int", "co", "N", "P", "co_N", "co_P", "NP", "co_NP", "ar", "sigma", "lp")

colnames(posttru) <- colnames(posttpu) <- colnames(posttri) <- colnames(posttpi) <- c("int", "conc", "N", "P", "N_t", "P_t", "co", "NP", "NP_t", "conc_N", "conc_P", "conc_N_t", "conc_P_t", "co_N", "co_P", "co_N_t", "co_P_t", "conc_NP", "conc_NP_t", "co_NP", "co_NP_t", "sd_round", "sd_time", "round_1", "round_2", "round_3", "round_4", "time_1", "time_2", "time_3", "time_4", "time_5", "time_6", "time_7", "time_8", "lp")

# category average

avgcpu <- postcpu %>%
  transmute(low = int,
            high_N = int + N,
            high_P = int + P,
            high_NP = int + N + P + NP,
            low_co =  int + co,
            N_co = int + N + co + co_N,
            P_co = int + P + co + co_P,
            NP_co = int + N + P + co + co_N + co_P + NP + co_NP) %>%
  gather(key = "treatment", value = "effect") %>%
  mutate(Inoculation = ifelse(grepl("co", treatment, fixed = T), "coinfection", "single"),
         Inoculation = factor(Inoculation, levels = c("single", "coinfection")),
         Nutrient = recode(treatment, high_N = "N", high_P = "P", high_NP = "N+P", low_co = "low", N_co = "N", P_co = "P", NP_co = "N+P"),
         Nutrient = factor(Nutrient, levels = c("low", "N", "P", "N+P"))) %>%
  as_tibble()

avgcpi <- postcpi %>%
  transmute(low = int,
            high_N = int + N,
            high_P = int + P,
            high_NP = int + N + P + NP,
            low_co =  int + co,
            N_co = int + N + co + co_N,
            P_co = int + P + co + co_P,
            NP_co = int + N + P + co + co_N + co_P + NP + co_NP) %>%
  gather(key = "treatment", value = "effect") %>%
  mutate(Inoculation = ifelse(grepl("co", treatment, fixed = T), "coinfection", "single"),
         Inoculation = factor(Inoculation, levels = c("single", "coinfection")),
         Nutrient = recode(treatment, high_N = "N", high_P = "P", high_NP = "N+P", low_co = "low", N_co = "N", P_co = "P", NP_co = "N+P"),
         Nutrient = factor(Nutrient, levels = c("low", "N", "P", "N+P"))) %>%
  as_tibble()

avgcru <- postcru %>%
  transmute(low = int,
            high_N = int + N,
            high_P = int + P,
            high_NP = int + N + P + NP,
            low_co =  int + co,
            N_co = int + N + co + co_N,
            P_co = int + P + co + co_P,
            NP_co = int + N + P + co + co_N + co_P + NP + co_NP) %>%
  gather(key = "treatment", value = "effect") %>%
  mutate(Inoculation = ifelse(grepl("co", treatment, fixed = T), "coinfection", "single"),
         Inoculation = factor(Inoculation, levels = c("single", "coinfection")),
         Nutrient = recode(treatment, high_N = "N", high_P = "P", high_NP = "N+P", low_co = "low", N_co = "N", P_co = "P", NP_co = "N+P"),
         Nutrient = factor(Nutrient, levels = c("low", "N", "P", "N+P"))) %>%
  as_tibble()

avgcri <- postcri %>%
  transmute(low = int,
            high_N = int + N,
            high_P = int + P,
            high_NP = int + N + P + NP,
            low_co =  int + co,
            N_co = int + N + co + co_N,
            P_co = int + P + co + co_P,
            NP_co = int + N + P + co + co_N + co_P + NP + co_NP) %>%
  gather(key = "treatment", value = "effect") %>%
  mutate(Inoculation = ifelse(grepl("co", treatment, fixed = T), "coinfection", "single"),
         Inoculation = factor(Inoculation, levels = c("single", "coinfection")),
         Nutrient = recode(treatment, high_N = "N", high_P = "P", high_NP = "N+P", low_co = "low", N_co = "N", P_co = "P", NP_co = "N+P"),
         Nutrient = factor(Nutrient, levels = c("low", "N", "P", "N+P"))) %>%
  as_tibble()

combpu <- posttpu %>%
  transmute(l_l = exp(int)/(1 + exp(int)),
            l_l_co = exp(int + co)/(1 + exp(int + co)),
            l_n = exp(int + N_t)/(1 + exp(int + N_t)),
            l_n_co = exp(int + N_t + co + co_N_t)/(1 + exp(int + N_t + co + co_N_t)),
            l_p = exp(int + P_t)/(1 + exp(int + P_t)),
            l_p_co = exp(int + P_t + co + co_P_t)/(1 + exp(int + P_t + co + co_P_t)),
            l_b = exp(int + N_t+ P_t + NP_t)/(1 + exp(int + N_t+ P_t + NP_t)),
            l_b_co = exp(int + N_t+ P_t + co + NP_t + co_N_t + co_P_t + co_NP_t)/(1 + exp(int + N_t+ P_t + co + NP_t + co_N_t + co_P_t + co_NP_t)),
            n_l = exp(int + N)/(1 + exp(int + N)),
            n_l_co = exp(int + N + co + co_N)/(1 + exp(int + N + co + co_N)),
            n_n = exp(int + N + N_t)/(1 + exp(int + N + N_t)),
            n_n_co = exp(int + N + N_t + co + co_N + co_N_t)/(1 + exp(int + N + N_t + co + co_N + co_N_t)),
            n_p = exp(int + N + P_t)/(1 + exp(int + N + P_t)),
            n_p_co = exp(int + N + P_t + co + co_N + co_P_t)/(1 + exp(int + N + P_t + co + co_N + co_P_t)),
            n_b = exp(int + N + N_t+ P_t + NP_t)/(1 + exp(int + N + N_t+ P_t + NP_t)),
            n_b_co = exp(int + N + N_t+ P_t + co + co_N + NP_t + co_N_t + co_P_t + co_NP_t)/(1 + exp(int + N + N_t+ P_t + co + co_N + NP_t + co_N_t + co_P_t + co_NP_t)),
            p_l = exp(int + P)/(1 + exp(int + P)),
            p_l_co = exp(int + P + co + co_P)/(1 + exp(int + P + co + co_P)),
            p_n = exp(int + P + N_t)/(1 + exp(int + P + N_t)),
            p_n_co = exp(int + P + N_t + co + co_P + co_N_t)/(1 + exp(int + P + N_t + co + co_P + co_N_t)),
            p_p = exp(int + P + P_t)/(1 + exp(int + P + P_t)),
            p_p_co = exp(int + P + P_t + co + co_P + co_P_t)/(1 + exp(int + P + P_t + co + co_P + co_P_t)),
            p_b = exp(int + P + N_t+ P_t + NP_t)/(1 + exp(int + P + N_t+ P_t + NP_t)),
            p_b_co = exp(int + P + N_t+ P_t + co + co_P + NP_t + co_N_t + co_P_t + co_NP_t)/(1 + exp(int + P + N_t+ P_t + co + co_P + NP_t + co_N_t + co_P_t + co_NP_t)),
            b_l = exp(int + N + P + NP)/(1 + exp(int + N + P + NP)),
            b_l_co = exp(int + N + P + co + NP + co_N + co_P + co_NP)/(1 + exp(int + N + P + co + NP + co_N + co_P + co_NP)),
            b_n = exp(int + N + P + NP + N_t)/(1 + exp(int + N + P + NP + N_t)),
            b_n_co = exp(int + N + P + co + NP + co_N + co_P + co_NP + N_t + co_N_t)/(1 + exp(int + N + P + co + NP + co_N + co_P + co_NP + N_t + co_N_t)),
            b_p = exp(int + N + P + NP + P_t)/(1 + exp(int + N + P + NP + P_t)),
            b_p_co = exp(int + N + P + co + NP + co_N + co_P + co_NP + P_t + co_P_t)/(1 + exp(int + N + P + co + NP + co_N + co_P + co_NP + P_t + co_P_t)),
            b_b = exp(int + N + P + NP + N_t+ P_t + NP_t)/(1 + exp(int + N + P + NP + N_t+ P_t + NP_t)),
            b_b_co = exp(int + N + P + co + NP + co_N + co_P + co_NP + N_t+ P_t + NP_t + co_N_t + co_P_t + co_NP_t)/(1 + exp(int + N + P + co + NP + co_N + co_P + co_NP + N_t+ P_t + NP_t + co_N_t + co_P_t + co_NP_t)))

avgtpu <- combpu %>%
  gather(key = "treatment", value = "effect") %>%
  mutate(Inoculation = ifelse(grepl("co", treatment, fixed = T), "coinfection", "single"),
         Inoculation = factor(Inoculation, levels = c("single", "coinfection")),
         Nutrient = case_when(substr(treatment, 1, 1) == "l" ~ "low",
                              substr(treatment, 1, 1) == "n" ~ "N",
                              substr(treatment, 1, 1) == "p" ~ "P",
                              substr(treatment, 1, 1) == "b" ~ "N+P"),
         Nutrient = factor(Nutrient, levels = c("low", "N", "P", "N+P")),
         Nutrient_t = case_when(substr(treatment, 3, 3) == "l" ~ "low",
                                substr(treatment, 3, 3) == "n" ~ "N",
                                substr(treatment, 3, 3) == "p" ~ "P",
                                substr(treatment, 3, 3) == "b" ~ "N+P"),
         Nutrient_t = factor(Nutrient_t, levels = c("low", "N", "P", "N+P"))) %>%
  as_tibble()

combpi <- posttpi %>%
  transmute(l_l = exp(int)/(1 + exp(int)),
            l_l_co = exp(int + co)/(1 + exp(int + co)),
            l_n = exp(int + N_t)/(1 + exp(int + N_t)),
            l_n_co = exp(int + N_t + co + co_N_t)/(1 + exp(int + N_t + co + co_N_t)),
            l_p = exp(int + P_t)/(1 + exp(int + P_t)),
            l_p_co = exp(int + P_t + co + co_P_t)/(1 + exp(int + P_t + co + co_P_t)),
            l_b = exp(int + N_t+ P_t + NP_t)/(1 + exp(int + N_t+ P_t + NP_t)),
            l_b_co = exp(int + N_t+ P_t + co + NP_t + co_N_t + co_P_t + co_NP_t)/(1 + exp(int + N_t+ P_t + co + NP_t + co_N_t + co_P_t + co_NP_t)),
            n_l = exp(int + N)/(1 + exp(int + N)),
            n_l_co = exp(int + N + co + co_N)/(1 + exp(int + N + co + co_N)),
            n_n = exp(int + N + N_t)/(1 + exp(int + N + N_t)),
            n_n_co = exp(int + N + N_t + co + co_N + co_N_t)/(1 + exp(int + N + N_t + co + co_N + co_N_t)),
            n_p = exp(int + N + P_t)/(1 + exp(int + N + P_t)),
            n_p_co = exp(int + N + P_t + co + co_N + co_P_t)/(1 + exp(int + N + P_t + co + co_N + co_P_t)),
            n_b = exp(int + N + N_t+ P_t + NP_t)/(1 + exp(int + N + N_t+ P_t + NP_t)),
            n_b_co = exp(int + N + N_t+ P_t + co + co_N + NP_t + co_N_t + co_P_t + co_NP_t)/(1 + exp(int + N + N_t+ P_t + co + co_N + NP_t + co_N_t + co_P_t + co_NP_t)),
            p_l = exp(int + P)/(1 + exp(int + P)),
            p_l_co = exp(int + P + co + co_P)/(1 + exp(int + P + co + co_P)),
            p_n = exp(int + P + N_t)/(1 + exp(int + P + N_t)),
            p_n_co = exp(int + P + N_t + co + co_P + co_N_t)/(1 + exp(int + P + N_t + co + co_P + co_N_t)),
            p_p = exp(int + P + P_t)/(1 + exp(int + P + P_t)),
            p_p_co = exp(int + P + P_t + co + co_P + co_P_t)/(1 + exp(int + P + P_t + co + co_P + co_P_t)),
            p_b = exp(int + P + N_t+ P_t + NP_t)/(1 + exp(int + P + N_t+ P_t + NP_t)),
            p_b_co = exp(int + P + N_t+ P_t + co + co_P + NP_t + co_N_t + co_P_t + co_NP_t)/(1 + exp(int + P + N_t+ P_t + co + co_P + NP_t + co_N_t + co_P_t + co_NP_t)),
            b_l = exp(int + N + P + NP)/(1 + exp(int + N + P + NP)),
            b_l_co = exp(int + N + P + co + NP + co_N + co_P + co_NP)/(1 + exp(int + N + P + co + NP + co_N + co_P + co_NP)),
            b_n = exp(int + N + P + NP + N_t)/(1 + exp(int + N + P + NP + N_t)),
            b_n_co = exp(int + N + P + co + NP + co_N + co_P + co_NP + N_t + co_N_t)/(1 + exp(int + N + P + co + NP + co_N + co_P + co_NP + N_t + co_N_t)),
            b_p = exp(int + N + P + NP + P_t)/(1 + exp(int + N + P + NP + P_t)),
            b_p_co = exp(int + N + P + co + NP + co_N + co_P + co_NP + P_t + co_P_t)/(1 + exp(int + N + P + co + NP + co_N + co_P + co_NP + P_t + co_P_t)),
            b_b = exp(int + N + P + NP + N_t+ P_t + NP_t)/(1 + exp(int + N + P + NP + N_t+ P_t + NP_t)),
            b_b_co = exp(int + N + P + co + NP + co_N + co_P + co_NP + N_t+ P_t + NP_t + co_N_t + co_P_t + co_NP_t)/(1 + exp(int + N + P + co + NP + co_N + co_P + co_NP + N_t+ P_t + NP_t + co_N_t + co_P_t + co_NP_t)))

avgtpi <- combpi %>%
  gather(key = "treatment", value = "effect") %>%
  mutate(Inoculation = ifelse(grepl("co", treatment, fixed = T), "coinfection", "single"),
         Inoculation = factor(Inoculation, levels = c("single", "coinfection")),
         Nutrient = case_when(substr(treatment, 1, 1) == "l" ~ "low",
                              substr(treatment, 1, 1) == "n" ~ "N",
                              substr(treatment, 1, 1) == "p" ~ "P",
                              substr(treatment, 1, 1) == "b" ~ "N+P"),
         Nutrient = factor(Nutrient, levels = c("low", "N", "P", "N+P")),
         Nutrient_t = case_when(substr(treatment, 3, 3) == "l" ~ "low",
                                substr(treatment, 3, 3) == "n" ~ "N",
                                substr(treatment, 3, 3) == "p" ~ "P",
                                substr(treatment, 3, 3) == "b" ~ "N+P"),
         Nutrient_t = factor(Nutrient_t, levels = c("low", "N", "P", "N+P"))) %>%
  as_tibble()

combru <- posttru %>%
  transmute(l_l = exp(int)/(1 + exp(int)),
            l_l_co = exp(int + co)/(1 + exp(int + co)),
            l_n = exp(int + N_t)/(1 + exp(int + N_t)),
            l_n_co = exp(int + N_t + co + co_N_t)/(1 + exp(int + N_t + co + co_N_t)),
            l_p = exp(int + P_t)/(1 + exp(int + P_t)),
            l_p_co = exp(int + P_t + co + co_P_t)/(1 + exp(int + P_t + co + co_P_t)),
            l_b = exp(int + N_t+ P_t + NP_t)/(1 + exp(int + N_t+ P_t + NP_t)),
            l_b_co = exp(int + N_t+ P_t + co + NP_t + co_N_t + co_P_t + co_NP_t)/(1 + exp(int + N_t+ P_t + co + NP_t + co_N_t + co_P_t + co_NP_t)),
            n_l = exp(int + N)/(1 + exp(int + N)),
            n_l_co = exp(int + N + co + co_N)/(1 + exp(int + N + co + co_N)),
            n_n = exp(int + N + N_t)/(1 + exp(int + N + N_t)),
            n_n_co = exp(int + N + N_t + co + co_N + co_N_t)/(1 + exp(int + N + N_t + co + co_N + co_N_t)),
            n_p = exp(int + N + P_t)/(1 + exp(int + N + P_t)),
            n_p_co = exp(int + N + P_t + co + co_N + co_P_t)/(1 + exp(int + N + P_t + co + co_N + co_P_t)),
            n_b = exp(int + N + N_t+ P_t + NP_t)/(1 + exp(int + N + N_t+ P_t + NP_t)),
            n_b_co = exp(int + N + N_t+ P_t + co + co_N + NP_t + co_N_t + co_P_t + co_NP_t)/(1 + exp(int + N + N_t+ P_t + co + co_N + NP_t + co_N_t + co_P_t + co_NP_t)),
            p_l = exp(int + P)/(1 + exp(int + P)),
            p_l_co = exp(int + P + co + co_P)/(1 + exp(int + P + co + co_P)),
            p_n = exp(int + P + N_t)/(1 + exp(int + P + N_t)),
            p_n_co = exp(int + P + N_t + co + co_P + co_N_t)/(1 + exp(int + P + N_t + co + co_P + co_N_t)),
            p_p = exp(int + P + P_t)/(1 + exp(int + P + P_t)),
            p_p_co = exp(int + P + P_t + co + co_P + co_P_t)/(1 + exp(int + P + P_t + co + co_P + co_P_t)),
            p_b = exp(int + P + N_t+ P_t + NP_t)/(1 + exp(int + P + N_t+ P_t + NP_t)),
            p_b_co = exp(int + P + N_t+ P_t + co + co_P + NP_t + co_N_t + co_P_t + co_NP_t)/(1 + exp(int + P + N_t+ P_t + co + co_P + NP_t + co_N_t + co_P_t + co_NP_t)),
            b_l = exp(int + N + P + NP)/(1 + exp(int + N + P + NP)),
            b_l_co = exp(int + N + P + co + NP + co_N + co_P + co_NP)/(1 + exp(int + N + P + co + NP + co_N + co_P + co_NP)),
            b_n = exp(int + N + P + NP + N_t)/(1 + exp(int + N + P + NP + N_t)),
            b_n_co = exp(int + N + P + co + NP + co_N + co_P + co_NP + N_t + co_N_t)/(1 + exp(int + N + P + co + NP + co_N + co_P + co_NP + N_t + co_N_t)),
            b_p = exp(int + N + P + NP + P_t)/(1 + exp(int + N + P + NP + P_t)),
            b_p_co = exp(int + N + P + co + NP + co_N + co_P + co_NP + P_t + co_P_t)/(1 + exp(int + N + P + co + NP + co_N + co_P + co_NP + P_t + co_P_t)),
            b_b = exp(int + N + P + NP + N_t+ P_t + NP_t)/(1 + exp(int + N + P + NP + N_t+ P_t + NP_t)),
            b_b_co = exp(int + N + P + co + NP + co_N + co_P + co_NP + N_t+ P_t + NP_t + co_N_t + co_P_t + co_NP_t)/(1 + exp(int + N + P + co + NP + co_N + co_P + co_NP + N_t+ P_t + NP_t + co_N_t + co_P_t + co_NP_t)))

avgtru <- combru %>%
  gather(key = "treatment", value = "effect") %>%
  mutate(Inoculation = ifelse(grepl("co", treatment, fixed = T), "coinfection", "single"),
         Inoculation = factor(Inoculation, levels = c("single", "coinfection")),
         Nutrient = case_when(substr(treatment, 1, 1) == "l" ~ "low",
                              substr(treatment, 1, 1) == "n" ~ "N",
                              substr(treatment, 1, 1) == "p" ~ "P",
                              substr(treatment, 1, 1) == "b" ~ "N+P"),
         Nutrient = factor(Nutrient, levels = c("low", "N", "P", "N+P")),
         Nutrient_t = case_when(substr(treatment, 3, 3) == "l" ~ "low",
                                substr(treatment, 3, 3) == "n" ~ "N",
                                substr(treatment, 3, 3) == "p" ~ "P",
                                substr(treatment, 3, 3) == "b" ~ "N+P"),
         Nutrient_t = factor(Nutrient_t, levels = c("low", "N", "P", "N+P"))) %>%
  as_tibble()

combri <- posttri %>%
  transmute(l_l = exp(int)/(1 + exp(int)),
            l_l_co = exp(int + co)/(1 + exp(int + co)),
            l_n = exp(int + N_t)/(1 + exp(int + N_t)),
            l_n_co = exp(int + N_t + co + co_N_t)/(1 + exp(int + N_t + co + co_N_t)),
            l_p = exp(int + P_t)/(1 + exp(int + P_t)),
            l_p_co = exp(int + P_t + co + co_P_t)/(1 + exp(int + P_t + co + co_P_t)),
            l_b = exp(int + N_t+ P_t + NP_t)/(1 + exp(int + N_t+ P_t + NP_t)),
            l_b_co = exp(int + N_t+ P_t + co + NP_t + co_N_t + co_P_t + co_NP_t)/(1 + exp(int + N_t+ P_t + co + NP_t + co_N_t + co_P_t + co_NP_t)),
            n_l = exp(int + N)/(1 + exp(int + N)),
            n_l_co = exp(int + N + co + co_N)/(1 + exp(int + N + co + co_N)),
            n_n = exp(int + N + N_t)/(1 + exp(int + N + N_t)),
            n_n_co = exp(int + N + N_t + co + co_N + co_N_t)/(1 + exp(int + N + N_t + co + co_N + co_N_t)),
            n_p = exp(int + N + P_t)/(1 + exp(int + N + P_t)),
            n_p_co = exp(int + N + P_t + co + co_N + co_P_t)/(1 + exp(int + N + P_t + co + co_N + co_P_t)),
            n_b = exp(int + N + N_t+ P_t + NP_t)/(1 + exp(int + N + N_t+ P_t + NP_t)),
            n_b_co = exp(int + N + N_t+ P_t + co + co_N + NP_t + co_N_t + co_P_t + co_NP_t)/(1 + exp(int + N + N_t+ P_t + co + co_N + NP_t + co_N_t + co_P_t + co_NP_t)),
            p_l = exp(int + P)/(1 + exp(int + P)),
            p_l_co = exp(int + P + co + co_P)/(1 + exp(int + P + co + co_P)),
            p_n = exp(int + P + N_t)/(1 + exp(int + P + N_t)),
            p_n_co = exp(int + P + N_t + co + co_P + co_N_t)/(1 + exp(int + P + N_t + co + co_P + co_N_t)),
            p_p = exp(int + P + P_t)/(1 + exp(int + P + P_t)),
            p_p_co = exp(int + P + P_t + co + co_P + co_P_t)/(1 + exp(int + P + P_t + co + co_P + co_P_t)),
            p_b = exp(int + P + N_t+ P_t + NP_t)/(1 + exp(int + P + N_t+ P_t + NP_t)),
            p_b_co = exp(int + P + N_t+ P_t + co + co_P + NP_t + co_N_t + co_P_t + co_NP_t)/(1 + exp(int + P + N_t+ P_t + co + co_P + NP_t + co_N_t + co_P_t + co_NP_t)),
            b_l = exp(int + N + P + NP)/(1 + exp(int + N + P + NP)),
            b_l_co = exp(int + N + P + co + NP + co_N + co_P + co_NP)/(1 + exp(int + N + P + co + NP + co_N + co_P + co_NP)),
            b_n = exp(int + N + P + NP + N_t)/(1 + exp(int + N + P + NP + N_t)),
            b_n_co = exp(int + N + P + co + NP + co_N + co_P + co_NP + N_t + co_N_t)/(1 + exp(int + N + P + co + NP + co_N + co_P + co_NP + N_t + co_N_t)),
            b_p = exp(int + N + P + NP + P_t)/(1 + exp(int + N + P + NP + P_t)),
            b_p_co = exp(int + N + P + co + NP + co_N + co_P + co_NP + P_t + co_P_t)/(1 + exp(int + N + P + co + NP + co_N + co_P + co_NP + P_t + co_P_t)),
            b_b = exp(int + N + P + NP + N_t+ P_t + NP_t)/(1 + exp(int + N + P + NP + N_t+ P_t + NP_t)),
            b_b_co = exp(int + N + P + co + NP + co_N + co_P + co_NP + N_t+ P_t + NP_t + co_N_t + co_P_t + co_NP_t)/(1 + exp(int + N + P + co + NP + co_N + co_P + co_NP + N_t+ P_t + NP_t + co_N_t + co_P_t + co_NP_t)))

avgtri <- combri %>%
  gather(key = "treatment", value = "effect") %>%
  mutate(Inoculation = ifelse(grepl("co", treatment, fixed = T), "coinfection", "single"),
         Inoculation = factor(Inoculation, levels = c("single", "coinfection")),
         Nutrient = case_when(substr(treatment, 1, 1) == "l" ~ "low",
                              substr(treatment, 1, 1) == "n" ~ "N",
                              substr(treatment, 1, 1) == "p" ~ "P",
                              substr(treatment, 1, 1) == "b" ~ "N+P"),
         Nutrient = factor(Nutrient, levels = c("low", "N", "P", "N+P")),
         Nutrient_t = case_when(substr(treatment, 3, 3) == "l" ~ "low",
                                substr(treatment, 3, 3) == "n" ~ "N",
                                substr(treatment, 3, 3) == "p" ~ "P",
                                substr(treatment, 3, 3) == "b" ~ "N+P"),
         Nutrient_t = factor(Nutrient_t, levels = c("low", "N", "P", "N+P"))) %>%
  as_tibble()


#### figure of category averages ####

pcpu <- avgcpu %>%
  group_by(treatment, Nutrient, Inoculation) %>%
  mean_hdi() %>%
  ggplot(aes(x = Nutrient, y = effect,  color = Nutrient)) +
  geom_pointinterval(aes(shape = Inoculation), fatten_point = 2.5, size_range = c(0.4, 0.6), position = position_dodge(0.3), fill = "white", show.legend = F) +
  theme_bw() +
  theme(axis.title.y = element_text(color = "black", size = lg_txt),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_blank(),
        plot.title = element_text(color = "black", size = lg_txt, hjust = 0.5),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.key.width = unit(1, "cm")) +
  scale_colour_manual(values = col_pal) +
  scale_shape_manual(values = c(19, 21)) +
  xlab("Nutrient") +
  ylab("Est. ln(PAV density)") +
  ggtitle("Uninformative")

pcpi <- avgcpi %>%
  group_by(treatment, Nutrient, Inoculation) %>%
  mean_hdi() %>%
  ggplot(aes(x = Nutrient, y = effect,  color = Nutrient)) +
  geom_point(aes(shape = Inoculation), position = position_dodge(0.3)) +
  geom_pointinterval(aes(shape = Inoculation), fatten_point = 2.5, size_range = c(0.4, 0.6), position = position_dodge(0.3), fill = "white", show.legend = F) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_blank(),
        plot.title = element_text(color = "black", size = lg_txt, hjust = 0.5),
        legend.title = element_text(color = "black", size = lg_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.spacing.y = unit(-0.1, "mm"), 
        legend.key = element_rect(size = 0.5, color = "white"),
        legend.key.size = unit(0.7, 'lines')) +
  scale_colour_manual(values = col_pal) +
  scale_shape_manual(values = c(19, 21), guide = F) +
  ggtitle("Informative")

pcru <- avgcru %>%
  group_by(treatment, Nutrient, Inoculation) %>%
  mean_hdi() %>%
  ggplot(aes(x = Nutrient, y = effect,  color = Nutrient)) +
  geom_pointinterval(aes(shape = Inoculation), fatten_point = 2.5, size_range = c(0.4, 0.6), position = position_dodge(0.3), fill = "white", show.legend = F) +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_blank(),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.key.width = unit(1, "cm")) +
  scale_colour_manual(values = col_pal) +
  scale_shape_manual(values = c(19, 21)) +
  ylab("Est. ln(RPV density)") +
  xlab("Nutrient")

pcri <- avgcri %>%
  mutate(Inoculation = recode(Inoculation, coinfection = "co")) %>%
  group_by(treatment, Nutrient, Inoculation) %>%
  mean_hdi() %>%
  ggplot(aes(x = Nutrient, y = effect,  color = Nutrient)) +
  geom_point(aes(shape = Inoculation), position = position_dodge(0.3)) +
  geom_pointinterval(aes(shape = Inoculation), fatten_point = 2.5, size_range = c(0.4, 0.6), position = position_dodge(0.3), fill = "white", show.legend = F) +
  theme_bw() +
  theme(axis.title.x = element_text(color = "black", size = lg_txt),
        axis.title.y = element_blank(),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_blank(),
        legend.title = element_text(color = "black", size = lg_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.spacing.y = unit(-0.1, "mm"), 
        legend.key = element_rect(size = 0.5, color = "white"),
        legend.key.size = unit(0.7, 'lines')) +
  scale_colour_manual(values = col_pal, guide = F) +
  scale_shape_manual(values = c(19, 21), name = "Inoculation") +
  xlab("Nutrient")

ptpu <- avgtpu %>%
  group_by(treatment, Nutrient, Nutrient_t, Inoculation) %>%
  mean_hdi() %>%
  ggplot(aes(x = Nutrient, y = effect)) +
  geom_pointinterval(aes(shape = Inoculation,  color = Nutrient), fatten_point = 2.5, size_range = c(0.3, 0.4), position = position_dodge(0.5), fill = "white") +
  facet_wrap(~Nutrient_t, nrow = 1, strip.position = "bottom") +
  theme_bw() +
  theme(axis.title.y = element_text(color = "black", size = lg_txt),
        axis.title.x = element_blank(),
        axis.text.y = element_text(color = "black", size = sm_txt),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(color = "black", size = sm_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.85, 0.15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(0, "lines")) +
  scale_colour_manual(values = col_pal, guide = F) +
  scale_shape_manual(values = c(19, 21), guide = F) +
  ylab("Est. PAV transmission") 

ptpi <- avgtpi %>%
  group_by(treatment, Nutrient, Nutrient_t, Inoculation) %>%
  mean_hdi() %>%
  ggplot(aes(x = Nutrient, y = effect)) +
  geom_pointinterval(aes(shape = Inoculation,  color = Nutrient), fatten_point = 2.5, size_range = c(0.3, 0.4), position = position_dodge(0.5), fill = "white") +
  facet_wrap(~Nutrient_t, nrow = 1, strip.position = "bottom") +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(color = "black", size = sm_txt),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(color = "black", size = sm_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.85, 0.15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(0, "lines")) +
  scale_colour_manual(values = col_pal, guide = F) +
  scale_shape_manual(values = c(19, 21), guide = F) 

ptru <- avgtru %>%
  group_by(treatment, Nutrient, Nutrient_t, Inoculation) %>%
  mean_hdi() %>%
  ggplot(aes(x = Nutrient, y = effect)) +
  geom_pointinterval(aes(shape = Inoculation,  color = Nutrient), fatten_point = 2.5, size_range = c(0.3, 0.4), position = position_dodge(0.5), fill = "white") +
  facet_wrap(~Nutrient_t, nrow = 1, strip.position = "bottom") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "black", size = lg_txt),
        axis.title.y = element_text(color = "black", size = lg_txt, hjust = 0.3),
        axis.text.y = element_text(color = "black", size = sm_txt),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(color = "black", size = sm_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.85, 0.15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(0, "lines")) +
  scale_colour_manual(values = col_pal, guide = F) +
  scale_shape_manual(values = c(19, 21), guide = F) +
  ylab("Est. RPV transmission") +
  xlab("Receiving plant nutrient") 

ptri <- avgtri %>%
  group_by(treatment, Nutrient, Nutrient_t, Inoculation) %>%
  mean_hdi() %>%
  ggplot(aes(x = Nutrient, y = effect)) +
  geom_pointinterval(aes(shape = Inoculation,  color = Nutrient), fatten_point = 2.5, size_range = c(0.3, 0.4), position = position_dodge(0.5), fill = "white") +
  facet_wrap(~Nutrient_t, nrow = 1, strip.position = "bottom") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(color = "black", size = lg_txt),
        axis.text.y = element_text(color = "black", size = sm_txt),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(color = "black", size = sm_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.85, 0.15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(0, "lines")) +
  scale_colour_manual(values = col_pal, guide = F) +
  scale_shape_manual(values = c(19, 21), guide = F) +
  xlab("Receiving plant nutrient") 


#### combine plots ####

# extract legends
nleg <- get_legend(pcpi)
ileg <- get_legend(pcri)

# combine
plot1 <- cowplot::plot_grid(pcpu, pcpi + theme(legend.position = "none"), nleg, 
                  labels = c("a", "b"), 
                  nrow = 1,
                  label_size = lg_txt, 
                  rel_widths = c(1, 1, 0.3),
                  label_x = c(0, -0.03))

plot2 <- cowplot::plot_grid(pcru, pcri + theme(legend.position = "none"), ileg,
                            labels = c("c", "d"), 
                            nrow = 1,
                            label_size = lg_txt, 
                            rel_widths = c(1, 1, 0.3),
                            label_x = c(0, -0.03))

plot3 <- cowplot::plot_grid(ptpu, ptpi, 
                  labels = c("e", "f"), 
                  nrow = 1,
                  label_size = lg_txt, 
                  rel_heights = c(0.76, 0.76),
                  label_x = c(0, -0.01))

plot4 <- cowplot::plot_grid(ptru, ptri, 
                            labels = c("g", "h"), 
                            nrow = 1,
                            label_size = lg_txt, 
                            label_x = c(0, -0.01))

plot <- cowplot::plot_grid(plot1, plot2, plot3, plot4, ncol = 1)

# print
pdf("./output/figure_S1_prior_comparison.pdf", width = 6, height = 7)
plot
dev.off()


