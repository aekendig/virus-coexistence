## Goal: analyze qPCR data


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(nlme)
library(nlshelper) # anova_nlslist
library(cowplot)
library(brms)

# import data
dat <- read_csv("intermediate-data/qPCR_expt_data_cleaned.csv")


#### edit data ####

# make wide by target
# for each invasion trial, determine whether resident established
# for each uninvaded inoculation, determine whether there was contamination
# identify samples missing one of the viruses
# make column for roles
dat2 <- dat %>%
  select(-c(quant_adj, quant.ul)) %>% # remove other quantity calculations
  pivot_wider(names_from = target,
              values_from = quant.mg,
              names_glue = "{target}_quant.mg") %>% # make wide by target
  mutate(resident_est = case_when(invasion == "I" & first_inoculation == "PAV" & PAV_quant.mg > 0  ~ 1,
                                  invasion == "I" & first_inoculation == "RPV" & RPV_quant.mg > 0  ~ 1,
                                  TRUE ~ 0),
         uninvaded_cont = case_when(invasion == "S" & first_inoculation == "PAV" & RPV_quant.mg > 0  ~ 1,
                                 invasion == "S" & first_inoculation == "RPV" & PAV_quant.mg > 0  ~ 1,
                                 TRUE ~ 0),
         missing_quant = case_when(is.na(RPV_quant.mg) | is.na(PAV_quant.mg) ~ 1,
                                   TRUE ~ 0),
         PAV_role = case_when(invasion == "S" & first_inoculation == "PAV" ~ "uninvaded",
                              invasion == "I" & first_inoculation == "PAV" ~ "resident",
                              invasion == "I" & first_inoculation == "RPV" ~ "invader",
                              TRUE ~ NA_character_),
         RPV_role = case_when(invasion == "S" & first_inoculation == "RPV" ~ "uninvaded",
                              invasion == "I" & first_inoculation == "RPV" ~ "resident",
                              invasion == "I" & first_inoculation == "PAV" ~ "invader",
                              TRUE ~ NA_character_),
         dpiUI = case_when(invasion == "S" ~ dpiR, # DPI that can be used for uninvaded and invaders
                           invasion == "I" ~ dpiI)) # incorrect for residents

# samples with each issue
dat2 %>%
  filter(missing_quant == 1)
# 56 missing one of the viruses

dat2 %>%
  filter(invasion == "I" & resident_est == 0) 
# 32 failed invasions

dat2 %>%
  filter(invasion == "S" & uninvaded_cont == 1)
# 92 contaminated uninvaded infections out of 197

dat2 %>%
  filter(invasion == "S") %>%
  ggplot(aes(x = log10(RPV_quant.mg), color = as.factor(uninvaded_cont))) +
  geom_density()
# clear separation between contamination and inoculation

dat2 %>%
  filter(invasion == "S") %>%
  ggplot(aes(x = log10(PAV_quant.mg), color = as.factor(uninvaded_cont))) +
  geom_density()
# no clear separation

# remove failed invasions
# contaminated single inoculations
# and missing quantities
dat3 <- dat2 %>%
  filter(!(invasion == "I" & resident_est == 0) & 
           !(invasion == "S" & uninvaded_cont == 1) &
           missing_quant != 1)

# sample sizes
dat3 %>%
  group_by(nutrient, first_inoculation, second_inoculation, time, invasion) %>%
  count() %>%
  ggplot(aes(x = n)) +
  geom_histogram()
# may need to remove some treatments for analyses

# PAV quantities visualization
dat3 %>%
  filter(!is.na(PAV_role)) %>%
  ggplot(aes(x = time, y = PAV_quant.mg, color = PAV_role)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  facet_wrap(~ nutrient)
# little uninvaded PAV data
# PAV growth may be higher with N

# RPV quantities visualization
dat3 %>%
  filter(!is.na(RPV_role)) %>%
  ggplot(aes(x = time, y = RPV_quant.mg, color = RPV_role)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  facet_wrap(~ nutrient)
# uninvaded and resident are very close
# RPV growth may be higher with N

# PAV uninvaded - why are they missing?
(PAV_uninvaded_summary <- dat2 %>%
  filter(PAV_role == "uninvaded") %>%
  summarise(missing = sum(missing_quant),
            contaminated = sum(uninvaded_cont),
            total = n(),
            maxRPV = max(RPV_quant.mg, na.rm = T)) %>%
  mutate(log_maxRPV = log10(maxRPV)))
# 92 out of 100 had some RPV
# max RPV is 171937, consistent with density figure

# same for RPV?
dat2 %>%
  filter(RPV_role == "uninvaded") %>%
  summarise(missing = sum(missing_quant),
            contaminated = sum(uninvaded_cont),
            total = n())
# no contamination

# replace contaminated uninvaded inoculations (use dat2)
# remove missing quantities
# make RPV quant NA below threshold
dat4 <- dat2 %>%
  filter(missing_quant != 1) %>%
  mutate(RPV_quant.mg = case_when(RPV_quant.mg < PAV_uninvaded_summary$maxRPV ~ NA_real_,
                                  TRUE ~ RPV_quant.mg),
         resident_est = case_when(invasion == "I" & first_inoculation == "PAV" & PAV_quant.mg > 0  ~ 1,
                                  invasion == "I" & first_inoculation == "RPV" & RPV_quant.mg > 0  ~ 1,
                                  TRUE ~ 0))

# failed invasions removed with missing quantity
dat2 %>%
  filter(invasion == "I" & resident_est == 0 & missing_quant == 1)
# 27

# new failed invasion count
dat4 %>%
  filter(invasion == "I" & resident_est == 0) 
# 6 (only added one additional)

# remove failed invasions
dat5 <- dat4 %>%
  filter(!(invasion == "I" & resident_est == 0))
# recovered almost 100 samples

# sample sizes
dat5 %>%
  group_by(nutrient, first_inoculation, second_inoculation, time, invasion) %>%
  count() %>%
  ggplot(aes(x = n)) +
  geom_histogram()
# still have some small sample sizes

# PAV quantities visualization
dat5 %>%
  filter(!is.na(PAV_role)) %>%
  ggplot(aes(x = time, y = PAV_quant.mg, color = PAV_role)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  facet_wrap(~ nutrient)

# RPV quantities visualization
dat5 %>%
  filter(!is.na(RPV_role)) %>%
  ggplot(aes(x = time, y = RPV_quant.mg, color = RPV_role)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  facet_wrap(~ nutrient)
# invasions still visible

# divide datasets by virus
PAVdat <- dat5 %>%
  filter(!is.na(PAV_role) & !is.na(PAV_quant.mg)) %>%
  mutate(quant.mg = PAV_quant.mg,
         log_quant.mg = log(quant.mg + 1),
         role = PAV_role)

RPVdat <- dat5 %>%
  filter(!is.na(RPV_role) & !is.na(RPV_quant.mg)) %>%
  mutate(quant.mg = RPV_quant.mg,
         log_quant.mg = log(quant.mg + 1),
         role = RPV_role)


#### estimate r, K, and N0 parameters ####

# select invader and uninvaded
PAVUIdat <- PAVdat %>%
  filter(PAV_role %in% c("uninvaded", "invader"))

RPVUIdat <- RPVdat %>%
  filter(RPV_role %in% c("uninvaded", "invader"))

# visualize growth of invaders and single as one process
PAVUIdat %>%
  ggplot(aes(dpiUI, PAV_quant.mg, color = PAV_role)) +
  stat_summary(geom = "line", fun = "mean") +
  geom_point() +
  facet_wrap(~ nutrient)
# invaders don't reach K

RPVUIdat %>%
  ggplot(aes(dpiUI, RPV_quant.mg, color = RPV_role)) +
  stat_summary(geom = "line", fun = "mean") +
  geom_point() +
  facet_wrap(~ nutrient)
# invaders don't reach K

# simulated data/starting values
PAVUIdat %>%
  filter(dpiUI == 5) %>%
  group_by(nutrient) %>%
  summarise(minVal = min(quant.mg),
            meanVal = mean(quant.mg))

PAVUIsim <- PAVUIdat %>%
  select(dpiUI, PAV_role) %>%
  unique() %>%
  mutate(N0 = 1,
         r = 0.4,
         K = 1500,
         quant.mg = K*N0/(N0 + (K-N0) * exp(-r * dpiUI)))

# visualize
ggplot(PAVUIdat, aes(dpiUI, quant.mg)) +
  geom_point() +
  geom_line(data = PAVUIsim) +
  facet_wrap(~nutrient)

# simulated data/starting values
RPVUIdat %>%
  filter(dpiUI == 5) %>%
  group_by(nutrient) %>%
  summarise(minVal = min(quant.mg),
            meanVal = mean(quant.mg))

RPVUIsim <- RPVUIdat %>%
  select(dpiUI) %>%
  unique() %>%
  mutate(N0 = 1e4,
         r = 0.4,
         K = 5e6,
         quant.mg = K*N0/(N0 + (K-N0) * exp(-r * dpiUI)))

# visualize
ggplot(RPVUIdat, aes(dpiUI, quant.mg)) +
  geom_point() +
  geom_line(data = RPVUIsim) +
  facet_wrap(~nutrient)

# distributions
x <- seq(0, 10, length.out = 100)
y <- dgamma(x, shape = 2, scale = 1) # note that this scale is 1/(stan scale)
plot(x, y, type = "l")

# fit models
PAVUImod1 <- brm(data = PAVUIdat, family = gaussian,
               formula = bf(quant.mg ~ K*N0/(N0 + (K-N0) * exp(-r * dpiUI)),
                            N0 ~ 1,
                            r ~ 1,
                            K ~ 1,
                            nl = T),
               prior <- c(prior(gamma(2, 1), nlpar = "N0", lb = 0),
                          prior(normal(0, 1), nlpar = "r"),
                          prior(normal(1500, 100), nlpar = "K"),
                          prior(cauchy(0, 1), class = sigma)),
               iter = 6000, warmup = 1000, chains = 1, cores = 1)

summary(PAVUImod1)

PAVUImod2 <- brm(data = PAVUIdat, family = gaussian,
                 formula = bf(quant.mg ~ K*N0/(N0 + (K-N0) * exp(-r * dpiUI)),
                              N0 ~ 1,
                              r ~ highN*highP,
                              K ~ highN*highP,
                              nl = T),
                 prior <- c(prior(gamma(2, 1), nlpar = "N0", lb = 0),
                            prior(normal(0, 1), nlpar = "r"),
                            prior(normal(1500, 100), nlpar = "K"),
                            prior(cauchy(0, 1), class = sigma)),
                 iter = 6000, warmup = 1000, chains = 1, cores = 1)

summary(PAVUImod2)

# visualize
PAVUIsim1 <- PAVUIdat %>%
  select(nutrient, highN, highP, dpiUI) %>%
  unique() %>%
  mutate(quant.mg = fitted(PAVUImod2, newdata = .)[, "Estimate"])

ggplot(PAVUIdat, aes(dpiUI, quant.mg)) +
  geom_point() +
  geom_line(data = PAVUIsim1) +
  facet_wrap(~nutrient)

# issue is that invader points are pulling down the value of uninvaded so that the population doesn't reach K
# uninvaded is influencing r and invaded is influencing K


#### estimate r and N0 parameters ####

# subset data
PAVIdat <- PAVdat %>%
  filter(role == "invader")

RPVIdat <- RPVdat %>%
  filter(role == "invader")

# simulated data/starting values
PAVIsim <- PAVIdat %>%
  select(dpiI) %>%
  unique() %>%
  mutate(log_N0 = log(1),
         r = 0.3,
         log_quant.mg = log_N0 + r * dpiI)

# visualize
ggplot(PAVIdat, aes(dpiI, log_quant.mg)) +
  geom_point() +
  geom_line(data = PAVIsim) +
  facet_wrap(~nutrient)

# simulated data/starting values
RPVIsim <- RPVIdat %>%
  select(dpiI) %>%
  unique() %>%
  mutate(log_N0 = 10,
         r = 0.3,
         log_quant.mg = log_N0 + r * dpiI)

# visualize
ggplot(RPVIdat, aes(dpiI, log_quant.mg)) +
  geom_point() +
  geom_line(data = RPVIsim) +
  facet_wrap(~nutrient)

# model 1: each nutrient treatment has a different growth rate
PAVImod1 <- lm(log_quant.mg ~ dpiI + highN:dpiI + highP:dpiI + highN:highP:dpiI, data = PAVIdat)
summary(PAVImod1)
RPVImod1 <- lm(log_quant.mg ~ dpiI + highN:dpiI + highP:dpiI + highN:highP:dpiI, data = RPVIdat)
summary(RPVImod1)

# model 2: N and P each have a different growth rate, no interaction
PAVImod2 <- lm(log_quant.mg ~ dpiI + highN:dpiI + highP:dpiI, data = PAVIdat)
summary(PAVImod2)
RPVImod2 <- lm(log_quant.mg ~ dpiI + highN:dpiI + highP:dpiI, data = RPVIdat)
summary(RPVImod2)

# model 3: only N causes a different growth rate
PAVImod3 <- lm(log_quant.mg ~ dpiI + highN:dpiI, data = PAVIdat)
summary(PAVImod3)
RPVImod3 <- lm(log_quant.mg ~ dpiI + highN:dpiI, data = RPVIdat)
summary(RPVImod3)

# model 4: only P causes a different growth rate
PAVImod4 <- lm(log_quant.mg ~ dpiI + highP:dpiI, data = PAVIdat)
summary(PAVImod4)
RPVImod4 <- lm(log_quant.mg ~ dpiI + highP:dpiI, data = RPVIdat)
summary(RPVImod4)

# model 5: all growth rates are the same
PAVImod5 <- lm(log_quant.mg ~ dpiI, data = PAVIdat)
summary(PAVImod5)
RPVImod5 <- lm(log_quant.mg ~ dpiI, data = RPVIdat)
summary(RPVImod5)

# compare
anova(PAVImod1, PAVImod2) # no effect of losing interaction
anova(PAVImod3, PAVImod5) # no effect of losing N
anova(PAVImod4, PAVImod5) # no effect of losing P

anova(RPVImod1, RPVImod2) # no effect of losing interaction
anova(RPVImod3, RPVImod5) # no effect of losing N
anova(RPVImod4, RPVImod5) # no effect of losing P

# visualize
PAVIsim1 <- PAVIdat %>%
  select(nutrient, highN, highP, dpiI) %>%
  unique() %>%
  mutate(quant.mg = exp(predict(PAVImod5, newdata = .)))

ggplot(PAVIdat, aes(dpiI, quant.mg)) +
  geom_point() +
  geom_line(data = PAVIsim1) +
  facet_wrap(~nutrient)

RPVIsim1 <- RPVIdat %>%
  select(nutrient, highN, highP, dpiI) %>%
  unique() %>%
  mutate(quant.mg = exp(predict(RPVImod5, newdata = .)))

ggplot(RPVIdat, aes(dpiI, quant.mg)) +
  geom_point() +
  geom_line(data = RPVIsim1) +
  facet_wrap(~nutrient)


#### compare resident to single ####

# select data
PAVURdat <- PAVdat %>%
  filter(role %in% c("uninvaded", "resident")) %>%
  mutate(dpiUR = dpiR - min(dpiR))

RPVURdat <- RPVdat %>%
  filter(role %in% c("uninvaded", "resident")) %>%
  mutate(dpiUR = dpiR - min(dpiR))

# visualize
ggplot(PAVURdat, aes(dpiUR, log_quant.mg, color = nutrient)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ role)
# role affects intercept and linearity

ggplot(RPVURdat, aes(dpiUR, log_quant.mg, color = nutrient)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ role)
# nutrient affects peak

# no effects of nutrients
PAVURmod1 <- lm(log_quant.mg ~ role  * (dpiUR + I(dpiUR^2)), data = PAVURdat)
summary(PAVURmod1)
RPVURmod1 <- lm(log_quant.mg ~ role  * (dpiUR + I(dpiUR^2)), data = RPVURdat)
summary(RPVURmod1)

#### start here: test above models with nutrients ####


# fit PAV models 
PAVURmod1 <- lm(log_quant.mg ~ lowN * highP * role * (dpiSR + I(dpiSR^2)), data = PAVSRdat)
summary(PAVSRmod1)
# plot(PAVSRmod1)
PAVSRmod2 <- lm(log_quant.mg ~ lowN * highP * role * dpiSR, data = PAVSRdat)
summary(PAVSRmod2)
# plot(PAVSRmod2)
PAVSRmod3 <- lm(log_quant.mg ~ dpiSRF + lowN * highP * role, data = PAVSRdat)
summary(PAVSRmod3)
# plot(PAVSRmod3)

# compare models
AIC(PAVSRmod1, PAVSRmod2, PAVSRmod3)
# separate intercept model is the best

PAVSRmod4 <- lme(log_quant.mg ~ highN * highP * role, random = ~  1|dpiSR, data = PAVSRdat)
summary(PAVSRmod4)
plot(PAVSRmod4)

# check autocorrelation
PAVSRmod5 <- lm(log_quant.mg ~ highN * highP * role, data = PAVSRdat)
PAVSRacf <- PAVSRdat %>%
  mutate(resid = residuals(PAVSRmod5)) %>%
  arrange(set_rep, nutrient, role, dpiSR) %>%
  select(set_rep, nutrient, role, dpiSR, resid) %>%
  group_by(set_rep, nutrient, role) %>%
  complete(dpiSR = 0:21) %>%
  ungroup()

PAVSRacfN = map_df(1:7, 
                   ~PAVSRdat %>%
                     group_by(set_rep, nutrient, role) %>%
                     arrange(set_rep, nutrient, role, dpiSR) %>%
                     summarise(lag = list(diff(dpiSR, lag = .x )))) %>%
  unnest(lag) %>%
  group_by(lag) %>%
  summarise(n = n())

acf(PAVSRacf$resid, lag.max = 7, na.action = na.pass, ci = 0, ylim = c(-0.4, 1))
lines(PAVSRacfN$lag, -qnorm(1-.025)/sqrt(PAVSRacfN$n), lty = 2)
lines(PAVSRacfN$lag, qnorm(1-.025)/sqrt(PAVSRacfN$n), lty = 2)
# no sig autocorrelation

# fit RPV models
RPVSRmod1 <- lm(log_quant.mg ~ lowN * highP * role * (dpiSR + I(dpiSR^2)), data = RPVSRdat)
summary(RPVSRmod1)
# plot(RPVSRmod1)
RPVSRmod2 <- lm(log_quant.mg ~ lowN * highP * role * dpiSR, data = RPVSRdat)
summary(RPVSRmod2)
# plot(RPVSRmod2)
RPVSRmod3 <- lm(log_quant.mg ~ dpiSRF + lowN * highP * role, data = RPVSRdat)
summary(RPVSRmod3)
# plot(RPVSRmod3)

# compare models
AIC(RPVSRmod1, RPVSRmod2, RPVSRmod3)
# separate intercept model is the best
# this might be the cause because the quadratic doesn't need separate intercept, slope, and peaks for every treatment
# model simiplification with quadratic?

RPVSRmod4 <- lme(log_quant.mg ~ highN * highP * role, random = ~  1|dpiSR, data = RPVSRdat)
summary(RPVSRmod4)
plot(RPVSRmod4)

# check autocorrelation
RPVSRmod5 <- lm(log_quant.mg ~ highN * highP * role, data = RPVSRdat)
RPVSRacf <- RPVSRdat %>%
  mutate(resid = residuals(RPVSRmod5)) %>%
  arrange(set_rep, nutrient, role, dpiSR) %>%
  select(set_rep, nutrient, role, dpiSR, resid) %>%
  group_by(set_rep, nutrient, role) %>%
  complete(dpiSR = 0:21) %>%
  ungroup()

RPVSRacfN = map_df(1:7, 
                   ~RPVSRdat %>%
                     group_by(set_rep, nutrient, role) %>%
                     arrange(set_rep, nutrient, role, dpiSR) %>%
                     summarise(lag = list(diff(dpiSR, lag = .x )))) %>%
  unnest(lag) %>%
  group_by(lag) %>%
  summarise(n = n())

acf(RPVSRacf$resid, lag.max = 7, na.action = na.pass, ci = 0, ylim = c(-0.4, 1))
lines(RPVSRacfN$lag, -qnorm(1-.025)/sqrt(RPVSRacfN$n), lty = 2)
lines(RPVSRacfN$lag, qnorm(1-.025)/sqrt(RPVSRacfN$n), lty = 2)
# 7 day lag is negative (hump-shape)

# values
PAVSRdat %>%
  select(nutrient, highN, highP, role) %>%
  unique() %>%
  mutate(titer = exp(predict(PAVSRmod4, newdata = ., level = 0)))

RPVSRdat %>%
  select(nutrient, highN, highP, role) %>%
  unique() %>%
  mutate(titer = exp(predict(RPVSRmod4, newdata = ., level = 0)))


#### visualize ####

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
        legend.position = "none",
        legend.key.size = unit(5, "mm"),
        legend.spacing.y = unit(0.05, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        strip.placement = "outside",
        plot.title = element_text(size = 12, vjust = 0))

# dodge points
dodge_size = 1

### invasion figure ###
PAVSIsim <- PAVdat %>%
  filter(PAV_role == "invader") %>%
  select(PAV_role, dpiI, dpiSI, nutrient, highN, highP) %>%
  unique() %>%
  left_join(PAVcoef2) %>%
  left_join(PAVK) %>%
  mutate(N0 = as.numeric(coef(PAVmod1)[1]),
         PAV_quant.mg = K*N0/(N0 + (K-N0) * exp(-r * dpiI))) %>%
  full_join(PAVSRdat %>%
              filter(role == "uninvaded") %>%
              select(PAV_role, role, dpiI, dpiSI, dpiSR, nutrient, highN, highP) %>%
              unique() %>%
              mutate(PAV_quant.mg = exp(predict(PAVSRmod4, newdata = .)))) %>%
  mutate(Inoculation = case_when(PAV_role == "single" ~ "uninvaded",
                                 PAV_role == "invader" ~ "invader"),
         Nutrient = fct_recode(nutrient, "low" = "L", "+N+P" = "NP", "+N" = "N", "+P" = "P") %>%
           fct_relevel("low", "+N", "+P"))

PAVSIfig <- PAVdat %>%
  filter(PAV_role %in% c("single", "invader")) %>%
  mutate(Inoculation = case_when(PAV_role == "single" ~ "uninvaded",
                                 PAV_role == "invader" ~ "invader"),
         Nutrient = fct_recode(nutrient, "low" = "L", "+N+P" = "NP", "+N" = "N", "+P" = "P") %>%
           fct_relevel("low", "+N", "+P")) %>%
  ggplot(aes(dpiSI, PAV_quant.mg, color = Nutrient, group = interaction(Nutrient, Inoculation))) +
  stat_summary(geom = "point", fun = "mean", position = position_dodge(dodge_size), aes(shape = Inoculation)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, alpha = 0.7, position = position_dodge(dodge_size)) +
  geom_line(data = PAVSIsim, aes(linetype = Inoculation)) +
  scale_color_viridis_d() +
  ggtitle("(A) PAV") +
  xlab("Days post inoculation") +
  ylab(expression(paste("Virus titer (no. ", mg^-1, " tissue)", sep = ""))) +
  fig_theme +
  theme(legend.position = c(0.16, 0.6))

RPVSIsim <- RPVdat %>%
  filter(RPV_role == "invader") %>%
  select(RPV_role, dpiI, dpiSI, nutrient, highN, highP) %>%
  unique() %>%
  left_join(RPVcoef2) %>%
  left_join(RPVK) %>%
  mutate(N0 = as.numeric(coef(RPVmod1)[1]),
         RPV_quant.mg = K*N0/(N0 + (K-N0) * exp(-r * dpiI))) %>%
  full_join(RPVSRdat %>%
              filter(role == "uninvaded") %>%
              select(RPV_role, role, dpiI, dpiSI, dpiSR, nutrient, highN, highP) %>%
              unique() %>%
              mutate(RPV_quant.mg = exp(predict(RPVSRmod4, newdata = .)))) %>%
  mutate(Inoculation = case_when(RPV_role == "single" ~ "uninvaded",
                                 RPV_role == "invader" ~ "invader"),
         Nutrient = fct_recode(nutrient, "low" = "L", "+N+P" = "NP", "+N" = "N", "+P" = "P") %>%
           fct_relevel("low", "+N", "+P"))

RPVSIfig <- RPVdat %>%
  filter(RPV_role %in% c("single", "invader")) %>%
  mutate(Inoculation = case_when(RPV_role == "single" ~ "uninvaded",
                                 RPV_role == "invader" ~ "invader"),
         Nutrient = fct_recode(nutrient, "low" = "L", "+N+P" = "NP", "+N" = "N", "+P" = "P") %>%
           fct_relevel("low", "+N", "+P")) %>%
  ggplot(aes(dpiSI, RPV_quant.mg, color = Nutrient, group = interaction(Nutrient, Inoculation))) +
  stat_summary(geom = "point", fun = "mean", position = position_dodge(dodge_size), aes(shape = Inoculation)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, alpha = 0.7, position = position_dodge(dodge_size)) +
  geom_line(data = RPVSIsim, aes(linetype = Inoculation)) +
  scale_color_viridis_d() +
  ggtitle("(B) RPV") +
  xlab("Days post inoculation") +
  fig_theme +
  theme(axis.title.y = element_blank())

# combine
pdf("output/PAV_RPV_invasion_figure.pdf", width = 6.5, height = 3)
plot_grid(PAVSIfig, RPVSIfig,
          nrow = 1,
          rel_widths = c(1, 0.9))
dev.off()


### resident figure ###

# figure settings
fig_theme2 <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0,"line"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.background = element_blank(),
        legend.position = "none",
        legend.key.height = unit(3, "mm"),
        legend.key.width = unit(7, "mm"),
        legend.spacing.y = unit(0, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 12, vjust = 0))

# data for figure
PAVSRsim <- PAVSRdat %>%
  mutate(titer = predict(PAVSRmod4, newdata = .),
         Inoculation = case_when(PAV_role == "single" ~ "uninvaded",
                                 TRUE ~ PAV_role),
         Nutrient = fct_recode(nutrient, "low" = "L", "+N+P" = "NP", "+N" = "N", "+P" = "P") %>%
           fct_relevel("low", "+N", "+P"))

RPVSRsim <- RPVSRdat %>%
  mutate(titer = predict(RPVSRmod4, newdata = .),
         Inoculation = case_when(RPV_role == "single" ~ "uninvaded",
                                 TRUE ~ RPV_role),
         Nutrient = fct_recode(nutrient, "low" = "L", "+N+P" = "NP", "+N" = "N", "+P" = "P") %>%
           fct_relevel("low", "+N", "+P"))

# figures
PAVSRfig <- PAVSRsim %>%
  ggplot(aes(dpiR, log_quant.mg, color = Nutrient, group = interaction(Nutrient, Inoculation))) +
  stat_summary(geom = "point", fun = "mean", position = position_dodge(dodge_size), aes(shape = Inoculation)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, alpha = 0.7, position = position_dodge(dodge_size)) +
  geom_line(aes(y = titer, linetype = Inoculation)) +
  scale_x_continuous(breaks = c(17, 20, 24, 28, 31)) +
  facet_wrap(~Nutrient) +
  scale_color_viridis_d(guide = F) +
  ggtitle("(A) PAV") +
  xlab("Days post inoculation") +
  ylab(expression(paste("Virus titer (no. ", mg^-1, " tissue) [log]", sep = ""))) +
  fig_theme2 +
  theme(legend.position = c(0.16, 0.1))

RPVSRfig <- RPVSRsim %>%
  ggplot(aes(dpiR, log_quant.mg, color = Nutrient, group = interaction(Nutrient, Inoculation))) +
  stat_summary(geom = "point", fun = "mean", position = position_dodge(dodge_size), aes(shape = Inoculation)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, alpha = 0.7, position = position_dodge(dodge_size)) +
  geom_line(aes(y = titer, linetype = Inoculation)) +
  facet_wrap(~Nutrient) +
  scale_color_viridis_d() +
  ggtitle("(B) RPV") +
  xlab("Days post inoculation") +
  scale_x_continuous(breaks = c(17, 20, 24, 28, 31)) +
  scale_y_continuous(breaks = c(14, 15, 16)) +
  fig_theme2 +
  theme(axis.title.y = element_blank())

# combine
pdf("output/PAV_RPV_resident_figure.pdf", width = 6.5, height = 3)
plot_grid(PAVSRfig, RPVSRfig,
          nrow = 1,
          rel_widths = c(1, 0.9))
dev.off()
