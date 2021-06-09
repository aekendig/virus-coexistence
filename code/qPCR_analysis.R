## Goal: analyze qPCR data


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(nlme)
library(nlshelper) # anova_nlslist
library(cowplot)

# import data
dat <- read_csv("intermediate-data/qPCR_expt_data_cleaned.csv")


#### edit data ####

# make wide by target
# for each invasion trial, determine whether resident established
# for each single inoculation, determine whether there was contamination
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
         single_cont = case_when(invasion == "S" & first_inoculation == "PAV" & RPV_quant.mg > 0  ~ 1,
                                 invasion == "S" & first_inoculation == "RPV" & PAV_quant.mg > 0  ~ 1,
                                 TRUE ~ 0),
         missing_quant = case_when(is.na(RPV_quant.mg) | is.na(PAV_quant.mg) ~ 1,
                                   TRUE ~ 0),
         PAV_role = case_when(invasion == "S" & first_inoculation == "PAV" ~ "single",
                              invasion == "I" & first_inoculation == "PAV" ~ "resident",
                              invasion == "I" & first_inoculation == "RPV" ~ "invader",
                              TRUE ~ NA_character_),
         RPV_role = case_when(invasion == "S" & first_inoculation == "RPV" ~ "single",
                              invasion == "I" & first_inoculation == "RPV" ~ "resident",
                              invasion == "I" & first_inoculation == "PAV" ~ "invader",
                              TRUE ~ NA_character_),
         dpiSI = case_when(invasion == "S" ~ dpiR, # DPI that can be used for single and invaders
                           invasion == "I" ~ dpiI)) # incorrect for residents

# samples with each issue
dat2 %>%
  filter(missing_quant == 1)
# 56 missing one of the viruses

dat2 %>%
  filter(invasion == "I" & resident_est == 0) 
# 32 failed invasions

dat2 %>%
  filter(invasion == "S" & single_cont == 1)
# 92 contaminated single infections out of 197

dat2 %>%
  filter(invasion == "S") %>%
  ggplot(aes(x = log10(RPV_quant.mg), color = as.factor(single_cont))) +
  geom_density()
# clear separation between contamination and inoculation

dat2 %>%
  filter(invasion == "S") %>%
  ggplot(aes(x = log10(PAV_quant.mg), color = as.factor(single_cont))) +
  geom_density()
# no clear separation

# remove failed invasions
# contaminated single inoculations
# and missing quantities
dat3 <- dat2 %>%
  filter(!(invasion == "I" & resident_est == 0) & 
           !(invasion == "S" & single_cont == 1) &
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
# no single PAV data
# PAV growth may be higher with N

# RPV quantities visualization
dat3 %>%
  filter(!is.na(RPV_role)) %>%
  ggplot(aes(x = time, y = RPV_quant.mg, color = RPV_role)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  facet_wrap(~ nutrient)
# single and resident are very close
# RPV growth may be higher with N

# PAV single - why are they missing?
(PAV_single_summary <- dat2 %>%
  filter(PAV_role == "single") %>%
  summarise(missing = sum(missing_quant),
            contaminated = sum(single_cont),
            total = n(),
            maxRPV = max(RPV_quant.mg, na.rm = T)) %>%
  mutate(log_maxRPV = log10(maxRPV)))
# 92 out of 100 had some RPV
# max RPV is 171937, consistent with density figure

# same for RPV?
dat2 %>%
  filter(RPV_role == "single") %>%
  summarise(missing = sum(missing_quant),
            contaminated = sum(single_cont),
            total = n())
# no contamination

# replace contaminated single inoculations (use dat2)
# remove missing quantities
# make RPV quant NA below threshold
dat4 <- dat2 %>%
  filter(missing_quant != 1) %>%
  mutate(RPV_quant.mg = case_when(RPV_quant.mg < PAV_single_summary$maxRPV ~ NA_real_,
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
         log_quant.mg = log(quant.mg + 1))

RPVdat <- dat5 %>%
  filter(!is.na(RPV_role) & !is.na(RPV_quant.mg)) %>%
  mutate(quant.mg = RPV_quant.mg,
         log_quant.mg = log(quant.mg + 1))


#### compare resident to single ####

# select data
PAVSRdat <- PAVdat %>%
  filter(PAV_role %in% c("single", "resident")) %>%
  mutate(role = fct_recode(PAV_role, uninvaded = "single") %>%
           fct_rev(),
         dpiSR = dpiR - min(dpiR))

RPVSRdat <- RPVdat %>%
  filter(RPV_role %in% c("single", "resident")) %>%
  mutate(role = fct_recode(RPV_role, uninvaded = "single") %>%
           fct_rev(),
         dpiSR = dpiR - min(dpiR))

# visualize
ggplot(PAVSRdat, aes(dpiR, quant.mg, color = role)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ nutrient)

ggplot(RPVSRdat, aes(dpiR, quant.mg, color = role)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ nutrient)

# fit model 
PAVSRmod <- lme(quant.mg ~ highN * highP * role, random = ~  1|time, data = PAVSRdat)
summary(PAVSRmod)
plot(PAVSRmod)
PAVSRmod2 <- lme(log_quant.mg ~ highN * highP * role, random = ~  1|time, data = PAVSRdat)
summary(PAVSRmod2)
plot(PAVSRmod)

RPVSRmod <- lme(quant.mg ~ highN * highP * role, random = ~  1|time, data = RPVSRdat)
summary(RPVSRmod)
plot(RPVSRmod)
RPVSRmod2 <- lme(log_quant.mg ~ highN * highP * role, random = ~  1|time, data = RPVSRdat)
summary(RPVSRmod2)
plot(RPVSRmod)

# polynomial model?
RPVSRmod3 <- lm(log_quant.mg ~ highN * highP * role * (dpiSR + I(dpiSR^2)), data = RPVSRdat)
summary(RPVSRmod3)
plot(RPVSRmod3)
RPVSRmod4 <- lm(log_quant.mg ~ highN * highP * role * dpiSR, data = RPVSRdat)
summary(RPVSRmod4)
AIC(RPVSRmod4, RPVSRmod3)
# 3 is better than 4

# values
PAVSRval <- PAVSRdat %>%
  select(nutrient, highN, highP, role) %>%
  unique() %>%
  mutate(titer = exp(predict(PAVSRmod2, newdata = ., level = 0)))

RPVSRval <- RPVSRdat %>%
  select(nutrient, highN, highP, role) %>%
  unique() %>%
  mutate(titer = exp(predict(RPVSRmod2, newdata = ., level = 0)))


#### estimate growth parameters ####

# visualize growth of invaders and single as one process
PAVdat %>%
  filter(PAV_role %in% c("single", "invader")) %>%
  ggplot(aes(dpiSI, PAV_quant.mg, color = PAV_role)) +
  stat_summary(geom = "line", fun = "mean") +
  geom_point() +
  facet_wrap(~ nutrient)
# invaders don't reach K

RPVdat %>%
  filter(RPV_role %in% c("single", "invader")) %>%
  ggplot(aes(dpiSI, RPV_quant.mg, color = RPV_role)) +
  stat_summary(geom = "line", fun = "mean") +
  geom_point() +
  facet_wrap(~ nutrient)
# invaders don't reach K

# carrying capacity based on previous model
PAVK <- PAVSRval %>%
  filter(role == "uninvaded") %>%
  select(-role) %>%
  rename(K = titer)

RPVK <- RPVSRval %>%
  filter(role == "uninvaded") %>%
  select(-role) %>%
  rename(K = titer)

# select invader
PAVIdat <- PAVdat %>%
  filter(PAV_role == "invader") %>%
  select(set, nutrient, highN, highP, time, dpiI, replicate, quant.mg) %>%
  left_join(PAVK)

RPVIdat <- RPVdat %>%
  filter(RPV_role == "invader") %>%
  select(set, nutrient, highN, highP, time, dpiI, replicate, quant.mg) %>%
  left_join(RPVK)

# simulated data/starting values
PAVIdat %>%
  filter(time == 1) %>%
  group_by(nutrient) %>%
  summarise(minVal = min(quant.mg),
            meanVal = mean(quant.mg))

PAVIsim <- PAVIdat %>%
  select(dpiI) %>%
  unique() %>%
  expand_grid(PAVK) %>%
  mutate(N0 = 1,
         r = 0.3,
         quant.mg = K*N0/(N0 + (K-N0) * exp(-r * dpiI)))

# visualize
ggplot(PAVIdat, aes(dpiI, quant.mg)) +
  geom_point() +
  geom_line(data = PAVIsim) +
  facet_wrap(~nutrient)

# simulated data/starting values
RPVIdat %>%
  filter(time == 1) %>%
  group_by(nutrient) %>%
  summarise(minVal = min(quant.mg),
            meanVal = mean(quant.mg))

RPVIsim <- RPVIdat %>%
  select(dpiI) %>%
  unique() %>%
  expand_grid(RPVK) %>%
  mutate(N0 = 1e4,
         r = 0.3,
         quant.mg = K*N0/(N0 + (K-N0) * exp(-r * dpiI)))

# visualize
ggplot(RPVIdat, aes(dpiI, quant.mg)) +
  geom_point() +
  geom_line(data = RPVIsim) +
  facet_wrap(~nutrient)


### parameters don't vary by treatment ###

# formula
form1 <- formula(quant.mg ~ K*N0/(N0 + (K-N0) * exp(-r * dpiI)))

# model
(PAVmod1 <- nls(form1, data = PAVIdat,
                start = list(N0 = 1, r = 0.3)))

(RPVmod1 <- nls(form1, data = RPVIdat,
                start = list(N0 = 1e4, r = 0.3)))

# visualize
PAVIsim1 <- PAVIdat %>%
  select(dpiI) %>%
  unique() %>%
  expand_grid(PAVK) %>%
  mutate(quant.mg = predict(PAVmod1, newdata = .))

ggplot(PAVIdat, aes(dpiI, quant.mg)) +
  geom_point() +
  geom_line(data = PAVIsim1) +
  facet_wrap(~nutrient)

RPVIsim1 <- RPVIdat %>%
  select(dpiI) %>%
  unique() %>%
  expand_grid(RPVK) %>%
  mutate(quant.mg = predict(RPVmod1, newdata = .))

ggplot(RPVIdat, aes(dpiI, quant.mg)) +
  geom_point() +
  geom_line(data = RPVIsim1) +
  facet_wrap(~nutrient)

# add starting value to dataframe
PAVIdat2 <- PAVIdat %>%
  mutate(N0 = as.numeric(coef(PAVmod1)[1]))

RPVIdat2 <- RPVIdat %>%
  mutate(N0 = as.numeric(coef(RPVmod1)[1]))

# refit model with fixed N0 for comparison
(PAVmod1b <- nls(form1, data = PAVIdat2,
                start = list(r = 0.1)))

(RPVmod1b <- nls(form1, data = RPVIdat2,
                start = list(r = 0.1)))


### parameters vary by treatment (4) ###

# formula
form2 <- formula(quant.mg ~ K*N0/(N0 + (K-N0) * exp(-r * dpiI)) | nutrient)

# model
(PAVmod2 <- nlsList(form2, data = PAVIdat2, 
                    start = list(r = 0.1)))

(RPVmod2 <- nlsList(form2, data = RPVIdat2, 
                    start = list(r = 0.1)))

# compare
anova_nlslist(PAVmod2, PAVmod1b)
anova_nlslist(RPVmod2, RPVmod1b)

# coefficients
PAVcoef2 <- coef(PAVmod2) %>%
  mutate(nutrient = rownames(.))

RPVcoef2 <- coef(RPVmod2) %>%
  mutate(nutrient = rownames(.))

# visualize
PAVIsim2 <- PAVIdat2 %>%
  select(nutrient, dpiI, N0, K) %>%
  unique() %>%
  left_join(PAVcoef2) %>%
  mutate(quant.mg = K*N0/(N0 + (K-N0) * exp(-r * dpiI)))
# the predict function was giving multiple values
# for the same nutrient/dpI combinations

ggplot(PAVIdat2, aes(dpiI, quant.mg)) +
  geom_point() +
  geom_line(data = PAVIsim2) +
  facet_wrap(~nutrient)

RPVIsim2 <- RPVIdat2 %>%
  select(nutrient, dpiI, N0, K) %>%
  unique() %>%
  left_join(RPVcoef2) %>%
  mutate(quant.mg = K*N0/(N0 + (K-N0) * exp(-r * dpiI)))

ggplot(RPVIdat2, aes(dpiI, quant.mg)) +
  geom_point() +
  geom_line(data = RPVIsim2) +
  facet_wrap(~nutrient)


### parameters vary by N level (2) ###

# formula
form3 <- formula(quant.mg ~ K*N0/(N0 + (K-N0) * exp(-r * dpiI)) | highN)

# model
(PAVmod3 <- nlsList(form3, data = PAVIdat2, 
                    start = list(r = 0.1)))

(RPVmod3 <- nlsList(form3, data = RPVIdat2, 
                    start = list(r = 0.1)))

# compare
anova_nlslist(PAVmod3, PAVmod1b)
anova_nlslist(RPVmod3, RPVmod1b)


### parameters vary by P level (2) ###

# formula
form4 <- formula(quant.mg ~ K*N0/(N0 + (K-N0) * exp(-r * dpiI)) | highP)

# model
(PAVmod4 <- nlsList(form4, data = PAVIdat2, 
                    start = list(r = 0.1)))

(RPVmod4 <- nlsList(form4, data = RPVIdat2, 
                    start = list(r = 0.1)))

# compare
anova_nlslist(PAVmod4, PAVmod1b)
anova_nlslist(RPVmod4, RPVmod1b)

# grouping factor does not significantly improve fit for any of the models


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
  filter(PAV_role %in% c("single", "invader")) %>%
  select(PAV_role, dpiI, dpiSI, nutrient, highN, highP) %>%
  unique() %>%
  left_join(PAVcoef2) %>%
  left_join(PAVK) %>%
  mutate(N0 = as.numeric(coef(PAVmod1)[1]),
         PAV_quant.mg = case_when(PAV_role == "single" ~ K,
                                  PAV_role == "invader" ~ K*N0/(N0 + (K-N0) * exp(-r * dpiI))),
         Inoculation = case_when(PAV_role == "single" ~ "uninvaded",
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
  filter(RPV_role %in% c("single", "invader")) %>%
  select(RPV_role, dpiI, dpiSI, nutrient, highN, highP) %>%
  unique() %>%
  left_join(RPVcoef2) %>%
  left_join(RPVK) %>%
  mutate(N0 = as.numeric(coef(RPVmod1)[1]),
         RPV_quant.mg = case_when(RPV_role == "single" ~ K,
                                  RPV_role == "invader" ~ K*N0/(N0 + (K-N0) * exp(-r * dpiI))),
         Inoculation = case_when(RPV_role == "single" ~ "uninvaded",
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
  mutate(titer = predict(PAVSRmod2, newdata = ., level = 0),
         Inoculation = case_when(PAV_role == "single" ~ "uninvaded",
                                 TRUE ~ PAV_role),
         Nutrient = fct_recode(nutrient, "low" = "L", "+N+P" = "NP", "+N" = "N", "+P" = "P") %>%
           fct_relevel("low", "+N", "+P"))

RPVSRsim <- RPVSRdat %>%
  mutate(titer = predict(RPVSRmod2, newdata = ., level = 0),
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
