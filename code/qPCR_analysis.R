## Goal: analyze qPCR data


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
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
         role = PAV_role,
         rel_quant.mg = quant.mg/max(quant.mg, na.rm = T),
         rel_log_quant.mg = log_quant.mg/max(log_quant.mg))

RPVdat <- dat5 %>%
  filter(!is.na(RPV_role) & !is.na(RPV_quant.mg)) %>%
  mutate(quant.mg = RPV_quant.mg,
         log_quant.mg = log(quant.mg + 1),
         role = RPV_role,
         rel_quant.mg = quant.mg/max(quant.mg, na.rm = T),
         rel_log_quant.mg = log_quant.mg/max(log_quant.mg))


#### estimate r, K, and N0 parameters ####

# # select invader and uninvaded
# PAVUIdat <- PAVdat %>%
#   filter(PAV_role %in% c("uninvaded", "invader"))
# 
# RPVUIdat <- RPVdat %>%
#   filter(RPV_role %in% c("uninvaded", "invader"))
# 
# # visualize growth of invaders and single as one process
# PAVUIdat %>%
#   ggplot(aes(dpiUI, PAV_quant.mg, color = PAV_role)) +
#   stat_summary(geom = "line", fun = "mean") +
#   geom_point() +
#   facet_wrap(~ nutrient)
# # invaders don't reach K
# 
# RPVUIdat %>%
#   ggplot(aes(dpiUI, RPV_quant.mg, color = RPV_role)) +
#   stat_summary(geom = "line", fun = "mean") +
#   geom_point() +
#   facet_wrap(~ nutrient)
# # invaders don't reach K
# 
# # simulated data/starting values
# PAVUIdat %>%
#   filter(dpiUI == 5) %>%
#   group_by(nutrient) %>%
#   summarise(minVal = min(quant.mg),
#             meanVal = mean(quant.mg))
# 
# PAVUIsim <- PAVUIdat %>%
#   select(dpiUI, PAV_role) %>%
#   unique() %>%
#   mutate(N0 = 1,
#          r = 0.4,
#          K = 1500,
#          quant.mg = K*N0/(N0 + (K-N0) * exp(-r * dpiUI)))
# 
# # visualize
# ggplot(PAVUIdat, aes(dpiUI, quant.mg)) +
#   geom_point() +
#   geom_line(data = PAVUIsim) +
#   facet_wrap(~nutrient)
# 
# # simulated data/starting values
# RPVUIdat %>%
#   filter(dpiUI == 5) %>%
#   group_by(nutrient) %>%
#   summarise(minVal = min(quant.mg),
#             meanVal = mean(quant.mg))
# 
# RPVUIsim <- RPVUIdat %>%
#   select(dpiUI) %>%
#   unique() %>%
#   mutate(N0 = 1e4,
#          r = 0.4,
#          K = 5e6,
#          quant.mg = K*N0/(N0 + (K-N0) * exp(-r * dpiUI)))
# 
# # visualize
# ggplot(RPVUIdat, aes(dpiUI, quant.mg)) +
#   geom_point() +
#   geom_line(data = RPVUIsim) +
#   facet_wrap(~nutrient)
# 
# # distributions
# x <- seq(0, 10, length.out = 100)
# y <- dgamma(x, shape = 2, scale = 1) # note that this scale is 1/(stan scale)
# plot(x, y, type = "l")
# 
# # fit models
# PAVUImod1 <- brm(data = PAVUIdat, family = gaussian,
#                formula = bf(quant.mg ~ K*N0/(N0 + (K-N0) * exp(-r * dpiUI)),
#                             N0 ~ 1,
#                             r ~ 1,
#                             K ~ 1,
#                             nl = T),
#                prior <- c(prior(gamma(2, 1), nlpar = "N0", lb = 0),
#                           prior(normal(0, 1), nlpar = "r"),
#                           prior(normal(1500, 100), nlpar = "K"),
#                           prior(cauchy(0, 1), class = sigma)),
#                iter = 6000, warmup = 1000, chains = 1, cores = 1)
# 
# summary(PAVUImod1)
# 
# PAVUImod2 <- brm(data = PAVUIdat, family = gaussian,
#                  formula = bf(quant.mg ~ K*N0/(N0 + (K-N0) * exp(-r * dpiUI)),
#                               N0 ~ 1,
#                               r ~ highN*highP,
#                               K ~ highN*highP,
#                               nl = T),
#                  prior <- c(prior(gamma(2, 1), nlpar = "N0", lb = 0),
#                             prior(normal(0, 1), nlpar = "r"),
#                             prior(normal(1500, 100), nlpar = "K"),
#                             prior(cauchy(0, 1), class = sigma)),
#                  iter = 6000, warmup = 1000, chains = 1, cores = 1)
# 
# summary(PAVUImod2)
# 
# # visualize
# PAVUIsim1 <- PAVUIdat %>%
#   select(nutrient, highN, highP, dpiUI) %>%
#   unique() %>%
#   mutate(quant.mg = fitted(PAVUImod2, newdata = .)[, "Estimate"])
# 
# ggplot(PAVUIdat, aes(dpiUI, quant.mg)) +
#   geom_point() +
#   geom_line(data = PAVUIsim1) +
#   facet_wrap(~nutrient)

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
  mutate(log_quant.mg = predict(PAVImod5, newdata = .),
         quant.mg = exp(log_quant.mg))

ggplot(PAVIdat, aes(dpiI, quant.mg)) +
  geom_point() +
  geom_line(data = PAVIsim1) +
  facet_wrap(~nutrient)

ggplot(PAVIdat, aes(dpiI, log_quant.mg)) +
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

# r values
coef(PAVImod5)[2]
coef(RPVImod5)[2]

# N0 values
exp(coef(PAVImod5)[1])
exp(coef(RPVImod5)[1])


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

# linear model with no effects of nutrients
PAVURmod1 <- lm(log_quant.mg ~ dpiUR, data = PAVURdat)
summary(PAVURmod1)
RPVURmod1 <- lm(log_quant.mg ~ dpiUR, data = RPVURdat)
summary(RPVURmod1)

# polynomial model with no effects of nutrients
PAVURmod2 <- lm(log_quant.mg ~ dpiUR + I(dpiUR^2), data = PAVURdat)
summary(PAVURmod2)
RPVURmod2 <- lm(log_quant.mg ~ dpiUR + I(dpiUR^2), data = RPVURdat)
summary(RPVURmod2)

# compare shapes
anova(PAVURmod1, PAVURmod2) # not different
anova(RPVURmod1, RPVURmod2) # quadratic

# full model
PAVURmod3 <- lm(log_quant.mg ~ highN * highP * role  * dpiUR, data = PAVURdat)
summary(PAVURmod3)
RPVURmod3 <- lm(log_quant.mg ~ highN * highP * role  * (dpiUR + I(dpiUR^2)), data = RPVURdat)
summary(RPVURmod3)

drop1(PAVURmod3, test = "F")
drop1(RPVURmod3, test = "F")

# remove 4-way interactions
PAVURmod4 <- update(PAVURmod3, .~. -highN:highP:role:dpiUR)
summary(PAVURmod4)
RPVURmod4 <- update(RPVURmod3, .~. -highN:highP:role:dpiUR-highN:highP:role:I(dpiUR^2))
summary(RPVURmod4)

drop1(PAVURmod4, test = "F")
drop1(RPVURmod4, test = "F")

# remove 3-way interactions
PAVURmod5 <- update(PAVURmod4, .~. -highN:highP:role-highN:highP:dpiUR-highN:role:dpiUR-highP:role:dpiUR)
summary(PAVURmod5)
RPVURmod5 <- update(RPVURmod4, .~. -highN:highP:role-highN:highP:dpiUR-highN:role:dpiUR-highP:role:dpiUR-
                      highN:highP:I(dpiUR^2)-highN:role:I(dpiUR^2)-highP:role:I(dpiUR^2))
summary(RPVURmod5)

drop1(PAVURmod5, test = "F")
drop1(RPVURmod5, test = "F")

# remove 2-way interactions
PAVURmod6 <- lm(log_quant.mg ~ highN + highP + role  + dpiUR, data = PAVURdat)
summary(PAVURmod6)
RPVURmod6 <- lm(log_quant.mg ~ highN + highP + role  + dpiUR + I(dpiUR^2), data = RPVURdat)
summary(RPVURmod6)

drop1(PAVURmod6, test = "F") # none sig
drop1(RPVURmod6, test = "F") # keep time

# simplified model
PAVURmod7 <- lm(log_quant.mg ~ 1, data = PAVURdat)
summary(PAVURmod7)
summary(RPVURmod2)

add1(PAVURmod7, ~ highN + highP + role + dpiUR, test = "F") # no effect
add1(RPVURmod2, ~ .+ highN + highP + role, test = "F") # no effect

#plot(PAVURmod7)
#plot(RPVURmod2)

# PAV value
exp(coef(PAVURmod7)[1])


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
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.position = "none",
        legend.key.size = unit(5, "mm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        strip.placement = "outside",
        plot.title = element_text(size = 12, vjust = 0))

# dodge points
dodge_size = 2

# shapes
shape_pal = c(16, 17, 15, 4, 1, 18)

# max values
maxPAV_quant.mg <- max(PAVdat$quant.mg, na.rm = T)
maxRPV_quant.mg <- max(RPVdat$quant.mg, na.rm = T)

### invasion figure ###
legFig <- PAVdat %>%
  filter(role == "invader") %>%
  mutate(Nutrient = fct_recode(nutrient, "low" = "L", "+N+P" = "NP", "+N" = "N", "+P" = "P") %>%
           fct_relevel("low", "+N", "+P")) %>%
  ggplot(aes(dpiI, rel_quant.mg)) +
  stat_summary(geom = "point", fun = "mean", position = position_dodge(dodge_size), aes(color = Nutrient)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, alpha = 0.7, position = position_dodge(dodge_size), aes(color = Nutrient)) +
  scale_color_viridis_d() +
  theme_bw() +
  theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.background = element_blank(),
        legend.position = c(0.27, 0.7),
        legend.key.size = unit(5, "mm"))

PAVISimFig <- PAVdat %>%
  filter(role == "invader") %>%
  select(role, dpiI) %>%
  unique() %>%
  mutate(quant.mg = predict(PAVImod5, newdata = .),
         quant.se = predict(PAVImod5, newdata = ., se.fit = T)$se.fit,
         rel_quant.mg = exp(quant.mg)/maxPAV_quant.mg,
         rel_quant.max = exp(quant.mg + quant.se)/maxPAV_quant.mg,
         rel_quant.min = exp(quant.mg - quant.se)/maxPAV_quant.mg,
         Virus = "PAV") %>%
  full_join(tibble(dpiI = seq(min(RPVURdat$dpiI), max(RPVURdat$dpiI), length.out = 50),
                   dpiUR = seq(min(RPVURdat$dpiUR), max(RPVURdat$dpiUR), length.out = 50)) %>%
              mutate(role = "resident") %>%
              mutate(quant.mg = predict(RPVURmod2, newdata = .),
                     quant.se = predict(RPVURmod2, newdata = ., se.fit = T)$se.fit,
                     rel_quant.mg = exp(quant.mg)/maxRPV_quant.mg,
                     rel_quant.max = exp(quant.mg + quant.se)/maxRPV_quant.mg,
                     rel_quant.min = exp(quant.mg - quant.se)/maxRPV_quant.mg,
                     Virus = "RPV"))

PAVIfig <- PAVdat %>%
  filter(role == "invader") %>%
  mutate(Virus = "PAV") %>%
  full_join(RPVdat %>%
              filter(role == "resident") %>%
              mutate(Virus = "RPV")) %>%
  mutate(Nutrient = fct_recode(nutrient, "low" = "L", "+N+P" = "NP", "+N" = "N", "+P" = "P") %>%
           fct_relevel("low", "+N", "+P")) %>%
  ggplot(aes(dpiI, rel_quant.mg)) +
  geom_ribbon(data = PAVISimFig, aes(ymin = rel_quant.min, ymax = rel_quant.max, shape = Virus), fill = "black", color = NA, alpha= 0.3) +
  geom_line(data = PAVISimFig, aes(linetype = Virus)) +
  stat_summary(geom = "point", fun = "mean", position = position_dodge(dodge_size), aes(shape = Virus, color = Nutrient)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, alpha = 0.7, position = position_dodge(dodge_size), aes(shape = Virus, color = Nutrient)) +
  scale_color_viridis_d(guide = F) +
  scale_shape_manual(values = shape_pal[1:2], labels = c("PAV invader", "RPV resident")) +
  scale_linetype_manual(values = c("solid", "dashed"), labels = c("PAV invader", "RPV resident")) +
  ggtitle("(A) PAV invades RPV") +
  xlab("Days post invader inoculation") +
  ylab("Relative virus titer") +
  fig_theme +
  theme(legend.position = c(0.77, 0.87))

RPVISimFig <- RPVdat %>%
  filter(role == "invader") %>%
  select(role, dpiI) %>%
  unique() %>%
  mutate(quant.mg = predict(RPVImod5, newdata = .),
         quant.se = predict(RPVImod5, newdata = ., se.fit = T)$se.fit,
         rel_quant.mg = exp(quant.mg)/maxRPV_quant.mg,
         rel_quant.max = exp(quant.mg + quant.se)/maxRPV_quant.mg,
         rel_quant.min = exp(quant.mg - quant.se)/maxRPV_quant.mg,
         Virus = "RPV") %>%
  full_join(PAVdat %>%
              filter(role == "resident") %>%
              mutate(dpiUR = dpiR - min(dpiR)) %>%
              select(role, dpiI, dpiUR) %>%
              unique() %>%
              mutate(quant.mg = predict(PAVURmod7, newdata = .),
                     quant.se = predict(PAVURmod7, newdata = ., se.fit = T)$se.fit,
                     rel_quant.mg = exp(quant.mg)/maxPAV_quant.mg,
                     rel_quant.max = exp(quant.mg + quant.se)/maxPAV_quant.mg,
                     rel_quant.min = exp(quant.mg - quant.se)/maxPAV_quant.mg,
                     Virus = "PAV"))

RPVIfig <- RPVdat %>%
  filter(role == "invader") %>%
  mutate(Virus = "RPV") %>%
  full_join(PAVdat %>%
              filter(role == "resident") %>%
              mutate(Virus = "PAV")) %>%
  mutate(Nutrient = fct_recode(nutrient, "low" = "L", "+N+P" = "NP", "+N" = "N", "+P" = "P") %>%
           fct_relevel("low", "+N", "+P")) %>%
  ggplot(aes(dpiI, rel_quant.mg)) +
  geom_ribbon(data = RPVISimFig, aes(ymin = rel_quant.min, ymax = rel_quant.max, shape = Virus), fill = "black", color = NA, alpha= 0.3) +
  geom_line(data = RPVISimFig, aes(linetype = Virus)) +
  stat_summary(geom = "point", fun = "mean", position = position_dodge(dodge_size), aes(shape = Virus, color = Nutrient)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, alpha = 0.7, position = position_dodge(dodge_size), aes(shape = Virus, color = Nutrient)) +
  scale_color_viridis_d(guide = F) +
  scale_shape_manual(values = shape_pal[3:4], labels = c("PAV resident", "RPV invader")) +
  scale_linetype_manual(values = c("dashed", "solid"), labels = c("PAV resident", "RPV invader")) +
  ggtitle("(B) RPV invades PAV") +
  xlab("Days post invader inoculation") +
  fig_theme +
  theme(axis.title.y = element_blank(),
        legend.position = c(0.77, 0.87))


### resident figure ###

RPVRSimFig <- RPVURdat %>%
  select(role) %>%
  unique() %>%
  expand_grid(tibble(dpiR = seq(min(RPVURdat$dpiR), max(RPVURdat$dpiR), length.out = 50),
                     dpiUR = seq(min(RPVURdat$dpiUR), max(RPVURdat$dpiUR), length.out = 50))) %>%
  mutate(quant.mg = predict(RPVURmod2, newdata = .),
         quant.se = predict(RPVURmod2, newdata = ., se.fit = T)$se.fit,
         rel_quant.mg = exp(quant.mg)/maxRPV_quant.mg,
         rel_quant.max = exp(quant.mg + quant.se)/maxRPV_quant.mg,
         rel_quant.min = exp(quant.mg - quant.se)/maxRPV_quant.mg,
         Nutrient = "low")

RPVRfig <- RPVdat %>%
  filter(role != "invader") %>%
  mutate(Nutrient = fct_recode(nutrient, "low" = "L", "+N+P" = "NP", "+N" = "N", "+P" = "P") %>%
           fct_relevel("low", "+N", "+P")) %>%
  ggplot(aes(dpiR, rel_quant.mg, group = interaction(role, Nutrient), shape = role, linetype = role, color = Nutrient)) +
  geom_ribbon(data = RPVRSimFig, aes(ymin = rel_quant.min, ymax = rel_quant.max), fill = "black", color = NA, alpha= 0.3, show.legend = F) +
  geom_line(data = RPVRSimFig, color = "black") +
  stat_summary(geom = "point", fun = "mean", position = position_dodge(dodge_size)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, alpha = 0.7, position = position_dodge(dodge_size), linetype = "solid") +
  scale_color_viridis_d(guide = F) +
  scale_shape_manual(values = shape_pal[c(2,5)], name = "Invasion", labels = c("RPV resident", "RPV uninvaded")) +
  scale_linetype_manual(values = c("dashed", "dotted"), name = "Invasion", labels = c("RPV resident", "RPV uninvaded")) +
  ggtitle("(C) Established RPV") +
  xlab("Days post resident inoculation") +
  ylab("Relative virus titer") +
  fig_theme +
  theme(legend.position = c(0.77, 0.87))

PAVRSimFig <- PAVURdat %>%
  select(role, dpiR, dpiUR) %>%
  unique() %>%
  mutate(quant.mg = predict(PAVURmod7, newdata = .),
         quant.se = predict(PAVURmod7, newdata = ., se.fit = T)$se.fit,
         rel_quant.mg = exp(quant.mg)/maxPAV_quant.mg,
         rel_quant.max = exp(quant.mg + quant.se)/maxPAV_quant.mg,
         rel_quant.min = exp(quant.mg - quant.se)/maxPAV_quant.mg,
         Nutrient = "low")

# figures
PAVRfig <- PAVdat %>%
  filter(role != "invader") %>%
  mutate(Nutrient = fct_recode(nutrient, "low" = "L", "+N+P" = "NP", "+N" = "N", "+P" = "P") %>%
           fct_relevel("low", "+N", "+P")) %>%
  ggplot(aes(dpiR, rel_quant.mg, group = interaction(role, Nutrient), shape = role, linetype = role, color = Nutrient)) +
  geom_ribbon(data = PAVRSimFig, aes(ymin = rel_quant.min, ymax = rel_quant.max), fill = "black", color = NA, alpha= 0.3, show.legend = F) +
  geom_line(data = PAVRSimFig, color = "black") +
  stat_summary(geom = "point", fun = "mean", position = position_dodge(dodge_size)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, alpha = 0.7, position = position_dodge(dodge_size), linetype = "solid") +
  scale_color_viridis_d(guide = F) +
  scale_shape_manual(values = shape_pal[c(3,6)], name = "Invasion", labels = c("PAV resident", "PAV uninvaded")) +
  scale_linetype_manual(values = c("dashed", "dotted"), name = "Invasion", labels = c("PAV resident", "PAV uninvaded")) +
  ggtitle("(D) Established PAV") +
  xlab("Days post resident inoculation") +
  fig_theme +
  theme(legend.position = c(0.77, 0.87),
        axis.title.y = element_blank())

# add nutrient legend
nut_leg <- get_legend(legFig)
PAVIfig2 <- align_plots(PAVIfig, nut_leg, axis = "l")
PAVIfig3 <- ggdraw(PAVIfig2[[1]]) + draw_plot(PAVIfig2[[2]])

# combine
tiff("output/PAV_RPV_virus_titer_figure.tiff", width = 6.5, height = 6.5, units = "in", res = 300)
plot_grid(PAVIfig3, RPVIfig, RPVRfig, PAVRfig,
          nrow = 2,
          rel_widths = c(1, 0.9, 1, 0.9))
dev.off()

# outlier
filter(PAVdat, role == "uninvaded" & nutrient == "L" & time == 3) %>%
  select(tube_label, quant.mg, expt_notes, extraction_notes)
# looked at raw data from qPCR_data_processing
# one technical replicate was > 1000, the other two were not, their average was not
# averages below 1000 (or minimum standard) were converted to zero
