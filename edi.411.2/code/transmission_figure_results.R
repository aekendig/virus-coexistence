## Goal: figure of relationship between concentration and transmission from transmission_analysis.R

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse) # version used: 1.2.1
library(cowplot) # version used: 0.9.4
library(sjPlot)  # version used: 2.7.0
library(brms)  # version used: 2.7.0
library(tidybayes)  # version used: 1.0.4

# import data
datp <- read_csv("./output/transmission_analysis_pav_data.csv")
datr <- read_csv("./output/transmission_analysis_rpv_data.csv")

# load models
load("./output/transmission_pav_up_concentration_interaction_informative.rda")
load("./output/transmission_rpv_up_concentration_interaction_informative.rda")


#### edit data ####

# inoculation and nutrient columns
datp <- datp %>%
  mutate(inoculation = ifelse(co == 0, "single", "co"),
         inoculation = fct_relevel(inoculation, "single"),
         nutrient = fct_relevel(nutrient, "low", "N", "P"),
         nutrient_t = fct_relevel(nutrient_t, "low", "N", "P"),
         log_conc = log(conc))

datr <- datr %>%
  mutate(inoculation = ifelse(co == 0, "single", "co"),
         inoculation = fct_relevel(inoculation, "single"),
         nutrient = fct_relevel(nutrient, "low", "N", "P"),
         nutrient_t = fct_relevel(nutrient_t, "low", "N", "P"),
         log_conc = log(conc))

# posterior samples
postp <- posterior_samples(mpcii)
postr <- posterior_samples(mrcii)

# average concentration values by treatment
avgp_conc_s <- datp %>%
  group_by(inoculation, nutrient, nutrient_t, co, high_N, high_P, high_N_t, high_P_t) %>%
  summarize(conc_s = mean(conc_s))

avgr_conc_s <- datr %>%
  group_by(inoculation, nutrient, nutrient_t, co, high_N, high_P, high_N_t, high_P_t) %>%
  summarize(conc_s = mean(conc_s))

# merge average values and posterior samples
combp <- merge(avgp_conc_s, postp, all = T)
combr <- merge(avgr_conc_s, postr, all = T)

# calculate prevalence
prevfun <- function(dat){
  dat2 <- dat %>%
    mutate(exp_val = exp(b_Intercept + 
                           b_conc_s*conc_s + 
                           b_high_N*high_N + 
                           b_high_P*high_P + 
                           b_high_N_t*high_N_t +
                           b_high_P_t*high_P_t  + 
                           b_co*co + 
                           `b_high_N:high_P`*high_N*high_P + 
                           `b_high_N_t:high_P_t`*high_N_t*high_P_t +
                           `b_conc_s:high_N`*conc_s*high_N +
                           `b_conc_s:high_P`*conc_s*high_P +
                           `b_conc_s:high_N_t`*conc_s*high_N_t +
                           `b_conc_s:high_P_t`*conc_s*high_P_t +
                           `b_high_N:co`*high_N*co +
                           `b_high_P:co`*high_P*co +
                           `b_high_N_t:co`*high_N_t*co +
                           `b_high_P_t:co`*high_P_t*co + 
                           `b_conc_s:high_N:high_P`*conc_s*high_N*high_P +
                           `b_conc_s:high_N_t:high_P_t`*conc_s*high_N_t*high_P_t +
                           `b_high_N:high_P:co`*high_N*high_P*co +
                           `b_high_N_t:high_P_t:co`*high_N_t*high_P_t*co),
           prev = exp_val / (1 + exp_val))
  return(dat2)
} 

avgp <- prevfun(combp) %>% 
  select(inoculation, nutrient, nutrient_t, conc_s, prev) %>% 
  group_by(inoculation, nutrient, nutrient_t, conc_s) %>% 
  mean_hdi() %>%
  ungroup() %>%
  mutate(log_conc = log(conc_s * sd(datp$conc) + mean(datp$conc)))

avgr <- prevfun(combr) %>% 
  select(inoculation, nutrient, nutrient_t, conc_s, prev) %>% 
  group_by(inoculation, nutrient, nutrient_t, conc_s) %>% 
  mean_hdi() %>%
  ungroup() %>%
  mutate(log_conc = log(conc_s * sd(datr$conc) + mean(datr$conc)))

# percentage change
percfun <- function(dat){
  dat2 <- prevfun(dat) %>%
    select(inoculation, nutrient, nutrient_t, prev) %>%
    arrange(inoculation, nutrient, nutrient_t) %>%
    mutate(treatment = paste(inoculation, nutrient, nutrient_t, sep = "_"),
           count = rep(1:nrow(postp), nrow(avgp_conc_s))) %>%
    select(-c(inoculation, nutrient, nutrient_t)) %>%
    spread(key = treatment, value = prev) %>%
    transmute(P_low = (single_P_low - single_low_low)/single_low_low,
              N_low = (single_N_low - single_low_low)/single_low_low,
              b_low = (`single_N+P_low` - single_low_low)/single_low_low,
              low_P = (single_low_P - single_low_low)/single_low_low,
              low_N = (single_low_N - single_low_low)/single_low_low,
              low_b = (`single_low_N+P` - single_low_low)/single_low_low,
              co_P_low = (co_P_low - single_P_low)/single_P_low,
              co_N_low = (co_N_low - single_N_low)/single_N_low,
              co_b_low = (`co_N+P_low` - `single_N+P_low`)/`single_N+P_low`,
              co_low_P = (co_low_P - single_low_P)/single_low_P,
              co_low_N = (co_low_N - single_low_N)/single_low_N,
              co_low_b = (`co_low_N+P` - `single_low_N+P`)/`single_low_N+P`) %>%
    gather(key = "treatment", value = "perc") %>%
    as_tibble() %>%
    group_by(treatment) %>%
    mean_hdi() %>%
    ungroup() 
  
  return(dat2)
} 

percp <- percfun(combp %>% mutate(conc_s = 0))
percp_all <- percfun(combp)
percr <- percfun(combr %>% mutate(conc_s = 0))
percr_all <- percfun(combr)

# fitted values for concentration
datp_pred <- datp %>%
  group_by(inoculation, co, nutrient, nutrient_t, high_N, high_P, high_N_t, high_P_t) %>%
  summarize(min_conc = min(log_conc),
            max_conc = max(log_conc)) %>%
  expand_grid(i = 1:1000) %>%
  group_by(inoculation, co, nutrient, nutrient_t, high_N, high_P, high_N_t, high_P_t) %>%
  mutate(log_conc = seq(unique(min_conc), unique(max_conc), length.out = 1000),
         conc = exp(log_conc),
         conc_s = (conc - mean(datp$conc)) / sd(datp$conc),
         nut_t_title = paste("Recipient: ", nutrient_t, sep = "")) %>%
  ungroup() %>%
  select(-c(min_conc, max_conc, i))

datp_pred <- datp_pred %>%
  cbind(fitted(mpcii, newdata = datp_pred, re_formula = NA, nsamples = 100)) %>%
  rename(pred = Estimate,
         lower = Q2.5,
         upper = Q97.5)

datr_pred <- datr %>%
  group_by(inoculation, co, nutrient, nutrient_t, high_N, high_P, high_N_t, high_P_t) %>%
  summarize(min_conc = min(log_conc),
            max_conc = max(log_conc)) %>%
  expand_grid(i = 1:1000) %>%
  group_by(inoculation, co, nutrient, nutrient_t, high_N, high_P, high_N_t, high_P_t) %>%
  mutate(log_conc = seq(unique(min_conc), unique(max_conc), length.out = 1000),
         conc = exp(log_conc),
         conc_s = (conc - mean(datr$conc)) / sd(datr$conc),
         nut_t_title = paste("Recipient: ", nutrient_t, sep = "")) %>%
  ungroup() %>%
  select(-c(min_conc, max_conc, i))

datr_pred <- datr_pred %>%
  cbind(fitted(mrcii, newdata = datr_pred, re_formula = NA, nsamples = 100)) %>%
  rename(pred = Estimate,
         lower = Q2.5,
         upper = Q97.5)

# check that concentration scaling is consistent
datp %>%
  select(conc, conc_s) %>%
  mutate(type = "data") %>%
  full_join(datp_pred %>%
              select(conc, conc_s) %>%
              mutate(type = "pred")) %>%
  ggplot(aes(x = conc, y = conc_s, color = type)) +
  geom_line()

datr %>%
  select(conc, conc_s) %>%
  mutate(type = "data") %>%
  full_join(datr_pred %>%
              select(conc, conc_s) %>%
              mutate(type = "pred")) %>%
  ggplot(aes(x = conc, y = conc_s, color = type)) +
  geom_line()


#### figure settings ####

# palettes
col_pal = c("black", "darkgoldenrod2", "dodgerblue1", "palegreen4")
line_pal = c("solid", "dashed")
shape_pal = c(19, 21)

# text sizes
sm_txt = 6
lg_txt = 8

# base figure
base_theme <- theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        plot.title = element_text(color = "black", size = lg_txt, hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(color = "black", size = sm_txt),
        strip.background = element_blank(),
        panel.spacing = unit(0, "lines"),
        legend.position = "none")


#### raw data figures ####

plot_pRaw <- datp %>%
  ggplot(aes(x = dpi, y = t_up, color = nutrient)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1, position = position_dodge(0.6), alpha = 0.5, aes(size = inoculation)) +
  stat_summary(fun.y = "mean", geom = "point", size = 1.5, position = position_dodge(0.6), aes(shape = inoculation), fill = "white") +
  stat_summary(fun.y = "mean", geom = "line", position = position_dodge(0.6), aes(linetype = inoculation)) +
  base_theme +
  scale_size_manual(values = c(0.5, 0.5), guide = F) +
  scale_colour_manual(values = col_pal, name = "Source plant\n nutrient") +
  scale_shape_manual(values = shape_pal, name = "Source plant\ninfection") +
  scale_linetype_manual(values = line_pal, name = "Source plant\ninfection") +
  ylim(0, 1.07) + 
  xlab("Days post inoculation") +
  ylab("PAV transmission")

plot_rRaw <- plot_pRaw %+% datr %+%
  ylab("RPV transmission")

plot_pRawS <- plot_pRaw %+% aes(color = nutrient_t)

plot_rRawS <- plot_rRaw %+% aes(color = nutrient_t)


#### concentration-transmission figures ####
  
# overall average
pconc_mean <- log(mean(datp$conc))
rconc_mean <- log(mean(datr$conc))

# order factors
datp_pred$nut_t_title <- factor(datp_pred$nut_t_title, levels = c("Recipient: low", "Recipient: N", "Recipient: P", "Recipient: N+P"))
datr_pred$nut_t_title <- factor(datr_pred$nut_t_title, levels = c("Recipient: low", "Recipient: N", "Recipient: P", "Recipient: N+P"))

# PAV concentration
plot_pCon <- ggplot(datp_pred, aes(x = log_conc, color = nutrient)) +
  geom_vline(xintercept = pconc_mean, color = "gray40", linetype = "dotted") +
  geom_line(aes(y = pred, linetype = inoculation), size = 0.6) +
  facet_wrap(~nut_t_title, nrow = 2) +
  base_theme +
  xlab("ln(PAV density)") +
  ylab("Est. PAV transmission") +
  scale_colour_manual(values = col_pal) +
  scale_linetype_manual(values = line_pal) +
  scale_y_continuous(breaks = c(0, 0.3, 0.6, 0.9))

# RPV concentration
plot_rCon <- ggplot(datr_pred, aes(x = log_conc, color = nutrient)) +
  geom_vline(xintercept = rconc_mean, color = "gray40", linetype = "dotted") +
  geom_line(aes(y = pred, linetype = inoculation), size = 0.6) +
  facet_wrap(~nut_t_title, nrow = 2) +
  base_theme +
  xlab("ln(RPV density)") +
  ylab("Est. RPV transmission") +
  scale_colour_manual(values = col_pal) +
  scale_linetype_manual(values = line_pal) +
  scale_y_continuous(breaks = c(0, 0.3, 0.6, 0.9))


#### transmission estimates figures ####

plot_pTra <- ggplot(avgp, aes(x = nutrient, y = prev)) +
  geom_pointinterval(aes(shape = inoculation,  color = nutrient), fatten_point = 2.5, size_range = c(0.3, 0.4), position = position_dodge(0.5), fill = "white") +
  facet_wrap(~nutrient_t, nrow = 1, strip.position = "bottom") +
  base_theme +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ylab("Est. PAV transmission") +
  xlab("Recipient plant nutrient") +
  scale_colour_manual(values = col_pal) +
  scale_shape_manual(values = shape_pal) +
  scale_y_continuous(breaks = c(0, 0.3, 0.6, 0.9))

plot_rTra <- plot_pTra %+% avgr %+%
  ylab("Est. RPV transmission") 


#### legend ####

# legend theme
leg_theme <- theme(legend.position = "bottom",
                   legend.direction = "horizontal",
                   legend.title = element_text(color = "black", size = lg_txt),
                   legend.text = element_text(color = "black", size = sm_txt),
                   legend.background = element_blank(),
                   legend.key = element_rect(color = "white", size = 0.5),
                   legend.key.size = unit(1.5, "lines"),
                   legend.key.width = unit(1, "cm"),
                   legend.spacing.x = unit(0.001, "mm"),
                   legend.justification = c(0.5, 0.5),
                   plot.margin = unit(c(0, 0, 0, 0), units = "cm"),
                   legend.box.margin = margin(-10, -10, -10, -10))

# legends
leg <- get_legend(plot_pRaw %+% leg_theme)

leg_S <- get_legend(plot_pRaw %+% 
                      scale_colour_manual(values = col_pal, name = "Recipient plant\n nutrient")  %+% 
                      leg_theme)


#### combine plots ####
concplots <- cowplot::plot_grid(plot_pRaw, plot_rRaw,
                                plot_pCon, plot_rCon,
                                plot_pTra, plot_rTra,
                                labels = letters[1:6],
                                label_size = lg_txt,
                                nrow = 3,
                                rel_heights = c(0.6, 1, 0.6))

plot_all <- cowplot::plot_grid(concplots, leg,
                               rel_heights = c(1, 0.05),
                               nrow = 2)


pdf("./output/figure_3_transmission.pdf", height = 7, width = 6)
plot_all
dev.off()

supplot <- cowplot::plot_grid(plot_pRawS, plot_rRawS,
                              labels = letters[1:2],
                              label_size = lg_txt,
                              nrow = 1)

plot_allS <- cowplot::plot_grid(supplot, leg_S,
                               rel_heights = c(1, 0.1),
                               nrow = 2)

pdf("./output/figure_S3_transmission.pdf", width = 6, height = 2.5)
plot_allS
dev.off()


#### print model summaries ####

tab_model(mpcii)
summary(mpcii)
prior_summary(mpcii)
percp
percp_all

tab_model(mrcii)
summary(mrcii)
prior_summary(mrcii)
percr
percr_all


#### transmission values for model ####

# modify prevalence function for density-only effects
prevfun_den <- function(dat){
  dat2 <- dat %>%
    mutate(exp_val = exp(b_Intercept + 
                           b_conc_s*conc_s + 
                           `b_conc_s:high_N`*conc_s*high_N +
                           `b_conc_s:high_P`*conc_s*high_P +
                           `b_conc_s:high_N_t`*conc_s*high_N_t +
                           `b_conc_s:high_P_t`*conc_s*high_P_t +
                           `b_conc_s:high_N:high_P`*conc_s*high_N*high_P +
                           `b_conc_s:high_N_t:high_P_t`*conc_s*high_N_t*high_P_t),
           prev = exp_val / (1 + exp_val))
  return(dat2)
} 

# transmission function
transfun <- function(dat, density_only){
  
  # calculate prevalence
  ifelse(density_only == 0, dat2 <- prevfun(dat), dat2 <- prevfun_den(dat))
  
  # spread by infection type
  # calculate q parameter
  dat3 <- dat2 %>% 
    select(inoculation, nutrient, nutrient_t, prev) %>%
    arrange(inoculation, nutrient, nutrient_t) %>%
    mutate(count = rep(1:nrow(postp), nrow(avgp_conc_s))) %>%
    spread(key = inoculation, value = prev) %>%
    mutate(q = co / single,
           B = single) 
  
  # take mean and 95% CI for B
  dat_B <- dat3 %>%
    select(nutrient, nutrient_t, B) %>%
    group_by(nutrient, nutrient_t) %>% 
    mean_hdi() %>%
    ungroup() %>%
    rename(value = B) %>%
    mutate(parameter = "B")
  
  # take mean and 95% CI for q
  dat_q <- dat3 %>%
    select(nutrient, nutrient_t, q) %>%
    group_by(nutrient, nutrient_t) %>% 
    mean_hdi() %>%
    ungroup() %>%
    rename(value = q) %>%
    mutate(parameter = "q")
  
  # combine outputs
  dat4 <- full_join(dat_B, dat_q)

  return(dat4)
}

# apply to parameters
pav_trans_all <- transfun(combp, 0) %>% 
  mutate(mechanisms = "all")

rpv_trans_all <- transfun(combr, 0) %>% 
  mutate(mechanisms = "all")

pav_trans_non <- transfun(combp %>% mutate(conc_s = 0), 0) %>% 
  mutate(mechanisms = "non-density")

rpv_trans_non <- transfun(combr %>% mutate(conc_s = 0), 0) %>% 
  mutate(mechanisms = "non-density")

pav_trans_den <- transfun(combp, 1) %>% 
  mutate(mechanisms = "density")

rpv_trans_den <- transfun(combr, 1) %>% 
  mutate(mechanisms = "density")

# combine
pav_trans <- full_join(pav_trans_den, pav_trans_non) %>%
  full_join(pav_trans_all) %>%
  mutate(mechanisms = fct_relevel(mechanisms, "density", "non-density"))

rpv_trans <- full_join(rpv_trans_den, rpv_trans_non) %>%
  full_join(rpv_trans_all) %>%
  mutate(mechanisms = fct_relevel(mechanisms, "density", "non-density"))

# subset for same nutrient values
pav_trans_same <- pav_trans %>%
  filter(nutrient == nutrient_t)

rpv_trans_same <- rpv_trans %>%
  filter(nutrient == nutrient_t)

# figure (check against transmission figure - note that q is co/single)
ggplot(pav_trans_same, aes(x = nutrient, y = value)) +
  geom_pointinterval(aes(shape = parameter,  color = nutrient), fatten_point = 2.5, size_range = c(0.3, 0.4), position = position_dodge(0.5), fill = "white") +
  scale_shape_manual(values = shape_pal) +
  scale_color_manual(values = col_pal) +
  facet_wrap(~mechanisms, nrow = 3) +
  base_theme +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

ggplot(rpv_trans_same, aes(x = nutrient, y = value)) +
  geom_pointinterval(aes(shape = parameter,  color = nutrient), fatten_point = 2.5, size_range = c(0.3, 0.4), position = position_dodge(0.5), fill = "white") +
  scale_shape_manual(values = shape_pal) +
  scale_color_manual(values = col_pal) +
  facet_wrap(~mechanisms, nrow = 3) +
  base_theme +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

# save data
write_csv(pav_trans_same, "./output/pav_transmission_values_same_nuts.csv")
write_csv(rpv_trans_same, "./output/rpv_transmission_values_same_nuts.csv")
