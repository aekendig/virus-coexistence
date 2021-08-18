## Goal: summary results and figure of infection model estimates from infection_analysis.R and density model estimates from concentration_analysis.R

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse) # version used: 1.2.1
library(cowplot) # version used: 0.9.4
library(brms)  # version used: 2.7.0
library(tidybayes)  # version used: 1.0.4
library(sjPlot)  # version used: 2.7.0

# import data
pdati <- read_csv("./output/infection_analysis_pav_data.csv")
rdati <- read_csv("./output/infection_analysis_rpv_data.csv")
pdatd <- read_csv("./output/concentration_analysis_pav_data.csv")
rdatd <- read_csv("./output/concentration_analysis_rpv_data.csv")

# load models
load("./output/infection_analysis_uninformative_pav.rda")
load("./output/infection_analysis_uninformative_rpv.rda")
load("./output/concentration_analysis_informative_pav.rda")
load("./output/concentration_analysis_informative_rpv.rda")


#### edit data ####

# inoculation column, nutrient column, predicted values (tried including these in the time series)
pdati <- pdati %>%
  mutate(inoculation = ifelse(co == 0, "single", "co"),
         inoculation = fct_relevel(inoculation, "single"),
         nutrient = fct_relevel(nutrient, "low", "N", "P"))

rdati <- rdati %>%
  mutate(inoculation = ifelse(co == 0, "single", "co"),
         inoculation = fct_relevel(inoculation, "single"),
         nutrient = fct_relevel(nutrient, "low", "N", "P"))

pdatd <- pdatd %>%
  mutate(inoculation = ifelse(co == 0, "single", "co"),
         inoculation = fct_relevel(inoculation, "single"),
         nutrient = fct_relevel(nutrient, "low", "N", "P"))

rdatd <- rdatd %>%
  mutate(inoculation = ifelse(co == 0, "single", "co"),
         inoculation = fct_relevel(inoculation, "single"),
         nutrient = fct_relevel(nutrient, "low", "N", "P"))

# posterior samples
postpi <- posterior_samples(m.bu.p)
postri <- posterior_samples(m.bu.r)
postpd <- posterior_samples(m.li.p)
postrd <- posterior_samples(m.li.r)

# treatments
trt <- pdati %>%
  select(inoculation, nutrient, co, high_N, high_P) %>%
  unique()

# merge treatments and posterior samples
combpi <- merge(trt, postpi, all = T)
combri <- merge(trt, postri, all = T)
combpd <- merge(trt, postpd, all = T)
combrd <- merge(trt, postrd, all = T)

# category average for prevalence
avgfun_prev <- function(dat){
  dat2 <- dat %>%
    mutate(exp_val = exp(b_Intercept +
             b_high_N*high_N +
             b_high_P*high_P +
             b_co*co +
             `b_high_N:high_P`*high_N*high_P +
             `b_co:high_N`*co*high_N +
             `b_co:high_P`*co*high_P +
             `b_co:high_N:high_P`*co*high_N*high_P),
           prev = exp_val / (1 + exp_val))
  
  return(dat2)
} 

avgpi <- avgfun_prev(combpi) %>%
  select(nutrient, inoculation, prev)  %>% 
  group_by(nutrient, inoculation) %>%
  mean_hdi() %>%
  ungroup()

avgri <- avgfun_prev(combri) %>%
  select(nutrient, inoculation, prev)  %>% 
  group_by(nutrient, inoculation) %>%
  mean_hdi() %>%
  ungroup() 

# category average for density
avgfun_dens <- function(dat){
  dat2 <- dat %>%
    mutate(avg = b_Intercept +
             b_high_N*high_N +
             b_high_P*high_P +
             b_co*co +
             `b_high_N:high_P`*high_N*high_P +
             `b_co:high_N`*co*high_N +
             `b_co:high_P`*co*high_P +
             `b_co:high_N:high_P`*co*high_N*high_P) %>%
    select(nutrient, inoculation, avg) %>% 
    group_by(nutrient, inoculation) %>%
    mean_hdi() %>%
    ungroup()
  
  return(dat2)
} 

avgpd <- avgfun_dens(combpd) 
avgrd <- avgfun_dens(combrd) 

# percentage change of prevalence
percfun_prev <- function(dat){
  
  dat2 <- avgfun_prev(dat) %>%
    select(inoculation, nutrient, prev) %>%
    arrange(inoculation, nutrient) %>%
    mutate(treatment = paste(inoculation, nutrient, sep = "_"),
           count = rep(1:nrow(postpi), nrow(trt))) %>%
    select(-c(inoculation, nutrient)) %>%
    spread(key = treatment, value = prev) %>%
    transmute(high_N = (single_N - single_low)/single_low,
              high_P = (single_P - single_low)/single_low,
              high_NP = (`single_N+P` - single_low)/single_low,
              low_co = (co_low - single_low)/single_low,
              N_co = (co_N - single_N)/single_N,
              P_co = (co_P - single_P)/single_P,
              NP_co = (`co_N+P` - `single_N+P`)/`single_N+P`)%>%
    gather(key = "treatment", value = "perc") %>%
    mutate(inoculation = ifelse(grepl("co", treatment, fixed = T), "co", "single"),
           inoculation = factor(inoculation, levels = c("single", "co")),
           nutrient = recode(treatment, high_N = "N", high_P = "P", high_NP = "N+P", low_co = "low", N_co = "N", P_co = "P", NP_co = "N+P"),
           nutrient = factor(nutrient, levels = c("low", "N", "P", "N+P"))) %>%
    select(-treatment) %>%
    as_tibble() %>%
    group_by(inoculation, nutrient) %>%
    mean_hdi() %>%
    ungroup() 
  
  return(dat2)
} 

percpi <- percfun_prev(combpi)
percri <- percfun_prev(combri)

# percentage change of density
percfun_dens <- function(dat){
  dat2 <- dat %>%
    mutate(single = case_when(inoculation == "single" ~ 1,
                              TRUE ~ 0),
           perc = exp(b_high_N*high_N*single +
                        b_high_P*high_P*single +
                        `b_high_N:high_P`*high_N*high_P*single +
                        b_co*co +
                        `b_co:high_N`*high_N*co +
                        `b_co:high_P`*high_P*co +
                        `b_co:high_N:high_P`*high_N*high_P*co) - 1) %>%
    select(nutrient, inoculation, perc) %>% 
    group_by(nutrient, inoculation) %>%
    mean_hdi() %>%
    ungroup()
  
  return(dat2)
} 

percpd <- percfun_dens(combpd)
percrd <- percfun_dens(combrd)


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


#### figures of raw data ####

plotA <- ggplot(pdati, aes(x = dpi, y = present, color = nutrient)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1, position = position_dodge(0.6), aes(size = inoculation)) +
  stat_summary(fun.y = "mean", geom = "point", size = 1.5, position = position_dodge(0.6), aes(shape = inoculation), fill = "white") +
  stat_summary(aes(linetype = inoculation), fun.y = "mean", geom = "line", position = position_dodge(0.6)) +
  base_theme +
  scale_size_manual(values = c(0.5, 0.5)) +
  scale_colour_manual(values = col_pal) +
  scale_shape_manual(values = shape_pal) +
  scale_linetype_manual(values = line_pal) +
  xlab("Days post inoculation") +
  ylab("PAV establishment")

plotB <- plotA %+% rdati %+%
  ylab("RPV establishment")

plotE <- plotA %+% pdatd %+%
  aes(y = log_conc) %+%
  ylab("ln(PAV density)")

plotF <- plotE %+% rdatd %+%
  ylab("ln(RPV density)")


#### figures of category averages ####

plotC <- ggplot(avgpi, aes(x = nutrient, y = prev,  color = nutrient)) +
  geom_pointinterval(aes(shape = inoculation), fatten_point = 2.5, size_range = c(0.4, 0.6), position = position_dodge(0.3), fill = "white", show.legend = F) +
  base_theme +
  scale_colour_manual(values = col_pal) +
  scale_shape_manual(values = shape_pal) +
  xlab("Nutrient") +
  ylab("Est. PAV establishment")

plotD <- plotC %+% avgri %+%
  ylab("Est. RPV establishment") +
  theme(axis.title.y = element_text(color = "black", size = lg_txt, hjust = 0.2))

plotG <- plotC %+% avgpd %+%
  aes(y = avg) %+%
  ylab("Est. ln(PAV density)")

plotH <- plotG %+% avgrd %+%
  ylab("Est. ln(RPV density)")


#### legend ####

# save legend
leg <- get_legend(plotA %+% 
                    theme(legend.position = "bottom",
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
                          legend.box.margin = margin(-10, -10, -10, -10)) %+%
                    scale_size_manual(values = c(0.5, 0.5), name = "Inoculation") +
                    scale_colour_manual(values = col_pal, name = "Nutrient") +
                    scale_shape_manual(values = c(19, 21), name = "Inoculation") +
                    scale_linetype_manual(values = c("solid", "dashed"), name = "Inoculation"))


#### combine plots ####

# combine plots
plots <- cowplot::plot_grid(plotA, plotB, plotC, plotD, plotE, plotF, plotG, plotH,
                              labels = letters[1:8], 
                              label_size = lg_txt, 
                              ncol = 2)

# combine all
plot_all <- cowplot::plot_grid(plots, leg,
                               rel_heights = c(1, 0.05),
                               nrow = 2)

# print
pdf("./output/figure_2_infection_density.pdf", width = 6, height = 7)
plot_all
dev.off()


#### print model summaries ####

tab_model(m.bu.p)
summary(m.bu.p)
prior_summary(m.bu.p)
percpi

tab_model(m.bu.r)
summary(m.bu.r)
prior_summary(m.bu.r)
percri

tab_model(m.li.p)
summary(m.li.p)
prior_summary(m.li.p)
percpd

tab_model(m.li.r)
summary(m.li.r)
prior_summary(m.li.r)
percrd
