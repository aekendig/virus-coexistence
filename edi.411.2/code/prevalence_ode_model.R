## Goal: ordinary differential equation model of virus prevalence based on model estimated transmission values

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse) # version used: 1.2.1
library(cowplot) # version used: 0.9.4
library(brms)  # version used: 2.7.0
library(deSolve) # version used: 1.21

# import data
ptran <- read_csv("./output/pav_transmission_values_same_nuts.csv")
rtran <- read_csv("./output/rpv_transmission_values_same_nuts.csv")

#### edit data ####

# calculate parameters
ptran2 <- ptran %>%
  select(mechanisms, nutrient, parameter, value) %>%
  spread(key = parameter, value = value) %>%
  rename(BP = B,
         qP = q)

rtran2 <- rtran %>%
  select(mechanisms, nutrient, parameter, value) %>%
  spread(key = parameter, value = value) %>%
  rename(BR = B,
         qR = q)

# combine data
tran <- full_join(ptran2, rtran2)


#### model ####

simmod = function(parm_dat){
  
  parms = list(BP = parm_dat$BP, qP = parm_dat$qP, BR = parm_dat$BR, qR = parm_dat$qR, Sinit = Sinit, Pinit = Pinit, Rinit = Rinit, Cinit = Cinit, N = Ninit, simtime = simtime)
  
  mymodel = with(as.list(parms), function(t, x, parms){
    
    S = x["S"]
    P = x["P"]
    R = x["R"]
    C = x["C"]
    
    Sdot = -(BP*P + BR*R + qP*BP*(1 - qR*BR)*C + qR*BR*(1 - qP*BP)*C + qP*BP*qR*BR*C)*S/N 
    Pdot = BP*(P + qP*(1 - qR*BR)*C)*S/N - BR*(R + qR*C)*P/N
    Rdot = BR*(R + qR*(1 - qP*BP)*C)*S/N - BP*(P + qP*C)*R/N
    Cdot = BP*(P + qP*C)*R/N + BR*(R + qR*C)*P/N + qP*BP*qR*BR*C*S/N
    
    list(c(Sdot, Pdot, Rdot, Cdot))
  })
  
  xstart = c(S = Sinit, P = Pinit, R = Rinit, C = Cinit)
  
  times = seq(0, simtime, length = simtime)
  
  out = as.data.frame(lsoda(xstart, times, mymodel, parms, hmax=20))
  
  return(out)
}


#### constants ####

Ntot = 4000
Cinit = 0
days = 120
trans_time = 6
simtime = days / trans_time
time_dat <- tibble(time = 1:simtime)


#### PAV only model ####

Sinit = Ntot - 1
Pinit = 1
Rinit = 0
Ninit = Sinit + Pinit + Rinit + Cinit

# create columns
tranP <- expand_grid(tran, time_dat) %>%
  mutate(virus = "PAV",
         S = NA,
         P = NA,
         R = NA,
         C = NA)

# add output to each parameter set
pdf("./output/prevalence_ode_pav_simulations.pdf")
for(i in 1:nrow(tran)){
  
  mod <- simmod(tran[i, ])
  
  tranP[(simtime*(i-1) + 1):(simtime*i), c('S','P', 'R', 'C')] <- mod[ , c('S','P', 'R', 'C')]
  
  print(ggplot(mod, aes(x = time)) +
          geom_line(aes(y = S), color = "black") +
          geom_line(aes(y = P), color = "red") +
          geom_line(aes(y = R), color = "blue") +
          geom_line(aes(y = C), color = "purple") +
          theme_bw() +
          ggtitle(paste(tran[i,"nutrient"], tran[i, "mechanisms"], sep = ", ")))
}
dev.off()
tranP


#### RPV only model ####

Sinit = Ntot - 1
Pinit = 0
Rinit = 1
Ninit = Sinit + Pinit + Rinit + Cinit

# create columns
tranR <- expand_grid(tran, time_dat) %>%
  mutate(virus = "RPV",
         S = NA,
         P = NA,
         R = NA,
         C = NA)

# add output to each parameter set
pdf("./output/prevalence_ode_rpv_simulations.pdf")
for(i in 1:nrow(tran)){
  
  mod <- simmod(tran[i, ])
  
  tranR[(simtime*(i-1) + 1):(simtime*i), c('S','P', 'R', 'C')] <- mod[ , c('S','P', 'R', 'C')]
  
  print(ggplot(mod, aes(x = time)) +
          geom_line(aes(y = S), color = "black") +
          geom_line(aes(y = P), color = "red") +
          geom_line(aes(y = R), color = "blue") +
          geom_line(aes(y = C), color = "purple") +
          theme_bw() +
          ggtitle(paste(tran[i,"nutrient"], tran[i, "mechanisms"], sep = ", ")))
}
dev.off()
tranR


#### PAV and RPV model ####

Sinit = Ntot - 2
Pinit = 1
Rinit = 1
Ninit = Sinit + Pinit + Rinit + Cinit

# create columns
tranC <- expand_grid(tran, time_dat) %>%
  mutate(virus = "both",
         S = NA,
         P = NA,
         R = NA,
         C = NA)

# add output to each parameter set
pdf("./output/prevalence_ode_coinfection_simulations.pdf")
for(i in 1:nrow(tran)){
  
  mod <- simmod(tran[i, ])
  
  tranC[(simtime*(i-1) + 1):(simtime*i), c('S','P', 'R', 'C')] <- mod[ , c('S','P', 'R', 'C')]
  
  print(ggplot(mod, aes(x = time)) +
          geom_line(aes(y = S), color = "black") +
          geom_line(aes(y = P), color = "red") +
          geom_line(aes(y = R), color = "blue") +
          geom_line(aes(y = C), color = "purple") +
          theme_bw() +
          ggtitle(paste(tran[i,"nutrient"], tran[i, "mechanisms"], sep = ", ")))
}
dev.off()
tranC


#### format data ####

# data for comparing simulations with one or two viruses
prev_virus <- full_join(tranP, tranR) %>%
  full_join(tranC) %>%
  mutate(PAV = (P + C) / Ntot,
         RPV = (R + C) / Ntot) %>%
  select(mechanisms, virus, nutrient, time, PAV, RPV) %>%
  gather(key = infection, value = prev, -c(mechanisms, virus, nutrient, time)) %>%
  filter(!(virus == "PAV" & infection == "RPV") & !(virus == "RPV" & infection == "PAV")) %>%
  mutate(virus = recode(virus, PAV = "single", RPV = "single")) %>%
  spread(key = virus, value = prev) %>% 
  mutate(diff = both - single,
         time = time * trans_time,
         nutrient = fct_relevel(nutrient, "low", "N", "P"),
         mechanisms = recode(mechanisms, all = "all processes", density = "virus\ndensity-dependent", "non-density" = "virus\ndensity-independent"))

# data for looking at raw prevalence instead of diff
raw_prev_virus <- prev_virus %>%
  select(-diff) %>%
  rename(virus = infection) %>%
  gather(key = "infection", value = "prevalence", -c(mechanisms:virus))


#### figure settings ####

# palettes
col_pal = c("black", "darkgoldenrod2", "dodgerblue1", "palegreen4")
line_pal = c("solid", "dashed")
shape_pal = c(21, 22, 24)

# text sizes
sm_txt = 6
lg_txt = 8


#### figures ####

# text
text_dat <- tibble(mechanisms = rep(unique(prev_virus$mechanisms), 2),
                   infection = rep(c("PAV", "RPV"), each = 3),
                   text = letters[1:6]) %>%
  mutate(time = 0,
         diff = rep(c(0.63, 0.16), each = 3))

# effects of virus interactions on prevalence
simfig <- ggplot(prev_virus, aes(x = time, y = diff)) +
  geom_hline(yintercept = 0, color = 'gray', linetype = "dashed") +
  geom_line(aes(color = nutrient)) +
  geom_text(data = text_dat, aes(label = text), size = 3) +
  facet_grid(infection ~ mechanisms, scales = "free_y") +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        plot.title = element_text(color = "black", size = lg_txt, hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(color = "black", size = lg_txt),
        strip.background = element_blank(),
        panel.spacing = unit(0, "lines"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(color = "black", size = lg_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.background = element_blank(),
        legend.key = element_rect(color = "white", size = 0.5),
        legend.key.size = unit(1.5, "lines"),
        legend.key.width = unit(1, "cm"),
        legend.spacing.x = unit(0.001, "mm"),
        legend.justification = c(0.5, 0.5),
        legend.box.margin = margin(-10, -10, -10, -10)) +
  scale_color_manual(values = col_pal, name = "Nutrient") +
  xlab("Days") +
  ylab("Effect of among-virus interactions on infection prevalence")

pdf("./output/figure_4_simulation.pdf", width = 6, height = 4)
simfig
dev.off()

# text for raw data
raw_text_dat <- tibble(virus = c("PAV", "RPV"),
                   text = letters[1:2]) %>%
  mutate(time = 0,
         prevalence = 1)

# raw prevalence values
rawfig <- ggplot(filter(raw_prev_virus, mechanisms == "all processes" & infection == "both"), aes(x = time, y = prevalence)) +
  geom_line(aes(color = nutrient), alpha = 0.7) +
  geom_text(data = raw_text_dat, aes(label = text), size = 3) +
  facet_wrap(~virus) +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        plot.title = element_text(color = "black", size = lg_txt, hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(color = "black", size = lg_txt),
        strip.background = element_blank(),
        panel.spacing = unit(0, "lines"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(color = "black", size = lg_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.background = element_blank(),
        legend.key = element_rect(color = "white", size = 0.5),
        legend.key.size = unit(1.5, "lines"),
        legend.key.width = unit(1, "cm"),
        legend.spacing.x = unit(0.001, "mm"),
        legend.justification = c(0.5, 0.5),
        legend.box.margin = margin(-10, -10, -10, -10)) +
  scale_color_manual(values = col_pal, name = "Nutrient") +
  xlab("Days") +
  ylab("Predicted infection prevalence")

pdf("./output/figure_S4_simulation.pdf", width = 4, height = 2.5)
rawfig
dev.off()