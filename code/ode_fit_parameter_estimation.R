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


#### plant parameters ####

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

# initial values
R0_n <- a_n_lo * 7
R0_p <- a_p_lo * 7
Q0_n <- Qmin_n
Q0_p <- Qmin_p
B0 <- 1e-3

# Combine values for model
y0_plant <- c(R_n = R0_n, R_p = R0_p, Q_n = Q0_n, Q_p = Q0_p, B = B0)


#### plant model ####

plant.model = function (t, yy, parms) { 
  
  # supply rates
  a_n = parms[1];
  a_p = parms[2];
  g = parms[3]
  m = parms[4]

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
         time = dpi)

# initiate slider for ggplot
manipulate(plot(1:5, cex=size), size = slider(0.5,10,step=0.5))

# time 
# times <- seq(0, max(dat5$dpi), length.out = 100)
times <- seq(0, 100)

# wrapper function
plant_wrapper <- function(g, m){
  
  out_low <- ode(y0_plant, times, plant.model, c(a_n = a_n_lo, a_p = a_p_lo, g = g, m = m)) %>%
    as_tibble() %>%
    mutate(nutrient = "low",
           across(!nutrient, as.double),
           nutrient = as.character(nutrient))
  
  out_N <- ode(y0_plant, times, plant.model, c(a_n = a_n_hi, a_p = a_p_lo, g = g, m = m)) %>%
    as_tibble() %>%
    mutate(nutrient = "N",
           across(!nutrient, as.double),
           nutrient = as.character(nutrient))
  
  out_P <- ode(y0_plant, times, plant.model, c(a_n = a_n_lo, a_p = a_p_hi, g = g, m = m)) %>%
    as_tibble() %>%
    mutate(nutrient = "P",
           across(!nutrient, as.double),
           nutrient = as.character(nutrient))
  
  out_NP <- ode(y0_plant, times, plant.model, c(a_n = a_n_hi, a_p = a_p_hi, g = g, m = m)) %>%
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



###################################################
# Objective function 2 (modCost) - try with raw data
###################################################

# Subset for needed columns
rdatCH=subset(subset(rdatC,highP==0),select=c("highN","dpi","FullMass"))

# Re-value highN
rdatCH$highN=as.factor(rdatCH$highN)
rdatCH$highN=revalue(rdatCH$highN,c("0"="Hl","1"="Hh"))

# Rename columns 
colnames(rdatCH)=c("name","time","val")

# Combine data
rdat3=subset(rdatCH,name=="Hh"|name=="Ph")

sse.model <- function(params0){ 
  q = params0[1];
  out = ode(y=y0h,times=seq(0,max(rdat3$time),length=100),func=ERHP.model,parms=c(q=q))
  return(modCost(model=out,obs=rdat3,y="val"))   
  #return(out) #for troubleshooting
}


###################################################
### Optimization: Estimate parameters
###################################################

params0=c(0.00004)  #initial guess

fit = modFit(sse.model,params0,lower=c(0),upper=c(1))
summary(fit)
deviance(fit)
fit$ssr
fit$ms

#8.70e-6

###################################################
### Obtain model prediction
###################################################

q = fit$par[1];
mod.pred = ode(y0h,times=times,func=ERHP.model,parms=c(q))


###################################################
### Plot Model and Data
###################################################

par(mfrow=c(2,1))
plot(mod.pred[,1],mod.pred[,4],ylim=c(0,0.5),type="l",col="red",xlab="Time (days)",ylab="Host mass (g)")
points(rdath$dpi,rdath$meanMass,col="red")
arrows(rdath$dpi,(rdath$meanMass-rdath$seMass),rdath$dpi,(rdath$meanMass+rdath$seMass),length=0.05,angle=90,code=3,col="red")
plot(mod.pred[,1],mod.pred[,5],ylim=c(0,4000),type="l",col="red",xlab="Time",ylab="Pathogen concentration*10^10")
points(rdath$dpi,rdath$meanPath,col="red")
arrows(rdath$dpi,(rdath$meanPath-rdath$sePath),rdath$dpi,(rdath$meanPath+rdath$sePath),length=0.05,angle=90,code=3,col="red")


###################################################
### Plot Model Predictions
###################################################

par(mfrow=c(2,1))
plot(mod.pred[,1],mod.pred[,2],type="l",col="red",xlab="Time (days)",ylab="Environmental nutrients (g)")
plot(mod.pred[,1],mod.pred[,3],type="l",col="red",xlab="Time (days)",ylab="Tissue nutrient concentration")


###################################################
### Figure for paper - come back to
###################################################


setwd("~/Google Drive/Within Host Coexistence/Figures")
pdf("virusParameters_ModFits_071917.pdf",width=9,height=12)
par(mfrow=c(4,2),cex.lab=1.5,cex.axis=1.5,cex.main=1.5,mar=c(5, 6, 2, 2))
plot(mod.pred[,1],mod.pred[,5],ylim=c(0,2000),type="l",lwd=2,col="red",xlab="",ylab="Pathogen\nconcentration (pg/g)",main="Low N")
points(rdatl$dpi,rdatl$meanPath,col="red")
arrows(rdatl$dpi,(rdatl$meanPath-rdatl$sePath),rdatl$dpi,(rdatl$meanPath+rdatl$sePath),length=0.05,angle=90,code=3,col="red")
mtext("A",side=3,line=1,adj=0,font=2)
plot(mod.pred[,1],mod.pred[,9],ylim=c(0,2000),type="l",lwd=2,col="blue",xlab="Time",ylab="",main="High N")
points(rdath$dpi,rdath$meanPath,col="blue")
arrows(rdath$dpi,(rdath$meanPath-rdath$sePath),rdath$dpi,(rdath$meanPath+rdath$sePath),length=0.05,angle=90,code=3,col="blue")
mtext("B",side=3,line=1,adj=0,font=2)
plot(mod.pred[,1],mod.pred[,4],ylim=c(0,0.6),type="l",lwd=2,col="red",xlab="",ylab="Host mass (g)")
points(rdatl$dpi,rdatl$meanMass,col="red")
arrows(rdatl$dpi,(rdatl$meanMass-rdatl$seMass),rdatl$dpi,(rdatl$meanMass+rdatl$seMass),length=0.05,angle=90,code=3,col="red")
mtext("C",side=3,line=1,adj=0,font=2)
plot(mod.pred[,1],mod.pred[,8],ylim=c(0,0.6),type="l",lwd=2,col="blue",xlab="",ylab="")
points(rdath$dpi,rdath$meanMass,col="blue")
arrows(rdath$dpi,(rdath$meanMass-rdath$seMass),rdath$dpi,(rdath$meanMass+rdath$seMass),length=0.05,angle=90,code=3,col="blue")
mtext("D",side=3,line=1,adj=0,font=2)
plot(mod.pred[,1],mod.pred[,3],type="l",lwd=2,col="red",xlab="",ylab="Tissue nutrient\nconcentration")
mtext("E",side=3,line=1,adj=0,font=2)
plot(mod.pred[,1],mod.pred[,7],type="l",lwd=2,col="blue",xlab="",ylab="")
mtext("F",side=3,line=1,adj=0,font=2)
plot(mod.pred[,1],mod.pred[,2],type="l",lwd=2,col="red",xlab="Days post inoculation",ylab="Environmental\nnutrients (g)")
mtext("G",side=3,line=1,adj=0,font=2)
plot(mod.pred[,1],mod.pred[,6],type="l",lwd=2,col="blue",xlab="Days post inoculation",ylab="")
mtext("H",side=3,line=1,adj=0,font=2)
dev.off()
