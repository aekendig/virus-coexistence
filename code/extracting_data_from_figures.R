#### Goal: extract data from figures

# reference: Emmanuel Jjunju, https://www.r-bloggers.com/digitizing-jpeg-graphs-in-r/


#### Set-up ####

# clear all existing data
rm(list=ls())

# libraries
library(jpeg) # to use jpeg images in R
library(zoo)
library(tidyverse)

# import dray/wet weight conversion
chn <- read_csv("data/CHNweights_022514.csv")

# digitize functions
ReadAndCal = function(fname){
  ReadImg(fname)
  calpoints <- locator(n=4,type='p',pch=4,col='blue',lwd=2)
  return(calpoints)
}

ReadImg = function(fname){
  img <- readJPEG(fname)
  op <- par(mar=c(0,0,0,0))
  on.exit(par(op))
  plot.new()
  rasterImage(img,0,0,1,1)
}

DigitData = function(col='red',type='p',...){
  type <- ifelse(type=='b','o',type)
  type <- ifelse(type%in%c('l','o','p'),type,'p')
  locator(type=type,col=col,...)
}

Calibrate = function(data,calpoints,x1,x2,y1,y2){
  x 		<- calpoints$x[c(1,2)]
  y 		<- calpoints$y[c(3,4)]
  
  cx <- lm(formula = c(x1,x2) ~ c(x))$coeff
  cy <- lm(formula = c(y1,y2) ~ c(y))$coeff
  
  data$x <- data$x*cx[2]+cx[1]
  data$y <- data$y*cy[2]+cy[1]
  
  return(as.data.frame(data))
}


#### Steps ####

# ReadAndCal opens the jpeg in a plotting window and lets you define points on the x and y axes. You must start by clicking on the left-most x-axis point, then the right-most axis point, followed by the lower y-axis point and finally the upper y-axis point. You don’t need to choose the end points of the axis, only two points on the axis that you know the x or y value for. As you click on each of the 4 points, the coordinates are saved in the object cal.

# DigitData returns you to the figure window, and now you can click on each of the data points you’re interested in retrieving values for. The function will place a dot (colored red in this case) over each point you click on, and the raw x,y coordinates of that point will be saved to the data.points list. When you’re finished clicking points, you need to hit stop/Finish or right-click to stop the data point collection.

# Calibrate converts those raw x,y coordinates into the same scale as the original graph. Feed the function your data.point list, the ReadAndCal list that contains your 4 control points from the first step, and then 4 numeric values that represent the 4 original points you clicked on the x and y axes. These values should be in the original scale of the figure (i.e. read the values off the graph’s tick marks).


#### dry:wet conversion ####

wt_dat <- chn %>%
  mutate(wet_dry = WetWt/DryWt,
         dry_wet = DryWt/WetWt) %>%
  summarise(wet_dry = mean(wet_dry),
            dry_wet = mean(dry_wet))


#### Digitize figures ####

# Mattsson et al. 1991 4A
(cal_ma91_4a = ReadAndCal("data/Mattsson_etal_1991_4A.jpg"))
(data_ma91_4a = DigitData(col = 'red'))
df_ma91_4a = Calibrate(data_ma91_4a, cal_ma91_4a, 0, 0.2, 0, 100)
df_ma91_4a <- df_ma91_4a %>%
  rename(RA = x, Vmax = y) %>%
  mutate(variety = c("Laevigatum", "Golf", "Mette"),
         RA = round(RA, 2),
         Vmax = round(Vmax, 2))

# Mattsson et al. 1991 2B
(cal_ma91_2b = ReadAndCal("data/Mattsson_etal_1991_2B.jpg"))
(data_ma91_2b = DigitData(col = 'red'))
df_ma91_2b = Calibrate(data_ma91_2b, cal_ma91_2b, 0, 0.2, 0, 0.7)
df_ma91_2b <- df_ma91_2b %>%
  rename(RA = x, rootTot = y) %>%
  mutate(variety = c("Laevigatum", "Golf", "Mette"),
         RA = round(RA, 2),
         rootTot = round(rootTot, 2))

# Vmax for total biomass
# convert from mg to g
# use root:full from my experiment
df_ma91_vmax <- df_ma91_4a %>%
  full_join(df_ma91_2b) %>%
  mutate(VmaxTot = Vmax * rootTot, # mg N per g dry wt per day
         VmaxTot_g = VmaxTot * 0.001, # g N per g dry wt per day
         VmaxTot_g2 = Vmax * 0.07 * 0.001, # g N per g wet wt per day
         VmaxTot_g3 = VmaxTot_g * wt_dat$dry_wet) # g N per g wet wt per day

# save data
write_csv(df_ma91_vmax, "intermediate-data/Mattsson_etal_1991_vmax.csv")

# average across varieties
mean(df_ma91_vmax$VmaxTot_g)
mean(df_ma91_vmax$VmaxTot_g2)
mean(df_ma91_vmax$VmaxTot_g3)

# Mattsson et al. 1991 4B
(cal_ma91_4b = ReadAndCal("data/Mattsson_etal_1991_4B.jpg"))
(data_ma91_4b = DigitData(col = 'red'))
df_ma91_4b = Calibrate(data_ma91_4b, cal_ma91_4b, 0, 0.2, 0, 60)
df_ma91_4b <- df_ma91_4b %>%
  rename(RA = x, Km = y) %>%
  mutate(variety = c("Laevigatum", "Golf", "Mette"),
         RA = round(RA, 2),
         Km = round(Km, 2),
         Km_g_pot = Km * 0.0140067 / 1000 * 0.1233)

# save data
write_csv(df_ma91_4b, "intermediate-data/Mattsson_etal_1991_4b.csv")

# half-saturation constant
mean(df_ma91_4b$Km)
mean(df_ma91_4b$Km_g_pot)

# Mattsson et al. 1991 5A (Qmin for N)
(cal_ma91_5a = ReadAndCal("data/Mattsson_etal_1991_5A.jpg"))
(data_ma91_5a = DigitData(col = 'red'))
df_ma91_5a = Calibrate(data_ma91_5a, cal_ma91_5a, 0, 60, 0, 0.2)
df_ma91_5a <- df_ma91_5a %>%
  rename(N_conc_mg = x, plant_rg = y) %>%
  mutate(N_conc_g = N_conc_mg * 0.001,
         Q_recip = 1 / N_conc_g)

# save data
write_csv(df_ma91_5a, "intermediate-data/Mattsson_etal_1991_5A.csv")

# look at relationship
df_ma91_5a %>%
  ggplot(aes(Q_recip, plant_rg)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x)

# fit regression
mod_ma91_5a <- lm(plant_rg ~ Q_recip, data = df_ma91_5a)

# Qmin (when m is negligible)
# multiple by dry:wet mass to get in units of wet weight
-1 * coef(mod_ma91_5a)[2]/coef(mod_ma91_5a)[1] * wt_dat$dry_wet

# Ullrich-Eberius et al. 1981 6 (Qmin for P)
(cal_ue81_6 = ReadAndCal("data/Ullrich-Eberius_etal_1981_6.jpg"))
(data_ue81_6 = DigitData(col = 'red'))
df_ue81_6 = Calibrate(data_ue81_6, cal_ue81_6, 0, 16, 0, 30)
# y value is Qmin estimate

# Ullrich-Eberius et al. 1984
# half-sat constant = 7.9 uM (umol/L) P
# 1 umol = 10^-6 mol
# atomic mass = 30.97376 g/mol 
# 1 L water = 1 kg water
7.9 * 10^-6 * 30.97376 * 1 * 0.1233
# umol/L * (mol/umol) * (g/mol) * (L/kg) * kg
