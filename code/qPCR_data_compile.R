#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)


#### import ####

# make a list of the file names
file_names <- list.files(path = "./data", pattern="[.]csv")

# select qPCR data
file_names <- file_names[str_detect(file_names, "Ex2_Group") == T]

# import files
for(i in 1:length(file_names)){
  
  dat <- read.csv(paste0("./data/", file_names[i])) %>%
    mutate(q_group = i)
  
  assign(paste0("dat", i), dat)
  
}


#### edit data ####

# combine data files
dat <- rbind(dat1, dat2, dat3, dat4, dat5, dat6, dat7, dat8, dat9, dat10, dat11, dat12, dat13, dat14, dat15, dat16, dat17, dat18, dat19, dat20, dat21, dat22)

# rename
colnames(dat) = c("well",	"sample",	"target",	"task",	"reporter",	"quencher",	"cycle", "cycle_mean",	"cycle_sd",	"quantity",	"quantity_mean",	"quantity_sd",	"auto_threshold",	"threshold",	"auto_baseline",	"baseline_start",	"baseline_end",	"comments",	"high_sd",	"no_amp",	"outlier_rg",	"exp_fail",	"q_group")


#### save data ####
write_csv(dat, "./intermediate-data/qPCR_data_compiled.csv")


#### re-save plant data with simpler name ####
plant_dat <- read_csv("data/plant_data_071617.csv")
write_csv(plant_dat, "intermediate-data/plant_data.csv")
