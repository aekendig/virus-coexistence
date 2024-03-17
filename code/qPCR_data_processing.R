#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse) # version 2.0.0

# import data
qdat <- read_csv("intermediate-data/qPCR_data_compiled.csv") # qPCR data
edat <- read_csv("intermediate-data/plant_data.csv") # experiment data


#### edit data ####

# remove summarizing columns and empty wells
qdat2 <- qdat %>%
  select(-c(cycle_mean, cycle_sd, quantity_mean, quantity_sd, comments, high_sd, no_amp, outlier_rg, exp_fail)) %>%
  filter(!is.na(sample))

# examine sample names, tasks
unique(qdat2$sample)
unique(qdat2$task)
unique(qdat2$cycle)

# check notes
unique(edat$expt_notes)
filter(edat, expt_notes == "Planned to sample at T1, but sampled at T5. May be labelled as L PAV T1 S 1.") %>%
  data.frame()
unique(edat$extraction_notes)
# mis-labelling is handled later

# make cycle numeric
# update tasks
# fix label errors (manually checked the group-specific one)
qdat3 <- qdat2 %>%
  mutate(cycle = as.numeric(cycle),
         task = case_when(substr(sample, 1, 1) == "S" ~ "sample",
                          substr(sample, 1, 1) == "N" ~ "control",
                          substr(sample, 1, 1) == "P" & target == "RPV" ~ "nonTargetStandard",
                          substr(sample, 1, 1) == "R" & target == "PAV" ~ "nonTargetStandard",
                          task == "STANDARD" ~ "standard",
                          task == "TEST" ~ "standardTest"),
         sample = case_when(sample == "S1NPP2S1" ~ "S1 NPP2S1",
                            sample == "S1LR4S2" ~ "S1 LR4S2",
                            sample == "S35 NPR4S1" ~ "S3R NPR4S1",
                            sample == "S3 LR3I2" ~ "S2 LR3I2",
                            sample == "S2 LR2S2" & q_group == 7 ~"S1 LR2S2",
                            TRUE ~ sample),
         tube_label = case_when(task == "sample" ~ sample,
                                TRUE ~ NA_character_))

unique(qdat3$task)

# number of tube labels
length(unique(qdat3$tube_label))
length(unique(edat$tube_label))

# check for NA's
sum(is.na(qdat3$tube_label))
sum(is.na(edat$tube_label))

# examine edat NA's
filter(edat, is.na(tube_label)) %>%
  select(expt_notes, extraction_notes, extraction_date)
# 18 samples not extracted, all have justifications

# remove samples with errors from edat data
edat2 <- filter(edat, !is.na(tube_label))

# examine qdat NA's
filter(qdat3, is.na(tube_label)) %>%
  select(task) %>%
  unique()
# none are samples

# overlap in tube numbers
nrow(filter(qdat3, task == "sample" & !(tube_label %in% edat2$tube_label)))
nrow(filter(edat2, !(tube_label %in% qdat3$tube_label)))

# add edat information to samples
# create combine set and replicate column
qdat4 <- qdat3 %>%
  left_join(edat2 %>%
              select(-sample)) %>%
  mutate(set_rep = paste(set, replicate, sep = "_"))


#### error check ####

# see how many were added
nrow(edat2)
length(unique(filter(qdat4, !is.na(freezer_bag_ID))$tube_label))
# all

# check for missing edat data
qdat4 %>%
  filter(task == "sample" & (is.na(nutrient) | is.na(inoculation) | is.na(invasion))) %>% 
  select(sample) %>% 
  unique() # none

# check for duplicate qPCR labels
qdat4 %>%
  filter(task == "sample") %>%
  group_by(tube_label, q_group) %>%
  count() %>%
  filter(n != 6) # none duplicated within qPCR group

qdat4 %>%
  filter(task == "sample") %>%
  select(tube_label, q_group) %>%
  unique() %>%
  mutate(dup = duplicated(tube_label)) %>%
  filter(dup == T) %>%
  select(tube_label) %>%
  inner_join(qdat4 %>% 
               select(tube_label, q_group) %>%
               unique()) %>%
  arrange(tube_label) %>%
  data.frame() 
# all are in groups 21 and 22 for their second run

# summarize sample number
qdat4 %>%
  filter(task == "sample") %>%
  group_by(nutrient, inoculation, invasion, time, target) %>%
  summarize(reps = length(unique(set_rep))) %>%
  ggplot(aes(x = reps)) +
  geom_histogram(binwidth = 1) # 4-6

# examine non-samples
qdat4 %>%
  filter(task != "sample") %>%
  group_by(task, target, quantity) %>%
  summarize(reps = n()) %>%
  print.data.frame()

# select standards with different quantities
qdat4 %>%
  filter(task == "standardTest") %>%
  select(task, q_group, sample, target, quantity) %>%
  data.frame()
# all in group 1


#### standard curves ####

# group 1 standards
q_group_1_std <- qdat4 %>%
  filter(q_group == 1 & task %in% c("standard", "standardTest")) %>%
  mutate(assigned_quantity = as.numeric(substring(sample, 5)))

# visualize group 1 standards
ggplot(q_group_1_std, aes(x = log10(assigned_quantity), y = cycle, color = task)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ target)
# PAV tests are higher quantities than standards
# RPV tests and standards are overlapping

# wells of group 1 standards
q_group_1_std %>%
  group_by(task, target) %>%
  summarise(wells = paste(well, collapse = ","))
# tests on right side of plate
# location didn't affect RPV standards

# PAV standards function
stdPAVfun <- function(dat){
  
  # extract min, max, and curve characteristics
  stdPAVdf <- dat %>%
    subset(!is.na(cycle) & target == "PAV" & task == "standard") %>%
    group_by(q_group) %>%
    summarise(PAVmin = min(quantity), 
              PAVminRep = sum(quantity == min(quantity)),
              PAVmax = max(quantity), 
              PAVslope = lm(cycle ~ log10(quantity))$coefficients[2], 
              PAVint=lm(cycle ~ log10(quantity))$coefficients[1])
  
  # calculate efficiency
  stdPAVdf$PAVefficiency <- 100 * (10^(1/abs(stdPAVdf$PAVslope)) - 1)
  
  # plot
  # stdPAVplot <- dat %>%
  #  subset(!is.na(cycle) & target == "PAV" & task == "standard") %>%
  #  ggplot(aes(x = log10(quantity), y = cycle)) +
  #  geom_point(size = 2) +
  #  xlab("Quantity (log10)") +
  #  ylab("Detection cycle") +
  #  geom_smooth(method = lm, se = F)
  
  # print(stdPAVplot)
  
  return(stdPAVdf)
} 

# RPV standards function
stdRPVfun <- function(dat){
  
  # extract min, max, and curve characteristics
  stdRPVdf <- dat %>%
    subset(!is.na(cycle) & target == "RPV" & task == "standard") %>%
    group_by(q_group) %>%
    summarise(RPVmin = min(quantity), 
              RPVminRep = sum(quantity == min(quantity)),
              RPVmax = max(quantity), 
              RPVslope = lm(cycle ~ log10(quantity))$coefficients[2], 
              RPVint=lm(cycle ~ log10(quantity))$coefficients[1])
  
  # calculate efficiency
  stdRPVdf$RPVefficiency <- 100 * (10^(1/abs(stdRPVdf$RPVslope)) - 1)
  
  # plot
  # stdRPVplot <- dat %>%
  #  subset(!is.na(cycle) & target == "RPV" & task == "standard") %>%
  #  ggplot(aes(x = log10(quantity), y = cycle)) +
  #  geom_point(size = 2) +
  #   xlab("Quantity (log10)") +
  #   ylab("Detection cycle") +
  #   geom_smooth(method = lm, se = F)
  
  # print(stdRPVplot)
   
  return(stdRPVdf)
}

# examine standard curve figures 
# for(i in unique(qdat4$q_group)){
#   
#   stdPAVfun(filter(qdat4, q_group == i))
#   stdRPVfun(filter(qdat4, q_group == i))
#   
# }
# look pretty good - outliers aren't pulling the curve, most are linear

# contamination function
confun <- function(dat){
  
  # identify highest qPCR value in negative controls
  # add as columns to rest of data
  fulldf <- dat %>%
    filter(task == "control"|task == "nonTargetStandard") %>%
    group_by(q_group, target) %>%
    summarise(contamination = min(cycle, na.rm = T)) %>%
    ungroup() %>%
    spread(target, contamination) %>%
    rename(PAVcont = PAV, RPVcont = RPV) %>%
    merge(dat, ., all = T)
  
  # indicate whether standard should be removed because it's less 
  # concentrated than contamination
  stddf <- fulldf %>%
    group_by(q_group) %>%
    filter(task == "standard") %>%
    mutate(stdRem = case_when(target == "PAV" & cycle > unique(PAVcont) ~ 1,
                              target == "RPV" & cycle > unique(RPVcont) ~ 1,
                              TRUE ~ 0)) %>%
    ungroup() %>%
    group_by(q_group, target) %>%
    summarise(removed = sum(stdRem, na.rm = T)) %>%
    ungroup() %>%
    spread(target, removed) %>%
    rename(PAVstdRem = PAV, RPVstdRem = RPV) %>%
    merge(fulldf, ., all = T)
  
  # remove those standards
  # re-calculate standard curve
  PAVdf <- fulldf %>%
    filter(cycle < PAVcont) %>%
    stdPAVfun
  
  RPVdf <- fulldf %>%
    filter(cycle < RPVcont) %>%
    stdRPVfun
  
  # combine new dataframe
  findf <- stddf %>%
    merge(PAVdf, all = T) %>%
    merge(RPVdf, all = T)
  
  return(findf)
}

# identify contamination
# remove standards and re-calculate standard curve
qdat5 <- qdat4 %>%
  confun()
# warnings: Inf returned when no contamination
# this is okay because all cyclces < Inf

# check that confun worked
head(qdat5)
qdat5 %>%
  group_by(PAVstdRem, RPVstdRem) %>%
  summarise(groups = length(unique(q_group)))

# look at high concentration contamination
qdat5 %>%
  filter(RPVstdRem == 6 & PAVstdRem == 0 & target == "RPV") %>%
  ggplot(aes(cycle, log10(quantity), color = task)) +
  geom_point(alpha = 0.5) +
  geom_vline(aes(xintercept = RPVcont), size = 0.2)
# some samples below, but pretty low concentration

qdat5 %>%
  filter(RPVstdRem == 6 & PAVstdRem == 1 & target == "RPV") %>%
  ggplot(aes(cycle, log10(quantity), color = task)) +
  geom_point(alpha = 0.5) +
  geom_vline(aes(xintercept = RPVcont), size = 0.2)
# no samples below

qdat5 %>%
  filter(RPVstdRem == 17 & PAVstdRem == 0 & target == "RPV") %>%
  ggplot(aes(cycle, log10(quantity), color = task)) +
  geom_point(alpha = 0.5) +
  geom_vline(aes(xintercept = RPVcont), size = 0.2)
# all are below, group 22 (contains duplicates)

qdat5 %>%
  filter(q_group == 22 & task %in% c("control", "nonTargetStandard") & target == "RPV")
# two PAV samples have very high contamination

# modify samples with new standard curve
qdat6 <- qdat5 %>%
  mutate(quant_adj = case_when(task != "standard" & target == "PAV" ~ 10 ^ ((cycle - PAVint) / PAVslope),
                               task != "standard" & target == "RPV" ~ 10 ^ ((cycle - RPVint) / RPVslope),
                               task == "standard" ~ quantity))

# check that new standard curve conversion worked
qdat6 %>%
  ggplot(aes(x = log10(quantity), y = log10(quant_adj))) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0)


#### create sample dataset ####

samp <- qdat6 %>%
  filter(task == "sample") %>%
  as_tibble()

# no std curve info
samp %>%
  filter(is.na(PAVslope) | is.na(RPVslope)) %>%
  select(q_group) %>% unique()
# just 22 (high RPV contamination)

# minimum values
min(samp$RPVmin, na.rm = T)
min(samp$PAVmin, na.rm = T)

# check notes
unique(samp$extraction_notes)
filter(samp, extraction_notes == "May be PP4S2. Both tubes were labelled PP4S2")
filter(samp, extraction_notes == "May be PP4I2. Both tubes were labelled PP4S2")
# both contain PAV and RPV above the minimum detection

# identify which samples will be removed:
# standard curve efficiency outside of boundaries (80-120%)
# quantity greater than 1e3, but less than standard curve minimum (removes contamination)
# mis-labelled samples that can't be reconciled
samp2 <- samp %>%
  mutate(remove = case_when(is.na(PAVslope) | is.na(RPVslope) ~ 1, # 120 (group 22)
                            target == "RPV" & (round(RPVefficiency) < 80 | round(RPVefficiency) > 120) ~ 1, # 60 (group 4)
                            target == "PAV" & (round(PAVefficiency) < 80 | round(PAVefficiency) > 120) ~ 1, # 0
                            target == "RPV" & (quant_adj >= 1e3 & quant_adj < RPVmin) ~ 1, # 15 (groups 4 and 6)
                            target == "PAV" & (quant_adj >= 1e3 & quant_adj < PAVmin) ~ 1, # 23 (group 21)
                            sample %in% c("S1 PP4I2", "S1 PP4S2") ~ 1,
                            TRUE ~ 0))
sum(samp2$remove) 
# 227

# efficiencies of std curves
samp2 %>%
  select(q_group, RPVefficiency, PAVefficiency) %>%
  unique() %>%
  arrange(q_group) %>%
  data.frame()
# RPV group 4: 66%
# PAV group 11: 82%
# PAV group 21: 83%

# identify duplicate samples and select ones with detection or lower variance
samp_d <- samp2 %>%
  filter(remove == 0) %>% # only looking at ones that won't be removed
  select(q_group, target, sample) %>%
  distinct() %>%
  mutate(
    dup = duplicated(select(., target, sample))
  ) %>%
  filter(dup == T) %>% # select duplicates
  select(-c(q_group, dup)) %>%
  semi_join(filter(samp2, remove == 0), .) %>% # all samples in duplicate set
  mutate(detect = case_when(target == "RPV" & quant_adj >= RPVmin ~ 1, # detected?
                            target == "PAV" & quant_adj >= PAVmin ~ 1,
                            TRUE ~ 0)) %>%
  group_by(q_group, sample, target) %>%
  summarise(quant_var = sd(quant_adj, na.rm = T), # summary statistics within group
            quant_mean = mean(quant_adj, na.rm = T),
            detect = sum(detect) > 0) %>%
  ungroup() %>%
  group_by(sample, target) %>%
  mutate(min_var_group = case_when(quant_var == min(quant_var, na.rm = T) ~ 1,
                                   TRUE ~ 0),
         n_q_group = length(unique(q_group)),
         n_min = sum(min_var_group),
         n_det = sum(detect),
         sample_q_tar = paste(sample, q_group, target, sep = "_")) %>%
  ungroup() %>%
  arrange(sample, target)
# okay to ignore warnings

data.frame(samp_d)
# 4 not detected

ggplot(samp_d, aes(x = quant_mean, y = quant_var, color = target)) +
  geom_point()
cor.test(~ quant_mean + quant_var, data = samp_d)
# the high variance ones are also the high mean ones

# decided not to do below because:
# the undetected samples don't have detected replacements
# by removing the high variance samples, we're removing the high concentration samples

# indicate which samples should be removed
# not detected when another group was
# higher than minimum variance when multiple groups were detected
# samp_d_r <- samp_d %>%
#   mutate(remove = case_when(detect == F & n_det > 0 ~ 1,
#                             detect == T & n_det > 1 & min_var_group == 0 ~ 1,
#                             TRUE ~ 0)) %>%
#   filter(remove == 1)

# save full data
write_csv(samp2, "intermediate-data/qPCR_expt_data_compiled.csv")


#### clean for analyses ####

# NA values
samp2 %>%
  filter(is.na(quant_adj) & remove == 0) %>%
  select(cycle, quantity) %>%
  unique()

# no analysis column (from edat)
filter(samp2, no_analysis == 1) %>%
  select(expt_notes, extraction_notes) %>%
  unique()
unique(samp2$no_analysis)

# maximum values
unique(samp2$RPVmax)
samp2 %>%
  group_by(RPVmax) %>%
  count()
# most are 1e8
unique(samp2$PAVmax)

# time and dpi conversion
times = tibble(time = c(1,2,3,4,5),
               dpiI = c(5,8,12,16,19)) %>%
  mutate(dpiR = dpiI + 12,
         dpp = dpiR + 11,)

# clean data for analyses
samp3 <- samp2 %>%
  filter(remove == 0 & is.na(no_analysis)) %>% # removes issues identified above
  group_by(set, nutrient, inoculation, time, invasion, replicate, set_rep, sample, tube_label, extraction_mass.mg, shoot_mass.g, expt_notes, extraction_notes, target, q_group, PAVint, PAVslope, PAVmin, PAVmax, RPVint, RPVslope, RPVmin, RPVmax) %>%
  summarise(tech_cycle = mean(cycle, na.rm = T)) %>%
  ungroup() %>%
  mutate(quant = case_when(target == "PAV" ~ 10 ^ ((tech_cycle - PAVint) / PAVslope),
                           target == "RPV" ~ 10 ^ ((tech_cycle - RPVint) / RPVslope)),
         quant_adj = case_when(target == "PAV" & quant < PAVmin ~ 0, # below standard curve
                               target == "RPV" & quant < RPVmin ~ 0,
                               is.na(quant) ~ 0, # not detected
                               TRUE ~ quant)) %>%
  group_by(set, nutrient, inoculation, time, invasion, replicate, set_rep, sample, tube_label, extraction_mass.mg, shoot_mass.g, expt_notes, extraction_notes, target) %>%
  summarise(quant_adj = mean(quant_adj)) %>% # average across qPCR groups
  ungroup() %>%
  mutate(quant.ul = quant_adj / 7,
         quant.mg = quant.ul * 50 / extraction_mass.mg,
         highN = ifelse(nutrient %in% c("N", "NP"), 1, 0),
         highP = ifelse(nutrient %in% c("P", "NP"), 1, 0),
         second_inoculation = case_when(inoculation == "PAV" & invasion == "I" ~ "RPV",
                                        inoculation == "RPV" & invasion == "I" ~ "PAV",
                                        TRUE ~ NA_character_)) %>%
  left_join(times) %>%
  rename(first_inoculation = inoculation) %>%
  filter(!(target == "PAV" & quant_adj > 1e7) &
           !(target == "RPV" & quant_adj > 1e8)) # removes 30 RPV samples

# save data
write_csv(samp3, "intermediate-data/qPCR_expt_data_cleaned.csv")

