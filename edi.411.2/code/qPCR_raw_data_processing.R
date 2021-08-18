## Goal: combine raw qPCR data into a single dataset


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse) # version used: 1.2.1

# import qPCR data
qdat <- read_csv("./edi.411.2/data/source_plant_qPCR_data.csv")
#str(qlist)

# import experiment and molecular data
emdat <- read.csv("./edi.411.2/data/sample_exp_molc_data.csv", header = T, stringsAsFactors = F) %>%
  as_tibble()


#### edit full dataset ####

# remove summarizing column
qdat <- qdat %>%
  select(-c(cycle_mean, cycle_sd, quantity_mean, quantity_sd, comments, high_sd, no_amp, outlier_rg, exp_fail))

# make an inoculation treatment column
unique(qdat$sample)
# use first three letters if it's a standard or control
# use last two letters if it's a root sample
# use third character if it's a sample
qdat <- qdat %>%
  mutate(inoc = case_when(
    substr(sample, 1, 3) %in% 
      c("NTC", "NAC", "PAV", "RPV") ~ substr(sample, 1, 3),
    TRUE ~ substr(sample, 3, 3) 
  ))
unique(qdat$inoc)

# re-assign some levels
qdat$inoc <- recode(qdat$inoc, a ="blank", C ="coinfection", H="healthy", P ="PAV", R ="RPV")

# re-assign tasks for controls, standards, and samples.
unique(qdat$task)
qdat <- qdat %>%
  mutate(task = case_when(
    inoc %in% c("blank", "NTC", "NAC") ~ "control",
    (substr(sample, 1, 3) == "PAV" & target=="PAV") | 
      (substr(sample, 1, 3) =="RPV" & target=="RPV") ~ "standard",
    (substr(sample, 1, 3) == "PAV" & target=="RPV") | 
      (substr(sample, 1, 3) =="RPV" & target=="PAV") ~ "nonTargetStandard",
    TRUE ~ "sample"
  ))

# indicate source material
qdat <- qdat %>%
  mutate(material = case_when(
    task == "sample" & substr(sample, 7, 8) != "Rt" & substr(sample,8,9) != "Rt" ~ "shoot",
    task %in% c("control", "standard", "nonTargetStandard") ~ "RNA",
    substr(sample, 7, 8) == "Rt"| substr(sample,8,9) == "Rt" ~ "root"
  ))

# prepared standards assuming 2.5 ul would be added to wells (conc = 1eX/2.5 ul), but started adding 7 ul instead of 2.5 ul with group 3: replace quantity of standards from all groups except 1 and 2 with quantities divided by 2.5 and multiplied by 7 
qdat <- qdat %>%
  mutate(quantity = ifelse(task == "standard" & !(q_group %in% c("01","02")), quantity / 2.5 * 7, quantity))

# create a nutrient treatment column
qdat <- qdat %>%
  mutate(nutrient = case_when(
    task == "sample" & is.na(as.numeric(substr(sample, 5, 5))) == T ~ paste(substr(sample, 4, 4), substr(sample, 5, 5), sep = ""), 
    task == "sample" & !is.na(as.numeric(substr(sample, 5, 5))) == T ~ substr(sample, 4, 4)))
unique(qdat$nutrient)

# alter some nutrient values
filter(qdat, nutrient == "N/") %>% print.data.frame()
qdat$nutrient <- recode(qdat$nutrient, C = "low", L = "low", NP = "N+P", "N/" = "N/N+P")

# create columns for time, round, and replicate
qdat <- qdat %>%
  mutate(round = case_when(task == "sample" ~ substr(sample, 1, 1)),
         time = case_when(task == "sample" ~ substr(sample, 2, 2)),
         replicate = case_when(task == "sample" & !is.na(as.numeric(substr(sample, nchar(sample), nchar(sample)))) ~ as.numeric(substr(sample, nchar(sample), nchar(sample))),
                               task == "sample" & is.na(as.numeric(substr(sample, nchar(sample), nchar(sample)))) ~ as.numeric(substr(sample, 1, 1)) + 1))

qdat  %>%
  filter(task == "sample" & nutrient == "N+P") %>%
  select(sample, round, time, replicate)


#### edit individual samples ####

# incorrectly entered sample name
qdat$sample <- recode(qdat$sample,"RPV 1E5" = "RPV1E5")

# incorrectly entered target
unique(qdat$target)
filter(qdat, target == "") %>% print.data.frame()
qdat <- qdat %>%
  mutate(target = ifelse(target == "", "RPV", target))

# incorrectly entered sample name
# 16CN1 and 16CN2 were switched starting with extraction. Found out when  I did second extraction of 16CN2, but I kept the incorrect label to be consistent.
subset(qdat,sample %in% c("16CN1", "16CN2"))

qdat <- qdat %>%
  mutate(sample = case_when(
    sample == "16CN1" ~ "16CN2",
    sample == "16CN2" ~ "16CN1",
    TRUE ~ sample),
    replicate = case_when(
      sample == "16CN1" ~ 1,
      sample == "16CN2" ~ 2,
      TRUE ~ replicate
    )
  )

# look at notes
unique(qdat$qPCR_notes)

# remove column: low volume, contaminated
qdat %>%
  filter(qPCR_notes == "low volume") # remove the samples, leave the controls
qdat %>%
  filter(grepl("contaminated", qPCR_notes)) # just one sample
qdat <- qdat %>%
  mutate(remove = case_when(
    qPCR_notes == "low volume" & task == "sample" ~ 1,
    grepl("contaminated", qPCR_notes) ~ 1,
    TRUE ~ 0
    ))

# remove 12CN/NP2 and 12CN/NP2b - one is 12CN2 and one is 12CNP2, but don't know which
qdat <- qdat %>%
  mutate(remove = case_when(
    sample %in% c("12CN/NP2", "12CN/NP2b") ~ 1,
    TRUE ~ remove
  ))


#### error check ####

colnames(qdat)

# make sure samples have all their info
filter(qdat, task == "sample") %>%
  select(-c(quantity, qPCR_notes, cycle)) %>%
  filter(!complete.cases(.))

# summarize sample number
qdat %>%
  filter(task == "sample") %>%
  group_by(round, time, inoc, nutrient, target) %>%
  summarize(reps = length(unique(replicate))) %>%
  ggplot(aes(x = reps)) +
  geom_histogram(binwidth = 1)

# examine non-samples
qdat %>%
  filter(task != "sample") %>%
  group_by(task, target, quantity) %>%
  summarize(reps = n()) %>%
  print.data.frame()

# examine standards
qdat %>%
  filter(task == "standard") %>%
  ggplot(aes(x = log10(quantity))) +
  geom_histogram() +
  facet_wrap(~target)

# look at high values
qdat %>%
  filter(task == "standard" & log10(quantity) > 7.8) # 1e8 and 1e9 standards in group 1

  
#### standard curves ####

# PAV Function
stdPAVfun<-function(dat){
  
  # extract min, max, and curve characteristics
  stdPAVdf <- dat %>%
    subset(!is.na(cycle) & inoc == "PAV" & task == "standard") %>%
    group_by(q_group) %>%
    summarise(PAVmin = min(quantity), 
          PAVminRep = sum(quantity == min(quantity)),
          PAVmax = max(quantity), 
          PAVslope = lm(cycle ~ log10(quantity))$coefficients[2], 
          PAVint=lm(cycle ~ log10(quantity))$coefficients[1])
  
  # calculate efficiency
  stdPAVdf$PAVefficiency <- 100 * (10^(1/abs(stdPAVdf$PAVslope)) - 1)
  
  # plot - can add to return as part of list if needed
  #stdPAVplot <- dat %>%
  #  subset(!is.na(cycle) & inoc == "PAV" & task == "standard") %>%
  #  ggplot(aes(x = log10(quantity), y = cycle)) + 
  #  geom_point(size=2) + 
  #  xlab("quantity (log)") + 
  #  ylab("Detection Cycle") + 
  #  geom_smooth(method=lm,se=F)
  
  return(stdPAVdf)
}

# RPV Function
stdRPVfun<-function(dat){
  
  # extract min, max, and curve characteristics
  stdRPVdf <- dat %>%
    subset(!is.na(cycle) & inoc == "RPV" & task == "standard") %>%
    group_by(q_group) %>%
    summarise(RPVmin = min(quantity), 
          RPVminRep = sum(quantity == min(quantity)),
          RPVmax = max(quantity), 
          RPVslope = lm(cycle ~ log10(quantity))$coefficients[2], 
          RPVint=lm(cycle ~ log10(quantity))$coefficients[1])
  
  # calculate efficiency
  stdRPVdf$RPVefficiency <- 100 * (10^(1/abs(stdRPVdf$RPVslope)) - 1)
  
  # plot - can add to return as part of list if needed
  #stdRPVplot <- dat %>%
  #  subset(!is.na(cycle) & inoc == "RPV" & task == "standard") %>%
  #  ggplot(aes(x = log10(quantity), y = cycle)) + 
  #  geom_point(size=2) + 
  #  xlab("quantity (log)") + 
  #  ylab("Detection Cycle") + 
  #  geom_smooth(method=lm,se=F)
  
  return(stdRPVdf)
}

# contamination function
confun <- function(dat){
  
  fulldf <- dat %>%
    filter(task == "control"|task == "nonTargetStandard") %>%
    group_by(q_group, target) %>%
    summarise(contamination = min(cycle, na.rm = T)) %>%
    spread(target, contamination) %>%
    rename(PAVcont = PAV, RPVcont = RPV) %>%
    merge(dat, ., all = T)
  
  stddf <- fulldf %>%
    group_by(q_group) %>%
    filter(task == "standard") %>%
    mutate(stdRem = case_when(
      target == "PAV" & cycle > unique(PAVcont) ~ 1,
      target == "RPV" & cycle > unique(RPVcont) ~ 1)) %>%
    group_by(q_group, target) %>%
    summarise(removed = sum(stdRem, na.rm = T)) %>%
    spread(target, removed) %>%
    rename(PAVstdRem = PAV, RPVstdRem = RPV) %>%
    merge(fulldf, ., all = T)
  
  PAVdf <- fulldf %>%
    filter(cycle < PAVcont) %>%
    stdPAVfun
  
  RPVdf <- fulldf %>%
    filter(cycle < RPVcont) %>%
    stdRPVfun
  
  findf <- stddf %>%
    merge(PAVdf, all = T) %>%
    merge(RPVdf, all = T)
  
  return(findf)
}

# some standards appear to be prepared incorrectly - previously went through each run individually (see qPCR-group-comments.csv for notes)

# check notes
qdat %>%
  filter(q_group%in%c("04", "05", "07", "09") & target == "RPV") %>%
  select(qPCR_notes) %>%
  unique()

# add note about standards
qdat <- qdat %>%
  mutate(
    qPCR_notes = case_when(
      q_group%in%c("04", "05", "07", "09") & target == "RPV" & is.na(qPCR_notes) ~ "1E5 standards may have been prepared incorrectly",
      q_group%in%c("04", "05", "07", "09") & target == "RPV" & !is.na(qPCR_notes) ~ paste(qPCR_notes, "1E5 standards may have been prepared incorrectly", sep = ", "),
      TRUE ~ qPCR_notes)
  )

# analyze standard curves
qdat2 <- qdat %>%
  confun()

# check that confun worked
head(qdat2)

qdat %>%
  filter(q_group == "01") %>%
  filter(task == "control"|task == "nonTargetStandard") %>%
  group_by(target) %>%
  summarise(contamination = min(cycle, na.rm = T)) %>%
  spread(target, contamination) %>%
  rename(PAVcont = PAV, RPVcont = RPV) 

unique(qdat2$PAVcont)

qdat2 %>%
  filter(is.na(PAVcont)) %>%
  select(target) %>%
  unique()

unique(qdat2$RPVcont)

qdat2 %>%
  filter(is.na(RPVcont)) %>%
  select(target) %>%
  unique()

# modify samples with new standard curve
qdat2 <- qdat2 %>%
  mutate(quant_adj = case_when(
    task != "standard" & target == "PAV" ~ 10 ^ ((cycle - PAVint) / PAVslope),
    task != "standard" & target == "RPV" ~ 10 ^ ((cycle - RPVint) / RPVslope),
    task == "standard" ~ quantity))

# check that new standard curve conversion worked
qdat2 %>%
  ggplot(aes(x = log10(quantity), y = log10(quant_adj))) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0)

# 4 obvious outliers
qdat2 %>%
  filter(log10(quantity) < 0 & log10(quant_adj) > 0.5) # these are too low to be reliable

qdat2 %>%
  filter(log10(quantity) < 1 & log10(quant_adj) > 2) # from group 35, which had a high conncentration contamination, so lower concentration samples will be left out


#### create sample dataset ####

samp <- qdat2 %>%
  filter(task == "sample") %>%
  mutate(
    sample = paste(round, time, inoc, nutrient, replicate, material, sep = "")
  ) %>%
  as_tibble()

# identify which samples will be removed:
  # standard curve efficiency outside of boundaries
  # quantity greater than 1e3, but less than standard curve minimum
sum(samp$remove) # 15

samp <- samp %>%
  mutate(
    remove = case_when(
      target == "RPV" & (round(RPVefficiency) < 85 | round(RPVefficiency) > 115) ~ 1,
      target == "PAV" & (round(PAVefficiency) < 85 | round(PAVefficiency) > 115) ~ 1,
      target == "RPV" & (quant_adj >= 1e3 & quant_adj < RPVmin) ~ 1,
      target == "PAV" & (quant_adj >= 1e3 & quant_adj < PAVmin) ~ 1,
      TRUE ~ remove
    )
  )

sum(samp$remove) # 672

# identify duplicate samples and select ones with detection or lower variance
samp_d <- samp %>%
  filter(remove == 0) %>%
  select(q_group, target, sample) %>%
  distinct() %>%
  mutate(
    dup = duplicated(select(., target, sample))
    ) %>%
  filter(dup == T) %>%
  select(-c(q_group, dup)) %>%
  semi_join(filter(samp, remove == 0), .) %>%
  mutate(quant_adj = case_when(is.na(quant_adj) ~ 0,
                               TRUE ~ quant_adj),
         detect = case_when(target == "RPV" & quant_adj >= RPVmin ~ 1,
                            target == "PAV" & quant_adj >= PAVmin ~ 1,
                            TRUE ~ 0)) %>%
  group_by(q_group, sample, target) %>%
  summarise(quant_var = sd(quant_adj),
            quant_mean = mean(quant_adj),
            detect = sum(detect) > 0) %>%
  group_by(sample, target) %>%
  mutate(
    min_var_group = case_when(quant_var == min(quant_var, na.rm = T) ~ 1,
                              TRUE ~ 0),
    n_q_group = length(unique(q_group)),
    n_min = sum(min_var_group),
    n_det = sum(detect),
    sample_q_tar = paste(sample, q_group, target, sep = "_")) %>%
  ungroup()

# see if variance = 0, but mean doesn't
samp_d %>%
  filter(quant_var == 0 & quant_mean != 0) # no

# see how these criteria are related
samp_d %>%
  ggplot(aes(x = sample, y = log(quant_var+1), colour = detect)) +
  geom_point()
# all detected samples have higher variance than undetected samples

# look at groups with NA variance
samp_d %>% 
  filter(is.na(quant_var)) %>%
  select(min_var_group) %>%
  unique() # none are counted as the min var group

# indicate which samples should be removed
# not detected when another group was
# higher than minimum variance when multiple groups were detected
samp_d_r <- samp_d %>%
  mutate(
    remove = case_when(detect == F & n_det > 0 ~ 1,
                       detect == T & n_det > 1 & min_var_group == 0 ~ 1,
                       TRUE ~ 0
    )
  ) %>%
  filter(remove == 1)

# transfer this over to main dataset
samp <- samp %>%
  mutate(
    sample_q_tar = paste(sample, q_group, target, sep = "_"),
    remove = case_when(sample_q_tar %in% samp_d_r$sample_q_tar ~ 1,
                       TRUE ~ remove)) %>%
  select(-sample_q_tar)
sum(samp$remove) # 886

# merge with experiment and molecular information
dat <- samp %>%
  mutate(round = as.numeric(round),
         time = as.numeric(time)) %>%
  left_join(emdat)
