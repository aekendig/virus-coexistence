## Goal: evaluate treatment effects on virus concentration in Experiment 1


#### set up ####

# source data
source("./code/qPCR_raw_data_processing.R") # clears environment and loads tidyverse

# clear all except dataset
rm(list = setdiff(ls(), c("dat")))

# load libraries
library(brms) # version used: 2.7.0
library(tidybayes) # version used: 1.0.4
library(bayesplot) # version used: 1.6.0

# import data
sdat <- read_csv("./data/sample_exp_molc_data.csv")

# load models for priors
# load("./output/lacroix_output/lacroix_concentration_pav.rda")
# load("./output/lacroix_output/lacroix_concentration_rpv.rda")


#### edit data ####

# days post inoculation
dpi <- tibble(
  time = 1:8,
  dpi = c(5, 8, 12, 16, 19, 22, 26, 29)
)

# sample size in experiments (part of Table S2)
dat %>%
  filter(material == "shoot" & inoc != "healthy") %>%
  select(target, nutrient, inoc, round, time, replicate) %>%
  unique() %>%
  group_by(target, nutrient, inoc) %>%
  summarise(reps = n()) %>%
  data.frame()

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
                               TRUE ~ quant_adj),
         quant_zero = case_when(quant_adj == 0 ~ 1,
                                TRUE ~ 0)) %>%
  full_join(dpi)


#### examine distribution and edit data ####

# histogram of values
dat2 %>%
  ggplot() +
  geom_histogram(aes(x = quant_adj, fill = quant_zero)) 

# one value is very high
filter(dat2, quant_adj > 1e10) %>% data.frame() # a PAV value much higher than max

# look at all values above max
filter(dat2, (target == "RPV" & quant_adj > RPVmax) | (target == "PAV" & quant_adj > PAVmax)) # only 3, remove these too

# remove values above standard curve
dat3 <- dat2 %>%
  filter((target == "RPV" & quant_adj <= RPVmax) | (target == "PAV" & quant_adj <= PAVmax))

# re-do histogram
dat3 %>%
  ggplot() +
  geom_histogram(aes(x = quant_adj, fill = quant_zero))

# log-tranformed values
dat3 %>%
  filter(quant_zero == 0) %>%
  ggplot() +
  geom_histogram(aes(x = log(quant_adj)))


#### visualize sources of variation ####

# look at qPCR wells
dat3 %>%
  ggplot(aes(x = well, y = cycle, colour = target)) +
  geom_point(alpha = 0.5) +
  stat_summary(fun.data = "mean_cl_boot")
# not a strong trend due to wells and the ones that stick out have fewer replicates

# look at PAV and RPV concentrations in same well
wcdat <- dat3 %>%
  select(q_group, well, target, cycle) %>%
  spread(target, cycle) %>%
  group_by(well) %>%
  summarise(mPAV = mean(PAV, na.rm = T), seP = sd(PAV, na.rm = T)/sqrt(length(!is.na(PAV))), mRPV = mean(RPV, na.rm = T), seR = sd(RPV, na.rm = T)/sqrt(length(!is.na(RPV)))) %>%
  ungroup()

ggplot(wcdat, aes(x = mPAV, y = mRPV)) +
  geom_point() +
  geom_errorbar(aes(ymin = mRPV - seR, ymax = mRPV + seR)) +
  geom_errorbarh(aes(xmin = mPAV - seP, xmax = mPAV + seP)) +
  geom_smooth(method = "lm")
# wells that produce higher PAV values don't necessarily produce high RPV values and vice versa

cor.test(wcdat$mRPV, wcdat$mPAV)
# not a strong correlation - well does not seem to affect concentration

# correlation between well's distance from mean and number of reps
dat3 %>%
  mutate(mean_cycle = mean(cycle, na.rm = T)) %>%
  group_by(well, target, mean_cycle) %>%
  summarise(well_mean = mean(cycle, na.rm = T),
            well_reps = length(well)) %>%
  mutate(well_dist = well_mean - mean_cycle) %>%
  ggplot(aes(x = well_reps, y = well_dist)) +
  geom_point() 
# all are within about 5 cycles of the mean and those outside have very low replication

# look at standard deviation among technical replicates
dat3 %>%
  filter(inoc != "healthy") %>%
  group_by(q_group, target, dpi, inoc, nutrient, sample) %>%
  summarise(tech_sd = sd(cycle, na.rm = T),
            zero = ifelse(sum(quant_zero) > 0, 1, 0)) %>%
  ggplot(aes(x = dpi, y = tech_sd, colour = inoc, shape = as.factor(zero))) + 
  geom_point() +
  facet_grid(target ~ nutrient, scales = "free")
# all without zero's are less than 3 and there isn't a strong relationship to experimental treatment

# look at mean and standard deviation of technical replicates without zero's
dat3 %>%
  filter(quant_zero == 0 & inoc != "healthy") %>%
  group_by(q_group, target, dpi, inoc, nutrient, sample) %>%
  summarise(tech_sd = sd(cycle),
            tech_mean = mean(cycle)) %>%
  ggplot(aes(x = tech_mean, y = tech_sd, colour = target)) + 
  geom_point()
# higher sd samples generally have higher mean values

# look at experimental rounds
dat3 %>%
  ggplot(aes(x = round, y = cycle, colour = as.factor(replicate))) +
  geom_point(alpha = 0.5) +
  stat_summary(fun.data = "mean_cl_boot", size = 1.5, position = position_dodge(0.3)) +
    facet_wrap(~target)
# rounds are pretty similar, replicates within rounds are not necessarily more similar than across rounds (except 6 and 7 for RPV)

# look at rounds by treatment
dat3 %>%
  ggplot(aes(x = round, y = cycle, colour = as.factor(replicate))) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(0.3), alpha = 0.7) +
  geom_point(alpha = 0.3) +
  facet_grid(nutrient~target)

# look at round by day
dat3 %>%
  ggplot(aes(x = round, y = cycle, colour = as.factor(replicate))) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(0.3), alpha = 0.7) +
  geom_point(alpha = 0.3) +
  facet_grid(target~dpi)
# not all days are in all groups, may be causing group differences


#### average technical replicates ####

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
dups <-dat4 %>%
  group_by(target, sample) %>%
  mutate(dup = duplicated(sample)) %>%
  filter(dup == T) %>%
  select(sample, target) %>%
  ungroup()

dups %>%
  left_join(dat4) %>%
  select(sample, target, quant_ul) %>%
  data.frame() 
# all are zero's - combine samples if needed


#### visualize treatment effects ####

# look at mean values
dat4 %>%
  ggplot(aes(x = dpi, y = quant_ul)) + 
  stat_summary(data = filter(dat4, quant_zero == 0), fun.data = "mean_se") +
  stat_summary(data = filter(dat4, quant_zero == 0), fun.data = "mean_se", geom = "line") +
  geom_point(data = filter(dat4, quant_zero == 1), color = "blue", alpha = 0.5) + 
  facet_grid(target ~ inoc, scales = "free") # PAV growth is delayed by coinfection, RPV is enhanced (but more variable)

# look at mean values by nutrient - PAV
dat4 %>%
  ggplot(aes(x = dpi, y = quant_ul)) + 
  stat_summary(data = filter(dat4, quant_zero == 0 & target == "PAV"), fun.data = "mean_se") +
  stat_summary(data = filter(dat4, quant_zero == 0 & target == "PAV"), fun.data = "mean_se", geom = "line") +
  geom_point(data = filter(dat4, quant_zero == 1 & target == "PAV"), color = "blue", alpha = 0.5) + 
  facet_grid(inoc ~ nutrient, scales = "free") # the peak in coinfection is driven by low nutrients, but it is driven by N when PAV is alone

# look at PAV in coinfection
dat4 %>%
  ggplot(aes(x = dpi, y = quant_ul)) + 
  stat_summary(data = filter(dat4, quant_zero == 0 & target == "PAV" & inoc == "coinfection"), fun.data = "mean_se") +
  stat_summary(data = filter(dat4, quant_zero == 0 & target == "PAV" & inoc == "coinfection"), fun.data = "mean_se", geom = "line") +
  geom_point(data = filter(dat4, quant_zero == 1 & target == "PAV" & inoc == "coinfection"), color = "blue", alpha = 0.5) + 
  facet_wrap(~nutrient, nrow = 2, scales = "free") # peaks in the middle with high nutrients (like in single infection), especially P, peaks later with lower

# look at mean values by nutrient - log-transformed PAV
dat4 %>%
  filter(quant_zero == 0 & target == "PAV" & !(inoc %in% c("healthy", "RPV"))) %>%
  ggplot(aes(x = dpi, y = log(quant_ul))) + 
  stat_summary(fun.data = "mean_se") +
  stat_summary(fun.data = "mean_se", geom = "line") +
  facet_grid(inoc ~ nutrient, scales = "free")

# look at mean values by nutrient - RPV
dat4 %>%
  ggplot(aes(x = dpi, y = quant_ul)) + 
  stat_summary(data = filter(dat4, quant_zero == 0 & target == "RPV"), fun.data = "mean_se") +
  stat_summary(data = filter(dat4, quant_zero == 0 & target == "RPV"), fun.data = "mean_se", geom = "line") +
  geom_point(data = filter(dat4, quant_zero == 1 & target == "RPV"), color = "blue", alpha = 0.5) + 
  facet_grid(inoc ~ nutrient, scales = "free") # it looks like there are multiple peaks in the temporal dynamics and they occur at different times depending on the nutrient and inoculation treatment, highest peaks with P addition

# look at mean values by nutrient - log-transformed RPV
dat4 %>%
  filter(quant_zero == 0 & target == "RPV" & !(inoc %in% c("healthy", "PAV"))) %>%
  ggplot(aes(x = dpi, y = log(quant_ul))) + 
  stat_summary(fun.data = "mean_se") +
  stat_summary(fun.data = "mean_se", geom = "line") +
  facet_grid(inoc ~ nutrient, scales = "free")


#### format data for models ####

# edit data
# remove values that are too low to use
d.at <- dat4 %>%
  filter(quant_zero == 0 &
           inoc %in% c("PAV", "coinfection", "RPV")) %>%
  mutate(co = ifelse(inoc == "coinfection", 1, 0),
         log_quant = log(quant_ul),
         conc = quant_ul * 50 / mass_ext_mg,
         log_conc = log(conc),
         quant_rd = round(quant_ul))

# concentration values
d.at %>%
  ggplot(aes(x = quant_rd)) +
  geom_histogram() + 
  facet_wrap(~target, scales = "free")

d.at %>%
  ggplot(aes(x = log_conc)) +
  geom_histogram() + 
  facet_wrap(~target, scales = "free")

# replicates within same round and qPCR group
d.at %>%
  group_by(target, round, dpi, inoc, high_N, high_P, q_group) %>%
  summarise(reps = length(unique(sample))) %>%
  filter(reps >1) # 19

# select rows from sdat
sdat2 <- sdat %>%
  filter(material == "shoot") %>%
  select(round, time, inoc, nutrient, replicate, RTPCR_PAV, RTPCR_RPV)

# accidental RPV inoculations
(acc.r <- d.at %>%
  left_join(sdat2) %>%
  filter(inoc == "PAV" & (target == "RPV" | RTPCR_RPV == 1)) %>%
  select(sample, round, time, inoc, nutrient, replicate, quant_zero, conc, RTPCR_RPV))

d.at %>%
  filter(sample %in% acc.r$sample) %>%
  select(sample, nutrient, target, inoc, quant_zero, conc)
# 5 samples will be removed
# one has no PAV value anyway (probably too low)

# accidental PAV inoculations  
(acc.p <- d.at %>%
  left_join(sdat2) %>%
  filter(inoc == "RPV" & (target == "PAV" | RTPCR_PAV == 1)) %>%
  select(sample, round, time, inoc, nutrient, replicate, quant_zero, conc, RTPCR_PAV)) 

d.at %>%
  filter(sample %in% acc.p$sample) %>%
  select(sample, nutrient, target, inoc, quant_zero, conc)
# 7 samples will be removed

# accidental infections from all samples (Table S2)
dat2 %>%
  group_by(target, dpi, time, inoc, nutrient, round, replicate, sample, PAVmin, PAVint, PAVslope, RPVmin, RPVint, RPVslope, q_group, RTPCR_PAV, RTPCR_RPV) %>%
  summarise(tech_cycle = mean(cycle, na.rm = T)) %>%
  mutate(quant = case_when(target == "PAV" ~ 10 ^ ((tech_cycle - PAVint) / PAVslope),
                           target == "RPV" ~ 10 ^ ((tech_cycle - RPVint) / RPVslope))) %>%
  ungroup() %>%
  filter((inoc == "PAV" & ((target == "RPV" & quant > RPVmin) | (RTPCR_RPV == 1 & !is.na(quant))))|
           (inoc == "RPV" & ((target == "PAV" & quant > PAVmin) | (RTPCR_PAV == 1 & !is.na(quant))))) %>%
  select(inoc, sample, nutrient, quant, RTPCR_RPV, RTPCR_PAV) %>%
  arrange(inoc, nutrient) %>%
  data.frame()
# not going to exclude samples that were positive for RT-PCR, but had no detectable density (NA)

# data by virus
d.at.p <- d.at %>%
  filter(target == "PAV" & inoc != "RPV" & !(sample %in% acc.r$sample)) 

d.at.r <- d.at %>%
  filter(target == "RPV" & inoc != "PAV" & !(sample %in% acc.p$sample))

# new sample sizes (Table S2)
d.at.p %>%
  group_by(inoc, nutrient) %>%
  summarise(reps = n())

d.at.r %>%
  group_by(inoc, nutrient) %>%
  summarise(reps = n())


#### log-transformed models, uninformative priors ####

# PAV model
m.lu.p <- brm(data = d.at.p, family = gaussian,
              log_conc ~ co * high_N * high_P,
              autocor = cor_ar(~time),
              prior <- c(prior(normal(0, 10), class = Intercept),
                         prior(normal(0, 1), class = b)),
              iter = 6000, warmup = 1000, chains = 3, cores = 2)

summary(m.lu.p) # coinfection increases concentration
plot(m.lu.p) # convergence among chains
# plot(marginal_effects(m.lu.p), points = T)
pp_check(m.lu.p, nsamples = 100) # model distributions slightly higher

# save model
save(m.lu.p, file = "./output/concentration_analysis_uninformative_pav.rda")

# RPV model 
m.lu.r <- update(m.lu.p, newdata = d.at.r)

# inspect model
summary(m.lu.r) # no strong effects
plot(m.lu.r) # convergence among chains
# plot(marginal_effects(m.lu.r), points = T)
pp_check(m.lu.r, nsamples = 100) # pretty close

# save model
save(m.lu.r, file = "./output/concentration_analysis_uninformative_rpv.rda")


#### count models, uninformative priors ####

# check mean and variance
mean(d.at.p$quant_rd)
var(d.at.p$quant_rd) # much larger variance
mean(d.at.r$quant_rd)
var(d.at.r$quant_rd) # much larger variance

# PAV model (you can't currently implement autoregressive models, time included as a random effect)
# m.cu.p <- brm(data = d.at.p, family = negbinomial(),
#               quant_rd ~ offset(mass_ext_mg) + co * high_N * high_P + (1|time),
#               prior <- c(prior(normal(460, 100), class = Intercept),
#                          prior(normal(0, 1), class = b)),
#               iter = 6000, warmup = 1000, chains = 1, cores = 1,
#               control = list(adapt_delta = 0.99))
# 
# # inspect model
# summary(m.cu.p)
# yrep.cu.p <- posterior_predict(m.cu.p, draws = 100)
# ppc_dens_overlay(log(d.at.p$quant_rd), log(yrep.cu.p)) # not a very good fit
# 
# # RPV model
# m.cu.r <- update(m.cu.p, newdata = d.at.r)
# 
# # inspect model
# summary(m.cu.r) # positive effect of P
# yrep.cu.r <- posterior_predict(m.cu.r, draws = 100)
# ppc_dens_overlay(log(d.at.r$quant_rd), log(yrep.cu.r)) # not a very good fit


#### log-transformed models, specific priors ####

# PAV model
# summary(m.p)
m.li.p <- brm(data = d.at.p, family = gaussian,
              log_conc ~ co * high_N * high_P,
              autocor = cor_ar(~time),
              prior = c(prior(normal(0.25, 0.73), class = b, coef = co),
                        prior(normal(0.08, 0.76), class = b, coef = high_N),
                        prior(normal(-0.19, 0.83), class = b, coef = high_P),
                        prior(normal(0, 1), class = b, coef = co:high_N),
                        prior(normal(0, 1), class = b, coef = co:high_P),
                        prior(normal(0.17, 1.06), class = b, coef = high_N:high_P),
                        prior(normal(0, 1), class = b, coef = co:high_N:high_P),
                        prior(normal(0, 10), class = Intercept)),
              iter = 6000, warmup = 1000, chains = 3, cores = 2)

# save model
save(m.li.p, file = "./output/concentration_analysis_informative_pav.rda")

# inspect model
summary(m.li.p) # didn't change estimates too much
plot(m.li.p) # convergence among chains
# plot(marginal_effects(m.li.p), points = T)
pp_check(m.li.p, nsamples = 100) # model distributions slightly higher, similar to model with uninformative priors

# RPV model
# summary(m.r)
m.li.r <- update(m.li.p,
                 newdata = d.at.r,
                  prior = c(prior(normal(-0.35, 0.66), class = b, coef = co),
                            prior(normal(-0.21, 0.58), class = b, coef = high_N),
                            prior(normal(0.79, 0.97), class = b, coef = high_P),
                            prior(normal(-0.54, 0.84), class = b, coef = co:high_N),
                            prior(normal(-1.22, 1.39), class = b, coef = co:high_P),
                            prior(normal(-0.26, 1.15), class = b, coef = high_N:high_P),
                            prior(normal(0.74, 1.61), class = b, coef = co:high_N:high_P),
                            prior(normal(0, 10), class = Intercept)))

# save model
save(m.li.r, file = "./output/concentration_analysis_informative_rpv.rda")

# inspect model
summary(m.li.r) # similar estimates to uninformative priors, tighter confidence intervals
plot(m.li.r) # convergence among chains
# plot(marginal_effects(m.li.r), points = T)
pp_check(m.li.r, nsamples = 100) # similar to uninformative priors


# #### log-transformed models, uninformative priors, late dpi ####
# 
# # subset data
# d.at.p.sub <- d.at.p %>% filter(dpi > 10)
# d.at.r.sub <- d.at.r %>% filter(dpi > 10)
# 
# # PAV model
# m.lud.p <- update(m.lu.p, newdata = d.at.p.sub)
# 
# # save model
# save(m.lud.p, file = "./output/concentration_analysis_uninformative_late_pav.rda")
# 
# # inspect model
# summary(m.lud.p) # similar estimates to full dataset
# plot(m.lud.p)
# pp_check(m.lud.p, nsamples = 100) # estimates a little high
# 
# # RPV model
# m.lud.r <- update(m.lu.r, newdata = d.at.r.sub)
# 
# # save model
# save(m.lud.r, file = "./output/concentration_analysis_uninformative_late_rpv.rda")
# 
# # inspect model
# summary(m.lud.r)  # similar estimates to full dataset
# plot(m.lud.r)
# pp_check(m.lud.r, nsamples = 100) 
# 
# 
# #### log-transformed models, specific priors, late dpi ####
# 
# # PAV model
# m.lid.p <- update(m.li.p, newdata = d.at.p.sub)
# 
# # save model
# save(m.lid.p, file = "./output/concentration_analysis_informative_late_pav.rda")
# 
# # inspect model
# summary(m.lid.p) # similar estimates to uninformed
# plot(m.lid.p)
# pp_check(m.lid.p, nsamples = 100) # same issue as uninformed
# 
# # RPV model
# m.lid.r <- update(m.li.r, newdata = d.at.r.sub)
# 
# # save model
# save(m.lid.r, file = "./output/concentration_analysis_informative_late_rpv.rda")
# 
# # inspect model
# summary(m.lid.r) # similar estimates to uninformed
# plot(m.lid.r)
# pp_check(m.lid.r, nsamples = 100) 


#### save data for plotting ####

# save file
write_csv(d.at.p, "./output/concentration_analysis_pav_data.csv")
write_csv(d.at.r, "./output/concentration_analysis_rpv_data.csv")
write_csv(dat4, "./output/concentration_analysis_all_data.csv")
