
### Process qPCR Data from Ex2

# Updated from 7/10 because some of the mass data was missing, and I updated it
# Updated from 7/16 because I want the contamination values in the dataset
# Updated from 7/17 because I wanted samples with all reps below the standard curve to count as "not quantifiable" in Quant 2 and not be removed from the analysis
# Updated from 9/24 because I wanted to deal with contamination by cutting off the standard curve below it
# Updated from 10/12 because there were NA's for quant2 that shouldn't have been. NA's in low Samps for RPV group 22 and 0's for low samps when undetected
# Updated from 10/16 to double check code and make small tweaks
# Updated from 4/29 to remove duplicates (PAV that were tested in group 22 and other groups - removed all from 22)

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# import data
qdat <- read_csv("./intermediate-data/qPCR_data_compiled.csv")
base <- read_csv("./data/plant_data_071617.csv")


#### edit data ####

# remove summarizing columns and empty wells
qdat2 <- qdat %>%
  select(-c(cycle_mean, cycle_sd, quantity_mean, quantity_sd, comments, high_sd, no_amp, outlier_rg, exp_fail)) %>%
  filter(!is.na(sample))

# examine sample names, tasks
unique(qdat2$sample)
unique(qdat2$task)

# make cycle numeric
# update tasks
# fix label errors (manually checked the group-specific one)
# assigned quantities for standards
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

# number of tube labels
length(unique(qdat3$tube_label))
length(unique(base$tube_label))

# check for NA's
sum(is.na(qdat3$tube_label))
sum(is.na(base$tube_label))

# examine base NA's
filter(base, is.na(tube_label)) %>%
  select(expt_notes, extraction_notes, extraction_date)
# 18 samples not extracted, all have justifications

# remove samples with errors from base data
base2 <- filter(base, !is.na(tube_label))

# examine qdat NA's
filter(qdat3, is.na(tube_label)) %>%
  select(task) %>%
  unique()
# none are samples

# overlap in tube numbers
nrow(filter(qdat3, task == "sample" & !(tube_label %in% base2$tube_label)))
nrow(filter(base2, !(tube_label %in% qdat3$tube_label)))

# add base information to samples
# create combine set and replicate
qdat4 <- qdat3 %>%
  left_join(base2 %>%
              select(-sample)) %>%
  mutate(set_rep = paste(set, replicate, sep = "_"))


#### error check ####

# check for missing base data
qdat4 %>%
  filter(task == "sample" & is.na(nutrient)) %>% select(sample) %>% unique() # none

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
  mutate(dup =duplicated(tube_label)) %>%
  filter(dup == T) %>%
  data.frame() # all are in groups 21 and 22

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

# PAV Function
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

# RPV Function
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

# examine standard curves (uncomment within functions, comment return line)
# for(i in unique(qdat4$q_group)){
#   
#   stdPAVfun(filter(qdat4, q_group == i))
#   stdRPVfun(filter(qdat4, q_group == i))
#   
# }
# look pretty good - outliers aren't pulling the curve, most are linear

# contamination function
confun <- function(dat){
  
  fulldf <- dat %>%
    filter(task == "control"|task == "nonTargetStandard") %>%
    group_by(q_group, target) %>%
    summarise(contamination = min(cycle, na.rm = T)) %>%
    ungroup() %>%
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
    ungroup() %>%
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

# identify contamination as the highest concentrated control or non-target standard per target per run
# identify standards with lower concentration than these
# remove these standards and re-calculate standard curve
qdat5 <- qdat4 %>%
  confun()

# check that confun worked
head(qdat5)
qdat5 %>%
  group_by(PAVstdRem, RPVstdRem) %>%
  summarise(groups = length(unique(q_group)))

#### start here ####

# modify samples with new standard curve
qdat2 <- qdat2 %>%
  mutate(quant_adj = case_when(
    task != "standard" & target == "PAV" ~ 10 ^ ((cycle - PAVint) / PAVslope),
    task != "standard" & target == "RPV" ~ 10 ^ ((cycle - RPVint) / RPVslope),
    task == "standard" ~ quantity))








## Dealing with contamination

# Function for finding limit of standard curves, values below this limit, and contamination after this has been accounted for

contFun=function(dat){
	
	# Make CT numeric
	dat$CT=as.numeric(as.character(dat$CT))
	
	# Find the max CT and Quantity for each standard
	pmax=max(subset(dat,Task=="STANDARD"&substring(SampleName,1,3)=="PAV"&TargetName=="PAV"&!is.na(CT))$CT)
	rmax=max(subset(dat,Task=="STANDARD"&substring(SampleName,1,3)=="RPV"&TargetName=="RPV"&!is.na(CT))$CT)
	
	pmaxQ=max(subset(dat,Task=="STANDARD"&substring(SampleName,1,3)=="PAV"&TargetName=="PAV"&!is.na(CT))$Quantity)
	rmaxQ=max(subset(dat,Task=="STANDARD"&substring(SampleName,1,3)=="RPV"&TargetName=="RPV"&!is.na(CT))$Quantity)
	
	# Find the min CT and Quantity for each standard
	pmin=min(subset(dat,Task=="STANDARD"&substring(SampleName,1,3)=="PAV"&TargetName=="PAV"&!is.na(CT))$CT)
	rmin=min(subset(dat,Task=="STANDARD"&substring(SampleName,1,3)=="RPV"&TargetName=="RPV"&!is.na(CT))$CT)
	
	pminQ=min(subset(dat,Task=="STANDARD"&substring(SampleName,1,3)=="PAV"&TargetName=="PAV"&!is.na(CT))$Quantity)
	rminQ=min(subset(dat,Task=="STANDARD"&substring(SampleName,1,3)=="RPV"&TargetName=="RPV"&!is.na(CT))$Quantity)
	
	# Output dataframe 1
	out1=data.frame(pmax=pmax,rmax=rmax,pmin=pmin,rmin=rmin,pmaxQ=pmaxQ,rmaxQ=rmaxQ,pminQ=pminQ,rminQ=rminQ)
	
	# Check for contamination in standards or controls
	cont=subset(subset(subset(dat,substring(SampleName,1,3)=="NTC"|(substring(SampleName,1,3)=="PAV"&TargetName=="RPV")|(substring(SampleName,1,3)=="RPV"&TargetName=="PAV")),!is.na(CT)),select=c(SampleName,TargetName,CT,Quantity))
	
	# Make an empty df in case there are no contamination issues
	naDF=data.frame(SampleName=NA,TargetName=NA,CT=NA,Quantity=NA)
	
	# Merge with output df 1
	cont.out1=merge(cont,out1,all=T)
	cont.out2=merge(naDF,out1,all=T)
	
	# Mark contamination values within the standard curve as 1
	cont.out1$inSC<-with(cont.out1,ifelse(TargetName=="PAV"&CT<pmax,1,ifelse(TargetName=="RPV"&CT<rmax,1,0)))
	cont.out2$inSC=NA

	# Return
	ifelse(nrow(cont.out1)==0,return(cont.out2),return(cont.out1))
}

# Apply contFun to each dataset and store in a list
contOut=llply(qdatL,contFun)
contOut
# Groups with contamination inside the SC: 4, 6, 11, 21, 22 (all RPV except 4, which had some PAV)

unique(subset(contOut[[4]],select=c(TargetName,inSC))) # both
unique(subset(contOut[[6]],select=c(TargetName,inSC))) # RPV
unique(subset(contOut[[11]],select=c(TargetName,inSC))) # RPV
unique(subset(contOut[[21]],select=c(TargetName,inSC))) # RPV
unique(subset(contOut[[22]],select=c(TargetName,inSC)))# RPV


## Make dataset to analyze

# Add a qPCR Group ID and contamination information, make CT numeric
rows=c()
for(i in 1:length(ex2Files)){
	qdatL[[i]]$CT=as.numeric(as.character(qdatL[[i]]$CT))
	qdatL[[i]]$qPCRGroup=i
	qdatL[[i]]$pmax=contOut[[i]]$pmax[1]
	qdatL[[i]]$rmax=contOut[[i]]$rmax[1]
	qdatL[[i]]$pmin=contOut[[i]]$pmin[1]
	qdatL[[i]]$rmin=contOut[[i]]$rmin[1]
	qdatL[[i]]$pmaxQ=contOut[[i]]$pmaxQ[1]
	qdatL[[i]]$rmaxQ=contOut[[i]]$rmaxQ[1]
	qdatL[[i]]$pminQ=contOut[[i]]$pminQ[1]
	qdatL[[i]]$rminQ=contOut[[i]]$rminQ[1]
	qdatL[[i]]$pcon=min(subset(contOut[[i]],TargetName=="PAV")$CT)
	qdatL[[i]]$rcon=min(subset(contOut[[i]],TargetName=="RPV")$CT)
	rows=c(rows,nrow(qdatL[[i]]))
}
sum(rows)
# 4221 rows
head(qdatL[[1]]) 

# Make list into dataframe
dat=do.call(rbind.data.frame,qdatL)
head(dat)
nrow(dat)


## Assessing standard curves: Use dat (still has standards)

stdFun=function(df){
	
	# Dataframe for efficiency values
	stdEff=data.frame(group=unique(df$qPCRGroup),pEff=NA,rEff=NA,pInt=NA,pSlope=NA,rInt=NA,rSlope=NA)
	
	# Indexer (for when the groups aren't sequential, as used later)
	j=0
	

	# Cycle through groups and calculate curves and efficiency
	for(i in min(df$qPCRGroup):max(df$qPCRGroup)){
		
		j=j+1

		tempDat=subset(df,qPCRGroup==i)
	
		# Subset data for standards
		PAVdat=subset(tempDat,Task=="STANDARD"&substring(SampleName,1,3)=="PAV"&TargetName=="PAV")
		RPVdat=subset(tempDat,Task=="STANDARD"&substring(SampleName,1,3)=="RPV"&TargetName=="RPV")
	
		# Find slope and intercept of standard curve
		pInt=lm(CT~log10(Quantity),data=PAVdat)$coefficients[1]
		pSlope=lm(CT~log10(Quantity),data=PAVdat)$coefficients[2]
		rInt=ifelse(nrow(RPVdat)>0,lm(CT~log10(Quantity),data=RPVdat)$coefficients[1],NA)
		rSlope=ifelse(nrow(RPVdat)>0,lm(CT~log10(Quantity),data=RPVdat)$coefficients[2],NA)
	
		# Plot standard curves with data
	#print(ggplot(subset(tempDat,TargetName=="PAV"),aes(x=log10(Quantity),y=CT))+geom_point(size=2,aes(colour=Task))+xlab("Quantity (log 10)")+ylab("Detection Cycle")+geom_abline(intercept=pInt,slope=pSlope)+ggtitle(paste("PAV Plot. Group",tempDat$qPCRGroup[1],sep="")))
			#print(ggplot(subset(tempDat,TargetName=="RPV"),aes(x=log10(Quantity),y=CT))+geom_point(size=2,aes(colour=Task))+xlab("Quantity (log 10)")+ylab("Detection Cycle")+geom_abline(intercept=rInt,slope=rSlope)+ggtitle(paste("RPV Plot. Group",tempDat$qPCRGroup[1],sep="")))
		
		# Calculate efficiency
		pEff=100*(10^(1/abs(pSlope))-1)
		rEff=100*(10^(1/abs(rSlope))-1)
	
		# Add to dataframe
		stdEff$pEff[j]=pEff
		stdEff$rEff[j]=rEff
		stdEff$pInt[j]=pInt
		stdEff$rInt[j]=rInt
		stdEff$pSlope[j]=pSlope
		stdEff$rSlope[j]=rSlope
	}
	return(stdEff)
}

# Set up PDF with standard curves and data
#setwd("/Users/AmyKendig/Google Drive/Within Host Coexistence/qPCR/Ex2 Standards")
#pdf("Ex2_StdCurves_041017.pdf")
stdDat=stdFun(dat)
stdDat
#dev.off()


## Efficiencies outside of range 85%-115%
# Group 4: RPV eff=84
# Group 6: RPV eff=120
# Group 11: PAV eff=82
# Group 21: PAV eff=83
# Group 22: RPV eff=83
# Note: Groups with contamination: 4, 6, 11, 21, 22


## Check replication for each standard

(pavStdSum=ddply(subset(dat,Task=="STANDARD"&substring(SampleName,1,3)=="PAV"&TargetName=="PAV"),.(qPCRGroup,SampleName),summarise,meanCT=mean(CT,na.rm=T),meanQ=mean(Quantity),nReps=sum(!is.na(CT))))
# Groups with 0 1e3: 21
# Groups with 1 1e3: 11, 12, 13, 22
# Groups with 2 1e3: 15, 16, 18, 20

(rpvStdSum=ddply(subset(dat,Task=="STANDARD"&substring(SampleName,1,3)=="RPV"&TargetName=="RPV"),.(qPCRGroup,SampleName),summarise,meanCT=mean(CT,na.rm=T),meanQ=mean(Quantity),nReps=sum(!is.na(CT))))
# Groups with 2 1e3: 11, 22


## Re-calculate standard curves, max and min values, and sample values for runs with contamination inside the SC after curve has been modified


# Add a column indicating this
dat$contIn=with(dat,ifelse(qPCRGroup%in%c(6,11,21,22),"RPV",ifelse(qPCRGroup==4,"both","neither")))

reCFun=function(grps,df){
	
	for(i in 1:length(grps)){
		grpDat=subset(df,qPCRGroup==grps[i]&((Task=="STANDARD"&TargetName=="PAV"&substring(SampleName,1,3)=="PAV"&CT<pcon)|(Task=="STANDARD"&TargetName=="RPV"&substring(SampleName,1,3)=="RPV"&CT<rcon)))
		
	#stdReps=ddply(grpDat,.(SampleName),summarise,nreps=length(unique(Well)))
	
	#stdRepsKeep=subset(stdReps,nreps>1)
	
	#grpDat=subset(grpDat,SampleName%in%stdRepsKeep$SampleName)
		
	print(grpDat)
		
	newRMin=min(subset(grpDat,TargetName=="RPV"&!is.na(CT))$CT)
	newRMinQ=min(subset(grpDat,TargetName=="RPV"&!is.na(CT))$Quantity)
	newPMin=min(subset(grpDat,TargetName=="PAV"&!is.na(CT))$CT)
	newPMinQ=min(subset(grpDat,TargetName=="PAV"&!is.na(CT))$Quantity)
	
	newRMax=max(subset(grpDat,TargetName=="RPV"&!is.na(CT))$CT)
	newRMaxQ=max(subset(grpDat,TargetName=="RPV"&!is.na(CT))$Quantity)
	newPMax=max(subset(grpDat,TargetName=="PAV"&!is.na(CT))$CT)
	newPMaxQ=max(subset(grpDat,TargetName=="PAV"&!is.na(CT))$Quantity)
	
	df$rmin=with(df,ifelse(qPCRGroup==grps[i],newRMin,rmin))
	df$rminQ=with(df,ifelse(qPCRGroup==grps[i],newRMinQ,rminQ))
	df$pmin=with(df,ifelse(qPCRGroup==grps[i],newPMin,pmin))
	df$pminQ=with(df,ifelse(qPCRGroup==grps[i],newPMinQ,pminQ))
	
	df$rmax=with(df,ifelse(qPCRGroup==grps[i],newRMax,rmax))
	df$rmaxQ=with(df,ifelse(qPCRGroup==grps[i],newRMaxQ,rmaxQ))
	df$pmax=with(df,ifelse(qPCRGroup==grps[i],newPMax,pmax))
	df$pmaxQ=with(df,ifelse(qPCRGroup==grps[i],newPMaxQ,pmaxQ))
		
	stds=stdFun(grpDat)
	print(stds)
	df$QuantNew=with(df,ifelse(qPCRGroup==grps[i]&TargetName=="PAV"&Task!="STANDARD",10^((CT-stds$pInt)/stds$pSlope),ifelse(qPCRGroup==grps[i]&TargetName=="RPV"&Task!="STANDARD",10^((CT-stds$rInt)/stds$rSlope),QuantNew)))

	}
	return(df)
}

# Set up new column
dat$QuantNew=dat$Quantity

# Set up groups with contamination inside curve
contGrps=c(4,6,11,21,22)

# Apply function
dat2=reCFun(contGrps,dat)

# Efficiencies previously outside range
# Group 4: RPV eff=66 (worse)
# Group 6: RPV eff=96 (fixed)
# Group 11: PAV eff=82 (same, 76 when you remove the 1e3 replicate (code commented out above because it only really affects this group))
# Group 21: PAV eff=83 (same)
# Group 22: no value RPV because all standards are below the contamination

# Make sure it worked appropriately
unique(subset(subset(dat2,Quantity!=QuantNew),select=c(qPCRGroup,TargetName))) # correct groups, RPV is missing for 22
head(subset(dat2,qPCRGroup==6&TargetName=="PAV"&Quantity!=QuantNew)) # samples with PAV re-calculated are slightly different because of rounding of the slope and intercept

# Change PAV re-calculated PAV quantities back for runs without PAV contamination in the the SC (but did have RPV)
dat2$QuantNew=with(dat2,ifelse(contIn=="RPV"&TargetName=="PAV",Quantity,QuantNew))
unique(subset(subset(dat2,Quantity!=QuantNew),select=c(qPCRGroup,TargetName))) # Fixed 

# For RPV in group 22, RPV contamination was higher than standard curve (CT of 14)
subset(dat2,qPCRGroup==22&TargetName=="RPV"&CT<rcon) # none more conc.
nrow(subset(dat2,qPCRGroup==22&TargetName=="RPV"&CT>rcon)) # Lose 77 measurements
nrow(subset(dat2,qPCRGroup==22&TargetName=="RPV"&is.na(CT))) # 18 are NA
unique(subset(dat2,qPCRGroup==22&TargetName=="RPV"&CT>rcon)$QuantNew)
# Make a column to remove quantities when desired, currently, they're NA
dat2$remQuant=with(dat2,ifelse(qPCRGroup==22&TargetName=="RPV",1,0))

# Check that min and max re-calculations worked
unique(subset(subset(dat,qPCRGroup%in%c(4,6,11,21,22)),select=c(qPCRGroup,rmin,rmax,rminQ,rmaxQ,pmin,pmax,pminQ,pmaxQ)))
unique(subset(subset(dat2,qPCRGroup%in%c(4,6,11,21,22)),select=c(qPCRGroup,rmin,rmax,rminQ,rmaxQ,pmin,pmax,pminQ,pmaxQ)))

# Make 22's max and min NA instead of Inf (easier to deal with)
dat2$rmin=with(dat2,ifelse(qPCRGroup==22,NA,rmin))
dat2$rminQ=with(dat2,ifelse(qPCRGroup==22,NA,rminQ))
dat2$rmax=with(dat2,ifelse(qPCRGroup==22,NA,rmax))
dat2$rmaxQ=with(dat2,ifelse(qPCRGroup==22,NA,rmaxQ))



## Sample detected? 
# Don't worry about contamination at this point, because I didn't with the main experiment
dat2$det=as.numeric(with(dat2,ifelse(!is.na(CT),1,0)))
dat2$det


## Work with sample-only data 

# Subset for just samples
unique(dat2$SampleName)
dat3=subset(dat2,substring(SampleName,1,1)=="S")
unique(dat3$SampleName)

# Merge with base data
dat3=rename(dat3,c("SampleName"="TubeLabel"))
dat4=merge(dat3,base,all.x=T)

# Check that it worked
nrow(dat3) # 2586
nrow(dat4) # 2586
dat4$Set # There are some NA's
subset(dat4,is.na(Set))

#### removed renaming tubes (put at top) ####

# Merge again
dat4=merge(dat3,base,all.x=T)

# Check that it worked
nrow(dat3) # 2586
nrow(dat4) # 2586
sum(is.na(dat4$Set)) # No NA's
sum(is.na(dat4$SampleID)) # No NA's


## See how many samples exceed their highest standard
highSampsP=subset(dat4,TargetName=="PAV"&CT<pmin)
nrow(highSampsP) # None
highSampsR=subset(dat4,TargetName=="RPV"&CT<rmin)
nrow(highSampsR) # 105
highSampsRSum=ddply(highSampsR,.(qPCRGroup,TubeLabel),summarise,CTMean=mean(CT),n=length(TubeLabel),qMean=mean(Quantity))
nrow(highSampsRSum) # 45 samples
highSampsRSum2=subset(highSampsRSum,n>1)
nrow(highSampsRSum2)		#34 RPV samples have multiple reps with higher CT's than max
highSampsREarly=subset(highSampsRSum,qPCRGroup<4)
highSampsREarly				#14 samples were tested without 1e8. Add to last qPCR Run. If their values don't change much, no need to do the other 20.


## Check to see how the re-do of high samples went

subset(highSampsR,TubeLabel%in%highSampsREarly$TubeLabel&qPCRGroup>=4) # None of them showed up again

# Isolate the samples from the second and first runs
highSampsRReDo=subset(dat4,TubeLabel%in%highSampsREarly$TubeLabel&qPCRGroup>=4)
unique(highSampsRReDo$qPCRGroup) # All in group 21
highSampsRFirst=subset(dat4,TubeLabel%in%highSampsREarly$TubeLabel&qPCRGroup<4)

# Take the mean for each and merge
redoSum=ddply(highSampsRReDo,.(TubeLabel,TargetName),summarise,meanCT.redo=mean(CT,na.rm=T),nCT.redo=length(CT),meanQ.redo=mean(QuantNew,na.rm=T))
redoSum
firstSum=ddply(highSampsRFirst,.(TubeLabel,TargetName),summarise,meanCT.first=mean(CT,na.rm=T),nCT.first=length(CT),meanQ.first=mean(QuantNew,na.rm=T))
firstSum
rMaxSum=merge(redoSum,firstSum,all=T)
rMaxSum

# Compare values
ggplot(rMaxSum,aes(x=meanCT.redo,y=meanCT.first))+geom_point()+facet_wrap(~TargetName,scales="free") # The CT's are correlated, but higher for the redo runs than the first
ggplot(rMaxSum,aes(x=meanQ.redo,y=meanQ.first))+geom_point()+facet_wrap(~TargetName,scales="free") # For the redo, the values were estimated to be much smaller for RPV and higher for PAV
subset(subset(rMaxSum,TargetName=="RPV"),select=c(meanQ.redo,meanQ.first))

# Look at samples that weren't redone
subset(highSampsRSum,qPCRGroup>=4) # All are very close to 1e8
nrow(subset(highSampsRSum,qPCRGroup>=4)) #31 samples
nrow(subset(highSampsRSum,qPCRGroup>=4&n>1)) # 20 of them had more than one rep above the standard curve
# Leave them in the analysis. Can try without them

# Flag samples that are above the standard curve
nrow(highSampsR) # 105
dat4$highSamps=with(dat4,ifelse(CT<rmin&TargetName=="RPV",1,0))
sum(dat4$highSamps, na.rm=T) #105
sum(is.na(dat4$highSamps)) # 61 NA's
subset(subset(dat4,is.na(highSamps)),select=c(qPCRGroup,TubeLabel,TargetName,CT))# All in group 22, one that has an CT of NA that's in group 19 - change the NA to 0 (not too high) and the 22 to 1 (so that other reps decide whether all are high)
dat4$highSamps[is.na(dat4$CT)&dat4$qPCRGroup==19]=0
dat4$highSamps[dat4$qPCRGroup==22&dat4$TargetName=="RPV"]=1
# Changed above on 9/24/17 - used to be set to NA (427 samples)
# All RPV samps in qPCR group 22 are considered above curve (because there is no curve when contamination is taken into account)
sum(is.na(dat4$highSamps)) # 0
sum(subset(dat4,!is.na(highSamps))$highSamps) # 165
165-nrow(subset(dat4,qPCRGroup==22&TargetName=="RPV")) #105 (same as above)


## See how many samples are below the lowest standard
lowSampsP=subset(dat4,TargetName=="PAV"&CT>pmax)
nrow(lowSampsP) # 96
lowSampsR=subset(dat4,TargetName=="RPV"&CT>rmax)
nrow(lowSampsR) # 23

lowSampsPSum=ddply(lowSampsP,.(qPCRGroup,TubeLabel),summarise,CTMean=mean(CT),n=length(TubeLabel),qMean=mean(Quantity))
lowSampsPSum #57 samples affected
nrow(subset(lowSampsPSum,n==3)) #11 with 3 reps
lowSampsRSum=ddply(lowSampsR,.(qPCRGroup,TubeLabel),summarise,CTMean=mean(CT),n=length(TubeLabel),qMean=mean(Quantity))
lowSampsRSum #10 samples affected, 5 with 3 reps

# Flag samples with low values
dat4$lowSamps=with(dat4,ifelse((CT>rmax&TargetName=="RPV")|(CT>pmax&TargetName=="PAV"),1,0))
sum(is.na(dat4$lowSamps)) # 486 are NA's
unique(subset(subset(dat4,is.na(lowSamps)),select=c(qPCRGroup,TargetName))) # some are RPV group 22
head(subset(dat4,is.na(lowSamps)&TargetName=="RPV"&qPCRGroup==19)) # the rest are probably NA
dat4$lowSamps[is.na(dat4$CT)|(dat4$qPCRGroup==22&dat4$TargetName=="RPV")]=1
sum(is.na(dat4$lowSamps)) # 0
sum(dat4$lowSamps) # 605
sum(is.na(dat4$CT))+nrow(subset(dat4,qPCRGroup==22&TargetName=="RPV"))+96+23-nrow(subset(dat4,qPCRGroup==22&TargetName=="RPV"&is.na(CT))) # all are accounted for

# Create tube labe with target
dat4$TubeTarget=with(dat4,paste(TubeLabel,TargetName,sep="."))

## Mark duplicate testings

# See if each sample was once for PAV and once RPV
sampList=ddply(dat4,.(qPCRGroup,TubeLabel),summarise,targets=length(unique(TargetName)),highSamps=sum(highSamps,na.rm=T))
sampList 
sum(sampList$targets!=2)# Yes

# Mark duplicates
sampList$Dups=duplicated(sampList$TubeLabel)
sum(sampList$Dups==T) # 34
dupIDs=subset(sampList,Dups==T)$TubeLabel
dupIDs
unique(dupIDs)
dupList=ddply(subset(sampList,TubeLabel%in%dupIDs),.(TubeLabel),summarise,qGroups=paste(qPCRGroup,collapse=","),highS=paste(highSamps,collapse=","))
dupList # all in groups 21 or 22
# The group 7 one should be S1 LR2S2 - go back up and change - done
dat4$Dups=with(dat4,ifelse(TubeLabel%in%dupIDs,1,0))
# The duplicates in 21 were re-done because they were above the standard curve
# However, groups 4, 11, 21, and 22 all had poor efficiency curves and contamination (6 was fixed when potentially contaminated samples were removed)
# RPV for 22 is totally unreliable - quantities converted to NA below

# Group 22 PAV
# 22 and 11 only have one low PAV standard
# the efficiency for 11 is outside of the range
grp22Dups=subset(dat4,Dups==1&qPCRGroup==22&TargetName=="PAV")
length(unique(grp22Dups$TubeTarget)) #20 samples
grp22Dups
# some are low samples
subset(dat4,TubeTarget%in%grp22Dups$TubeTarget&qPCRGroup!=22)
# some are low samples
unique(subset(subset(dat4,TubeTarget%in%grp22Dups$TubeTarget&lowSamps==1),select=c(qPCRGroup,TubeTarget)))
# 1 only low in 22, the rest low in both runs
grp22Mean=ddply(subset(subset(dat4,TubeTarget%in%grp22Dups$TubeTarget),select=c(TubeTarget,QuantNew,qPCRGroup,Well)),.(TubeTarget,qPCRGroup),summarise,meanQ=mean(QuantNew,na.rm=T))
grp22Mean$qPCRGroup=as.factor(grp22Mean$qPCRGroup)
grp22Mean$qPCRGroup=revalue(grp22Mean$qPCRGroup,c("6"="non22","11"="non22","4"="non22"))
grp22Wide=reshape(grp22Mean,direction="wide",idvar="TubeTarget",v.names="meanQ",timevar="qPCRGroup")
grp22Wide$meanQ.non22[is.na(grp22Wide$meanQ.non22)]=0
grp22Wide$meanQ.22[is.na(grp22Wide$meanQ.22)]=0
ggplot(grp22Wide,aes(x=meanQ.non22,y=meanQ.22))+geom_point()+xlim(c(0,35000))+ylim(c(0,35000))
#values are all higher coming out of 22, but they are relatively correlated
#remove all from 22
dat4$remQuant=with(dat4,ifelse(TubeTarget%in%grp22Dups$TubeTarget&qPCRGroup==22,1,remQuant))

# Deicde how to deal with group 21 duplicates
# 21 is the only group wtih no PAV 1e3 (makes low values unreliable)
# poor efficiency for PAV, contamination with RPV that's within SC
grp21Dups=subset(dat4,Dups==1&qPCRGroup==21)$TubeTarget
subset(dat4,TubeTarget%in%grp21Dups&qPCRGroup!=21)
# some duplicates have no issues, some are high, and some are low

subset(subset(dat4,TubeTarget%in%grp21Dups&lowSamps==1),select=c(TubeLabel,TargetName,qPCRGroup,CT)) #the ones that were low samples in their other group were also low in 21 because they were not detected in both groups
# there are low samples in 21 that aren't low in the other groups
# the low samples should come from their non-21 group
subset(dat4,qPCRGroup==21&lowSamps==1&!(TubeTarget%in%grp21Dups))
# 4 samples have values below SC and are in 21 (no low SC values) do not have replacements 
# 2 of the 4 are not detected
# keep these 4 in

subset(subset(dat4,TubeTarget%in%grp21Dups&highSamps==1),select=c(TubeLabel,TargetName,qPCRGroup,CT)) #there are no high samples in 21
# From re-do analysis above, the RPV estimates were lower in grp 21 and PAV were higher
grp21HighSamps=subset(dat4,TubeTarget%in%grp21Dups&highSamps==1)$TubeTarget
subset(subset(dat4,TubeTarget%in%grp21HighSamps&qPCRGroup==21),select=c(TubeTarget,qPCRGroup,CT,rcon))
# all of the RPV samples have CT's more concentrated (lower) than the contamination
# I think I want to keep the group 21 samples in this case

# Duplicated samples:
# Keep from non 21 if there's no issue
# Keep from non 21 if they are below the SC in 21
# Keep from 21 if they're low (besides above) or high

remSamps21=subset(dat4,TubeTarget%in%grp21Dups&((qPCRGroup!=21&lowSamps==0&highSamps==0)|(qPCRGroup==21&lowSamps==1)))$TubeTarget
remSampsNon21=subset(dat4,TubeTarget%in%grp21Dups&(lowSamps==1|highSamps==1)&!(TubeTarget%in%remSamps21))$TubeTarget
length(unique(c(remSamps21,remSampsNon21)))
length(unique(grp21Dups))

# Mark samples for removal
dat4$remQuant=with(dat4,ifelse((TubeTarget%in%remSamps21&qPCRGroup==21)|(TubeTarget%in%remSampsNon21&qPCRGroup!=21),1,remQuant))

# Make sure all duplicates are represented
length(unique(subset(dat4,qPCRGroup!=22|TargetName=="PAV")$TubeTarget)) #794
length(unique(subset(dat4,remQuant==0)$TubeTarget)) #794


# Mark contamination and poor efficiency groups
dat4$PoorSCEff=with(dat4,ifelse((qPCRGroup==4&TargetName=="RPV")|(qPCRGroup%in%c(11,21)&TargetName=="PAV"),1,0))


# Check S2 PP5S1. The cage broke off the leaf after IAP, and I didn't mean to extract this one.
subset(dat4,TubeLabel=="S2 PP5S1")
# Remove to be consistent with analysis decisions. Also, it's supposed to be a singly infected PAV plant and RPV is present
dat5=subset(dat4,TubeLabel!="S2 PP5S1")

# See which samples have all three reps fitting criteria of low samps, high samps, or undetected
dat6<-ddply(dat5,.(TubeLabel,TargetName,qPCRGroup),mutate,nSamps=length(TubeLabel),nLow=sum(lowSamps),nHigh=sum(highSamps),nCTNA=sum(is.na(CT)))
# Did within qPCR group because I get rid of the duplicate groups later

# See if any of the above are NA
sum(is.na(dat6$nLow))
sum(is.na(dat6$nHigh))

# For samples with all three reps below the standard curve (lowSamps), want to change their values to 0 instead of NA for some of the quantities (below). 
dat6$AllBelow=with(dat6,ifelse(nLow==nSamps,1,0))
nrow(subset(dat6,AllBelow==1)) # 561 (include undetected)

# For samples with all three reps undetected, want to change their values to 0 instead of NA for all of the quantities (below).
dat6$AllCTNA=with(dat6,ifelse(nCTNA==nSamps,1,0))
nrow(subset(dat6,AllCTNA==1)) # 387

# How many have all three reps as high
dat6$AllHigh=with(dat6,ifelse(nHigh==nSamps,1,0))
nrow(subset(dat6,AllHigh==1)) # 138
nrow(subset(dat6,AllHigh==1&remQuant==0)) #39
length(unique(subset(dat6,AllHigh==1&remQuant==0)$TubeTarget)) #13 samples
subset(subset(dat6,AllHigh==1&remQuant==0),select=c(TubeTarget,qPCRGroup,CT,rmin,QuantNew,rmaxQ))
# all are very close to max quantity of SC, within 1/2 an order of magnitude larger
subset(subset(dat6,AllHigh==0&highSamps==1&remQuant==0),select=c(TubeTarget,qPCRGroup,CT,rmin,QuantNew,rmaxQ))
# these ones are also vary close
subset(subset(dat6,AllHigh==0&highSamps==1&remQuant==0),select=c(TubeTarget,qPCRGroup,CT,rmin,QuantNew,rmaxQ))
# so are these
someHigh=subset(dat6,AllHigh==0&highSamps==1&remQuant==0)$TubeTarget
ggplot(subset(dat6,TubeTarget%in%someHigh),aes(x=TubeTarget,y=QuantNew))+geom_point(aes(colour=as.factor(highSamps)))
# all the reps are close together - leave high samples in analysis
length(unique(subset(dat6,highSamps==1&remQuant==0)$TubeTarget)) #32 samples affected - keep them all in


### Add quantity columns

## Factors to consider
# contamination (contIn, both, RPV, neither)
# poor efficiency (PoorSCEff, 1/0)
# duplicates (Dups, 1/0)
# low or high samples relative to standard curve (lowSamps, highSamps, 1/0/NA)

## Quantity 1 
# Contamination accounted for (removed stds that could have been contaminated), poor eff okay, all detected quantities used, 0 if all reps undetected
dat6$Quant1=dat6$QuantNew
unique(subset(subset(dat6,AllCTNA==1),select=c(qPCRGroup,TargetName,Quant1)))
dat6$Quant1=with(dat6,ifelse(AllCTNA==1,0,Quant1))

# Remove RPV from group 22 and duplicates with better values
dat6$Quant1=with(dat6,ifelse(remQuant==1,NA,Quant1))

## Quantity 1.GE 
# Contamination accounted for, poor eff removed, all detected quantities used
dat6$Quant1.GE=with(dat6,ifelse(PoorSCEff==0,Quant1,NA))

## Quantity 2 
# Cont accounted for, poor eff okay, samples below curve set to NA, unless all reps are below, then they are zero, except for RPV grp 22 
dat6$Quant2=with(dat6,ifelse(lowSamps==1,NA,Quant1))
dat6$Quant2=with(dat6,ifelse(AllBelow==1,0,Quant2))

# Remove RPV from group 22 and duplicates with better values
dat6$Quant2=with(dat6,ifelse(remQuant==1,NA,Quant2))

## Quantity 2.GE 
# No cont, poor eff removed, samples outside curve set to NA
dat6$Quant2.GE=with(dat6,ifelse(PoorSCEff==0,Quant2,NA))

## Double check RPV group 22
subset(dat6,qPCRGroup==22)

## Export full dataset
setwd("/Users/AmyKendig/Google Drive/Within Host Coexistence/Data")
write.csv(dat6,"Expt2_qPCRDataEdited_050318.csv",row.names=F)