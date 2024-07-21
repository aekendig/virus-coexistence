## Goal: create EML metadata for project using the Environmental Data Initiative's EML Assembly Line (https://github.com/EDIorg/EMLassemblyline)

# Put data and code in a folder together to be grabbed by make_eml
# make an empty metadata folder
# Generate metadata files by editing current ones or generating them (see Github page for tutorial)
# Edit and run this script


#### set up ####

# clear environment
rm(list=ls())

# load libraries
library(EMLassemblyline)
library(tidyverse)
library(knitr)


#### import templates ####

# list of data files
dlist <- list.files(path = "./edi",
                    pattern = ".csv")

# create high level text files
template_core_metadata(path = "./metadata",
                       license = "CCBY")

# create an attribute file for each data table
template_table_attributes(path = "./metadata",
                          data.path = "./edi",
                          data.table = dlist)

# create a categorical value file for each data table
template_categorical_variables(path = "./metadata",
                               data.path = "./edi")

# look at units
view_unit_dictionary()


#### data file descriptors ####

dlist

# description list
ddlist <- NA

ddlist[1] <- "experiment information at the scale of a plant"
ddlist[2] <- "raw qPCR sample data"

# name list
dnlist <- c("Plant data",
            "qPCR data")

# print table
# dtable <- data.frame(data = dlist, description = ddlist)
# kable(dtable)


#### code descriptors ####

# list of code files
clist <- c(list.files(path = "./edi",
                      pattern = ".R"))

# re-arrange code list, omit data files
clist2 <- clist[c(6, 4, 2, 3, 1)]

# code descriptions
cdlist <- c(
  "code to process qPCR data",
  "code to analyze qPCR data",
  "code to estimate simulation model parameters",
  "code to run model simulations",
  "code for figures and model settings, referenced by other scripts"
)

# name list
cnlist <- c("Data processing",
            "Experiment data analysis",
            "Model parameters",
            "Model simulations",
            "Model settings")

# print table
# ctable <- data.frame(code = clist2, desription = cdlist)
# kable(ctable)


#### make eml ####

# copied data and code from the respective folders and put into edi folder

make_eml(path = "./metadata",
         data.path = "./edi",
         dataset.title = "Within-plant coexistence of viruses across nitrogen and phosphorus supply rates",
         data.table = dlist,
         data.table.name = dnlist,
         data.table.description = ddlist,
         data.table.quote.character = rep("\"", length(dlist)),
         other.entity = clist2,
         other.entity.name = cnlist,
         other.entity.description = cdlist,
         temporal.coverage = c("2015-08-29", "2015-12-21"),
         geographic.description = "St. Paul, MN, USA",
         geographic.coordinates = c(44.98, -93.18, 44.98, -93.18),
         maintenance.description = "completed", 
         user.id = "aekendig",
         user.domain = "EDI",
         package.id = "edi.1602.2")


#### check warnings ####

eml <- EML::read_eml("./metadata/edi.1602.2.xml")
EML::eml_validate(eml)
