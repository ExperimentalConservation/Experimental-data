## eme scheme trial 2
rm(list=ls())

library(devtools)

##install various packages
##devtools::install_github("Exp-Micro-Ecol-Hub/dmdScheme", ref = "master", build_opts = c("--no-resave-data"))
##devtools::install_github("Exp-Micro-Ecol-Hub/emeScheme", ref = "master", build_opts = c("--no-resave-data"))

library(dmdScheme)
library(emeScheme)

#open_new_spreadsheet(keepData=T, format=F)

##set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

x<-validate("EXP 1 meta.xlsx")

report('EXP 1 meta.xlsx', '~/Desktop/Work/Experimental-data/Clements et al. (2013) JAE', report = "html")

