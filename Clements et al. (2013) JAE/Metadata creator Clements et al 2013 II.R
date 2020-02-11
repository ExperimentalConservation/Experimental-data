## eme scheme trial 2
rm(list=ls())

library(devtools)

##install various packages
  ##devtools::install_github("Exp-Micro-Ecol-Hub/dmdScheme", ref = "master", build_opts = c("--no-resave-data"))
  ##devtools::install_github("Exp-Micro-Ecol-Hub/emeScheme", ref = "master", build_opts = c("--no-resave-data"))
devtools::install_github("Exp-Micro-Ecol-Hub/dmdScheme", ref = "prepareReleaseCRAN", build_opts = c("--no-resave-data"))

library(dmdScheme)
#library(emeScheme)

#open_new_spreadsheet(keepData=T, format=F)

##set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

read_excel("EXP 1 meta.xlsx", validate = T)

scheme_make("EXP 1 meta.xlsx", overwrite=T)

x<-validate("EXP 1 meta.xlsx", errorIfStructFalse=F)

report('EXP 1 meta.xlsx', '~/Desktop/Work/Experimental-data/Clements et al. (2013) JAE', report = "html", errorIfStructFalse=F)

