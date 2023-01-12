#####################----SEROFOI----   #################
#
#   Seroprevalence Force-of-Infection
#   Author:  Zulma M Cucunubá
#   Version 1.0.1. This version runs and compares three basic models

####################################################
rm(list = ls())


#installlibrary


library(rstan)
library(tidyverse)
library(reshape2)
library(bayesplot)
library(loo)
library(pracma)
library(cowplot)
library(grid)
library(gridExtra)
library(Hmisc)
library(dplyr)
#library(vscDebugger)
library(epitrix)
library(gsubfn) # obtain limits from applying cut function


source("r/modeling.R")
#source("r/fitting.R")
source("r/seroprevalence_data.R")
source("r/model_comparison.R")
source("r/visualisation.R")



#----- Read Data
dat <- readRDS("data/data.RDS")




run_model(data      = dat,
          name_model = "m0_constant")
run_model(data      = dat,
          name_model = "m1_normal")
run_model(data      = dat,
          name_model = "m2_normal_log")
check_model(run_modelobject)
plot_model(run_modelobject)


run_models (data = dat,
            names_models = c("m0_constant", "m1_normal_constant"))



