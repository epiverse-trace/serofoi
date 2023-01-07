#####################----SEROFOI----   #################
#
#   Seroprevalence Force-of-Infection
#   Author:  Zulma M Cucunubá
#   Version 1.0.1. This version runs and compares three basic models

####################################################


rm(list = ls())
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


# ---- Models
# For these to run well, I need to re-compile the RDSs files
# source('stanmodels-col/compile_stan_files_col.R')



#----- Read Data
dat0 <- readRDS("data/data.RDS")


# ----  Choose Models
Model0   <- readRDS("R/stanmodels/ConstantUniformFOI.RDS")
Model1   <- readRDS("R/stanmodels/ContinuousNormalFOI.RDS")
Model2   <- readRDS("R/stanmodels/ContinuousNormalLogFOI_lowt.RDS")




# Automated name of the folder where results will be stored
my_dir <- epitrix::clean_labels(paste0("tests_", Sys.time()))
dir_results(my_dir)


print(paste0("my results will be sortored at:_________test/", my_dir))



i <- dat0$survey[1]
run_save_models(my_dir    = my_dir,
              survey       = i,
              dat0      = dat0,
              n_iters   = 3000,
              Model0 = Model0, NameModel0 = "M0_Constant",
              Model1 = Model1, NameModel1 = "M1_Cont_Normal",
              Model2 = Model2, NameModel2 = "M2_Cont_NormalLog"
)
