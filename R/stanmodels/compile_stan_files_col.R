# compile_stan_files 

rm(list=ls())
library(rstan)
library(tidyverse)


ConstantUniformFOI     <- stan_model('stanmodels-col/constant_foi_Bi.stan')
saveRDS(ConstantUniformFOI, 'stanmodels-col/ConstantUniformFOI.RDS')

ContinuousNormalFOI    <- stan_model('stanmodels-col/continuous_foi_normal_Bi.stan')
saveRDS(ContinuousNormalFOI, 'stanmodels-col/ContinuousNormalFOI.RDS')

ContinuousNormalLogFOI_lowt  <- stan_model('stanmodels-col/continuous_foi_normal_log.stan')
saveRDS(ContinuousNormalLogFOI_lowt, 'stanmodels-col/ContinuousNormalLogFOI_lowt.RDS')

#For details, look into MS_word file