# compile_stan_files

rm(list=ls())
library(rstan)
library(tidyverse)


ConstantUniformFOI     <- stan_model('R/stanmodels/constant_foi_Bi.stan')
saveRDS(ConstantUniformFOI, 'R/stanmodels/ConstantUniformFOI.RDS')

ContinuousNormalFOI    <- stan_model('R/stanmodels/continuous_foi_normal_Bi.stan')
saveRDS(ContinuousNormalFOI, 'R/stanmodels/ContinuousNormalFOI.RDS')

ContinuousNormalLogFOI_lowt  <- stan_model('R/stanmodels/continuous_foi_normal_log.stan')
saveRDS(ContinuousNormalLogFOI_lowt, 'R/stanmodels/ContinuousNormalLogFOI_lowt.RDS')
