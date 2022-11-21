# compile_stan_files

rm(list = ls())
library(rstan)
library(tidyverse)


ConstantUniformFOI     <- rstan::stan_model('inst/extdata/stanmodels/constant_foi_Bi.stan')
saveRDS(ConstantUniformFOI, 'inst/extdata/stanmodels/ConstantUniformFOI.RDS')

ContinuousNormalFOI    <- rstan::stan_model('inst/extdata/stanmodels/continuous_foi_normal_Bi.stan')
saveRDS(ContinuousNormalFOI, 'inst/extdata/stanmodels/ContinuousNormalFOI.RDS')

ContinuousNormalLogFOI_lowt  <- rstan::stan_model('inst/extdata/stanmodels/continuous_foi_normal_log.stan')
saveRDS(ContinuousNormalLogFOI_lowt, 'inst/extdata/stanmodels/ContinuousNormalLogFOI_lowt.RDS')
