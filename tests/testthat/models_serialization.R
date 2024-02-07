library(serofoi)
library(testthat)

set.seed(1234) # For reproducibility

#----- Read and prepare data
data("simdata_large_epi")
simdata <- prepare_serodata(simdata_large_epi)
no_transm <- 0.0000000001
big_outbreak <- 1.5
foi_sim <- c(rep(no_transm, 32), rep(big_outbreak, 3), rep(no_transm, 15)) # 1 epidemics


#----- Run models
models_to_run <- c(
  "constant",
  "tv_normal",
  "tv_normal_log"
)

models_list <- lapply(models_to_run,
  run_seromodel,
  serodata = simdata,
  iter = 1000
)

saveRDS(
  models_list[[1]],
  testthat::test_path("extdata", "model_constant.RDS")
)

saveRDS(
  models_list[[2]],
  testthat::test_path("extdata", "model_tv_normal.RDS")
)

saveRDS(
  models_list[[3]],
  testthat::test_path("extdata", "model_tv_normal_log.RDS")
)
