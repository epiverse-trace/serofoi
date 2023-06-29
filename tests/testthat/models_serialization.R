library(devtools)
library(dplyr)
library(serofoi)
library(testthat)

set.seed(1234) # For reproducibility

#----- Read and prepare data
data("simdata_large_epi")
simdata <- simdata_large_epi %>% prepare_serodata()
no_transm <- 0.0000000001
big_outbreak   <- 1.5
foi_sim <- c(rep(no_transm, 32), rep(big_outbreak, 3), rep(no_transm, 15)) # 1 epidemics


#----- Run models
models_to_run <- c("constant",
                   "tv_normal",
                   "tv_normal_log")

models_list <- lapply(models_to_run,
                      run_seromodel,
                      serodata = simdata,
                      n_iters = 1000)

model_constant_json <- jsonlite::serializeJSON(models_list[[1]])
write_json(model_constant_json, testthat::test_path("extdata", "model_constant.json"))

model_tv_normal_json <- jsonlite::serializeJSON(models_list[[2]])
write_json(model_tv_normal_json, testthat::test_path("extdata", "model_tv_normal.json"))

model_tv_normal_log_json <- jsonlite::serializeJSON(models_list[[3]])
write_json(model_tv_normal_json, testthat::test_path("extdata", "model_tv_normal_log.json"))
