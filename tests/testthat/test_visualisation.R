# Test for the function plot_seroprev_fitted

library(testthat)

test_that("individual models", {
  # So far we are skipping tests on these platforms until
  # we find an efficient way to update rstan testthat snapshots on all of them
  skip_on_os(c("windows", "mac"))
  source("testing_utils.R")
  set.seed(1234) # For reproducibility

  library(devtools)
  library(dplyr)
  library(vdiffr)
  library(jsonlite)

  data("simdata_large_epi")
  simdata <- simdata_large_epi %>% prepare_serodata()
  no_transm <- 0.0000000001
  big_outbreak <- 1.5
  foi_sim <- c(rep(no_transm, 32), rep(big_outbreak, 3), rep(no_transm, 15)) # 1 epidemics


  #----- Results visualisation
  size_text <- 6
  max_lambda <- 1.55

  model_constant_json <- jsonlite::fromJSON(testthat::test_path("extdata", "model_constant.json"))
  model_constant <- jsonlite::unserializeJSON(model_constant_json)
  constant_plot <- plot_seromodel(model_constant,
    size_text = size_text,
    max_lambda = max_lambda,
    foi_sim = foi_sim
  )

  model_tv_normal_json <- fromJSON(testthat::test_path("extdata", "model_tv_normal.json"))
  model_tv_normal <- jsonlite::unserializeJSON(model_tv_normal_json)
  tv_normal_plot <- plot_seromodel(model_tv_normal,
    size_text = size_text,
    max_lambda = max_lambda,
    foi_sim = foi_sim
  )

  model_tv_normal_log_json <- fromJSON(testthat::test_path("extdata", "model_tv_normal_log.json"))
  model_tv_normal_log <- jsonlite::unserializeJSON(model_tv_normal_log_json)
  tv_normal_log_plot <- plot_seromodel(model_tv_normal_log,
    size_text = size_text,
    max_lambda = max_lambda,
    foi_sim = foi_sim
  )

  plot_arrange <- cowplot::plot_grid(constant_plot,
    tv_normal_plot,
    tv_normal_log_plot,
    ncol = 3, labels = "AUTO"
  )
  vdiffr::expect_doppelganger("plot_arrange_simdata_foi", plot_arrange)
})