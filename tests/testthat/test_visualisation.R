# Test for the function plot_seroprev_fitted

library(testthat)

test_that("individual models", {
  set.seed(1234) # For reproducibility

  library(devtools)
  library(vdiffr)
  library(jsonlite)

  data(simdata_large_epi)
  simdata <- prepare_serodata(simdata_large_epi)
  no_transm <- 0.0000000001
  big_outbreak <- 1.5
  foi_sim <- c(rep(no_transm, 32), rep(big_outbreak, 3), rep(no_transm, 15)) # 1 epidemics


  #----- Results visualisation
  size_text <- 6
  max_lambda <- 1.55

  model_constant <- readRDS(testthat::test_path("extdata", "model_constant.RDS"))
  constant_plot <- plot_seromodel(seromodel_object = model_constant,
                                  size_text = size_text,
                                  max_lambda = max_lambda,
                                  foi_sim = foi_sim
  )

  model_tv_normal <- readRDS(testthat::test_path("extdata", "model_tv_normal.RDS"))
  tv_normal_plot <- plot_seromodel(seromodel_object = model_tv_normal,
                                   size_text = size_text,
                                   max_lambda = max_lambda,
                                   foi_sim = foi_sim
  )

  model_tv_normal_log <- readRDS(testthat::test_path("extdata", "model_tv_normal_log.RDS"))
  tv_normal_log_plot <- plot_seromodel(seromodel_object = model_tv_normal_log,
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