test_that("individual models", {
  # So far we are skipping tests on these platforms until
  # we find an efficient way to update rstan testthat snapshots on all of them
  skip_on_os(c("windows", "mac"))
  skip_on_ci()
  source("testing_utils.R")
  set.seed(1234) # For reproducibility

  library(devtools)
  library(dplyr)
  library(vdiffr)

  #----- Read and prepare data
  data("serodata_simD")
  simdata <- serodata_simD %>% prepare_serodata()
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

  #----- Results visualisation
  size_text <- 6
  max_lambda <- 1.55
  constant_plot <- plot_seromodel(models_list[[1]],
                                  size_text = size_text,
                                  max_lambda = max_lambda,
                                  foi_sim = foi_sim)
  tv_normal_plot <- plot_seromodel(models_list[[2]],
                                   size_text = size_text,
                                   max_lambda = max_lambda,
                                   foi_sim = foi_sim)
  tv_normal_log_plot <- plot_seromodel(models_list[[3]],
                                       size_text = size_text,
                                       max_lambda = max_lambda,
                                       foi_sim = foi_sim)
  plot_arrange <- cowplot::plot_grid(constant_plot,
                                     tv_normal_plot,
                                     tv_normal_log_plot,
                                     ncol = 3, labels = "AUTO")
  vdiffr::expect_doppelganger("plot_arrange_simdata_foi", plot_arrange)

})
