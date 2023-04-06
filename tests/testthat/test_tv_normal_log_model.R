test_that("individual models", {
  # So far we are skipping tests on these platforms until
  # we find an efficient way to update rstan testthat snapshots on all of them
  skip_on_os(c("windows", "mac"))
  source("testing_utils.R")
  set.seed(1234) # For reproducibility

  library(devtools)
  library(dplyr)
  library(vdiffr)

  #----- Read and prepare data
  data("serodata")
  data_test <- serodata %>% prepare_serodata(alpha = 0.05)

  #----- Plot raw data
  data_test_plot <- plot_seroprev(data_test, size_text = 15)
  vdiffr::expect_doppelganger("serodata_plot", data_test_plot)

  #----- Generate plots for the constant model

  model_name <- "tv_normal_log"
  model <- run_seromodel(serodata = data_test,
                         foi_model = model_name,
                         n_iters = 1000)
  model_plot <- plot_seromodel(model, size_text = 6)
  vdiffr::expect_doppelganger(paste0(model_name, "_model_plot"), model_plot)

  model_seroprev_plot <- plot_seroprev_fitted(model, size_text = 15)
  vdiffr::expect_doppelganger(paste0(model_name, "_sp_fitted_plot"), model_seroprev_plot)

  model_foi_plot <- plot_foi(model, size_text = 15)
  vdiffr::expect_doppelganger(paste0(model_name, "_foi_plot"), model_foi_plot)

  model_rhats_plot <- plot_rhats(model, size_text = 15)
  vdiffr::expect_doppelganger(paste0(model_name, "_rhats_plot"), model_rhats_plot)

})
