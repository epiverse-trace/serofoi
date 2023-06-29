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


# #Test for the function plot_rhats
#
# library(testthat)
#
# # Define a test context
# context("Testing plot_rhats function")
#
# # Create a mock seromodel_object
# mock_seromodel_object <- list(fit = "model did not run")
#
# # Define a test case for the else statement
# test_that("plot_rhats function works for else statement", {
#
#   # Call the function with the mock seromodel_object
#   rhats_plot <- plot_rhats(mock_seromodel_object)
#
#   # Expect the output to be a ggplot object
#   expect_is(rhats_plot, "ggplot")
#
#   # Expect the plot to have a single point
#   expect_equal(length(rhats_plot$layers[[1]]$data), 1)
#
#   # Expect the plot to have a single label
#   expect_equal(length(rhats_plot$layers[[2]]$data), 1)
#
#   # Expect the label to be "errors"
#   expect_equal(rhats_plot$layers[[2]]$label, "errors")
# })
#
#
# #Test for the function plot_seromodel
#
# library(testthat)
#
# # Test for exception in else
# test_that("plot_seromodel prints an error message and returns an empty plot object when a model cannot be fitted", {
#   # Create a seromodel object with fit as a feature
#   seromodel_object <- list(fit = "no_fit", model = "my_model")
#   # Run the function and check that it returns an empty plot object
#   expect_silent(plot_seromodel(seromodel_object))
#   expect_equal(length(plot_seromodel(seromodel_object)$grobs), 5)
# })
#
# # We create a helper function that returns an unwrapped seromodel object
# create_dummy_seromodel <- function() {
#   # empty object
#   seromodel <- list()
#   seromodel$fit <- "dummy"
#   seromodel$serodata <- data.frame(age = c(0, 10, 20, 30, 40, 50, 60),
#                                    p_obs_bin = c(0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99),
#                                    bin_size = c(5, 10, 20, 30, 40, 50, 60))
#   return(seromodel)
# }
#
# # Unit tests for plot_seroprev_fitted()
# test_that("plot_seroprev_fitted() returns an empty plot for an unfitted seromodel object", {
#   # We create an unadjusted seromodel object
#   seromodel <- create_dummy_seromodel()
#   # We call the function plot_seroprev_fitted()
#   plot <- plot_seroprev_fitted(seromodel)
#   # We verify that the plot object is an empty ggplot
#   expect_true(class(plot) == "ggplot")
#   expect_true(length(plot$layers) == 0)
# })
#
#
# #Test for the function plot_foi
#
# library(testthat)
#
# # Define the test context
# context("Test of the function plot_foi")
#
# # Create a test to check the else block
# test_that("The plot_foi function should output an empty plot when the model is not running", {
#
#   # Create an empty object that simulates the output of the model that failed
#   empty_model <- list(fit = "failure")
#
#   # Run the plot_foi function with the empty model
#   plot <- plot_foi(empty_model)
#
#   # Check if the output is an empty graph
#   expect_identical(ggplot2::ggplot(), plot)
# })
#
