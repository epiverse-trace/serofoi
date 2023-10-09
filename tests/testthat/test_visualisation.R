library(testthat)

#Test for exception in plot_rhats
test_that("plot_rhats function works for else statement", {

  seromodel_object <- "model did not run"
  rhats_plot <- plot_rhats(seromodel_object)

  # Expect the output to be a ggplot object
  expect_is(rhats_plot, "ggplot")

  # Expect the plot to have no points
  expect_equal(length(rhats_plot$layers[[1]]$data), 0)

  # Expect the plot to have a single label per axis
  expect_equal(length(rhats_plot$layers[[2]]$data), 2)
})

# Unit tests for plot_seroprev_fitted()
test_that("plot_seroprev_fitted() returns an empty plot for an unfitted seromodel object", {
  # We create an unadjusted seromodel object
  seromodel <- "dummy"
  serodata <- data.frame(age = c(0, 10, 20, 30, 40, 50, 60),
                                   p_obs_bin = c(0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99),
                                   bin_size = c(5, 10, 20, 30, 40, 50, 60))
  # We call the function plot_seroprev_fitted()
  seroprev_fitted <- plot_seroprev_fitted(seromodel, serodata)
  # We verify that the plot object is an empty ggplot
  expect_true(class(seroprev_fitted)[2] == "ggplot")
  expect_true(length(seroprev_fitted$layers) == 2)
})

