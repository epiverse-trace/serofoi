library(testthat)
library(ggplot2)
library(rlang)


create_expected_serosurvey <- function(actual_serosurvey) {
  return(prepare_serosurvey_for_plotting(actual_serosurvey %>%
    add_age_group_to_serosurvey()))
}

expect_correct_serosurvey_plot <- function(actual_serosurvey, expected_serosurvey, plot) {
  # Test that the output is a ggplot object
  expect_s3_class(plot, "ggplot")

  # Test that the correct data is used
  expect_equal(plot$data, expected_serosurvey)

  # Checks that the right columns are used for ymin and yman in the error bars
  expect_equal(quo_name(get_expr(plot$layers[[1]]$mapping$ymin)), "seroprev_lower")
  expect_equal(quo_name(get_expr(plot$layers[[1]]$mapping$ymax)), "seroprev_upper")
  expect_equal(quo_name(get_expr(plot$layers[[2]]$mapping$y)), "seroprev")
  expect_equal(quo_name(get_expr(plot$layers[[2]]$mapping$size)), "n_sample")

  # Check that the plot uses the correct coordinate limits
  expect_equal(
    plot$coordinates$limits$x,
    c(min(expected_serosurvey$age_min), max(expected_serosurvey$age_max))
  )
  expect_equal(
    plot$coordinates$limits$y,
    c(min(expected_serosurvey$seroprev_lower), max(expected_serosurvey$seroprev_upper))
  )

  # Test labels
  expect_equal(plot$labels$y, "Seroprevalence")
  expect_equal(plot$labels$x, "Age")
}

test_that("plot_serosurvey creates a ggplot with correct structure", {
  data(veev2012)

  actual_serosurvey <- veev2012

  expected_serosurvey <- create_expected_serosurvey(actual_serosurvey)

  plot <- plot_serosurvey(actual_serosurvey)

  expect_correct_serosurvey_plot(actual_serosurvey, expected_serosurvey, plot)
})



# test_that("plot_serosurvey creates a binned ggplot with correct structure", {
#   data(veev2012)

#   actual_serosurvey <- veev2012

#   expected_serosurvey <- create_expected_serosurvey(actual_serosurvey)

#   plot <- plot_serosurvey(actual_serosurvey, bin_serosurvey = TRUE, bin_step = 10)

#   expect_correct_serosurvey_plot(actual_serosurvey, expected_serosurvey, plot)

#   expect_equal(plot$data$n_sample, tapply(expected_serosurvey$n_sample, expected_serosurvey$age_group, sum))
# })

