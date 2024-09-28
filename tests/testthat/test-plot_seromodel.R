library(testthat)
library(ggplot2)
library(rlang)
library(purrr)

extract_plot_data <- function(plot) {
  list(
    classes = class(plot),
    layers = extract_layers(plot$layers),
    coordinates = list(
      limits = plot$coordinates$limits
    ),
    labels = plot$labels
  )
}

extract_layers <- function(layers) {
  map(layers, function(layer) {
    # Extract mappings
    list(
      mapping =
        map(layer$mapping, function(x) {
          quo_name(get_expr(x))
        })
    )
  })
}


create_expected_serosurvey <- function(actual_serosurvey) {
  return(prepare_serosurvey_for_plotting(actual_serosurvey %>%
    add_age_group_to_serosurvey()))
}


test_that("plot_serosurvey creates a ggplot with correct structure", {
  data(veev2012)

  actual_serosurvey <- veev2012

  actual_plot <- extract_plot_data(plot_serosurvey(actual_serosurvey))

  expected_serosurvey <- create_expected_serosurvey(actual_serosurvey)

  expected_plot <- list(
    classes = c("gg", "ggplot"),
    layers = list(
      list(
        mapping = list(
          ymin = "seroprev_lower",
          ymax = "seroprev_upper"
        )
      ),
      list(
        mapping = list(
          y = "seroprev",
          size = "n_sample"
        )
      )
    ),
    coordinates = list(
      limits = list(
        x = c(min(expected_serosurvey$age_min), max(expected_serosurvey$age_max)),
        y = c(min(expected_serosurvey$seroprev_lower), max(expected_serosurvey$seroprev_upper))
      )
    ),
    labels = list(
      x = "Age",
      y = "Seroprevalence",
      ymin = "seroprev_lower",
      ymax = "seroprev_upper",
      size = "n_sample"
    )
  )

  expect_equal(expected_plot, actual_plot)
})



test_that("plot_serosurvey creates a binned ggplot with correct structure", {
  data(veev2012)

  actual_serosurvey <- veev2012

  actual_plot <- extract_plot_data(plot_serosurvey(actual_serosurvey, bin_serosurvey = TRUE, bin_step = 10))

  expected_serosurvey <- create_expected_serosurvey(actual_serosurvey)

  expected_plot <- list(
    classes = c("gg", "ggplot"),
    layers = list(
      list(
        mapping = list(
          ymin = "seroprev_lower",
          ymax = "seroprev_upper"
        )
      ),
      list(
        mapping = list(
          y = "seroprev",
          size = "n_sample"
        )
      )
    ),
    coordinates = list(
      limits = list(
        x = c(2, max(expected_serosurvey$age_max)),
        y = c(min(expected_serosurvey$seroprev_lower), max(expected_serosurvey$seroprev_upper))
      )
    ),
    labels = list(
      x = "Age",
      y = "Seroprevalence",
      ymin = "seroprev_lower",
      ymax = "seroprev_upper",
      size = "n_sample"
    )
  )

  expect_equal(expected_plot, actual_plot)
})
