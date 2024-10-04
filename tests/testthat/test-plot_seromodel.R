library(testthat)
library(ggplot2)
library(rlang)
library(purrr)

# NOTE: If the variable `expected_plot` needs to be updated because of changes in the plots
# You can run dput(actual_plot) to obtain the current structure of the plot

# Common data
data(veev2012)
serosurvey <- veev2012

seromodel <- fit_seromodel(
  serosurvey = serosurvey
)

create_prepared_serosurvey <- function(actual_serosurvey) {
  return(prepare_serosurvey_for_plotting(actual_serosurvey %>%
    add_age_group_to_serosurvey()))
}


test_that("plot_serosurvey creates a ggplot with correct structure", {
  actual_plot <- extract_plot_data(plot_serosurvey(serosurvey))

  prepared_serosurvey <- create_prepared_serosurvey(serosurvey)

  expected_plot <- list(
    classes = c("gg", "ggplot"),
    layers = list(
      list(
        mapping = list(
          ymin = "seroprev_lower", ymax = "seroprev_upper"
        ), geom_params = list(
          na.rm = FALSE, orientation = NA, width = 0.1
        )
      ),
      list(
        mapping = list(
          y = "seroprev", size = "n_sample"
        ), geom_params = list(na.rm = FALSE)
      )
    ),
    coordinates = list(limits = list(x = c(1, 60), y = c(
      0.3,
      1
    ))),
    labels = list(
      x = "Age", y = "Seroprevalence", ymin = "seroprev_lower",
      ymax = "seroprev_upper", size = "n_sample"
    )
  )

  expect_lists_equal_with_tolerance(expected_plot, actual_plot)
})



test_that("plot_serosurvey creates a binned ggplot with correct structure", {
  actual_plot <- extract_plot_data(plot_serosurvey(serosurvey, bin_serosurvey = TRUE, bin_step = 10))

  prepared_serosurvey <- create_prepared_serosurvey(serosurvey)

  expected_plot <- list(
    classes = c("gg", "ggplot"),
    layers = list(
      list(
        mapping = list(
          ymin = "seroprev_lower", ymax = "seroprev_upper"
        ),
        geom_params = list(
          na.rm = FALSE, orientation = NA, width = 0.1
        )
      ),
      list(
        mapping = list(
          y = "seroprev", size = "n_sample"
        ),
        geom_params = list(na.rm = FALSE)
      )
    ),
    coordinates = list(
      limits = list(x = c(2, 60), y = c(
        0.3,
        1
      ))
    ),
    labels = list(
      x = "Age", y = "Seroprevalence", ymin = "seroprev_lower",
      ymax = "seroprev_upper", size = "n_sample"
    )
  )

  expect_lists_equal_with_tolerance(expected_plot, actual_plot)
})

test_that("plot_summary creates a ggplot with correct structure", {
  loo_estimate_digits <- 1
  central_estimate_digits <- 2
  seroreversion_digits <- 2
  rhat_digits <- 2
  size_text <- 11

  plot <- plot_summary(
    seromodel = seromodel,
    serosurvey = serosurvey,
    loo_estimate_digits = loo_estimate_digits,
    central_estimate_digits = central_estimate_digits,
    rhat_digits = rhat_digits,
    size_text = size_text
  )
  actual_data <- plot$data

  actual_plot <- extract_plot_data(plot)

  expected_data <- data.frame(
    row = c(5, 4, 3, 2, 1),
    text = c(
      "model_name: constant_no_seroreversion",
      "elpd_loo: -23.9(se=12.2)",
      "foi: 0.082(95% CI, 0.063-0.1)",
      "foi_rhat: 1",
      "converged: yes"
    )
  )

  expected_plot <- list(
    classes = c("gg", "ggplot"),
    layers = list(list(
      mapping = list(
        label = "text"
      ),
      geom_params = list(
        parse = FALSE, check_overlap = FALSE,
        size.unit = "mm", na.rm = FALSE
      )
    )),
    coordinates = list(
      limits = list(
        x = NULL, y = NULL
      )
    ),
    labels = list(x = "x", y = "row", label = "text")
  )

  # Checks that the actual plot data are at most 20% disimilar to the expected data
  # using Levenshtein distance (adist)
  expect_true({
    all(map2_lgl(expected_data$text, actual_data$text, \(x, y) adist(x, y) / nchar(x) < 0.2))
  })
  expect_lists_equal_with_tolerance(expected_plot, actual_plot)
})


test_that("plot_seromodel creates a ggplot with correct structure", {
  plot <- plot_seromodel(
    seromodel = seromodel,
    serosurvey = serosurvey
  )
  actual_plot <- extract_plot_data(plot)

  expected_plot <- list(
    classes = c("gg", "ggplot"),
    layers = list(
      list(
        mapping = NULL,
        geom_params = list(
          xmin = 0, xmax = 1, ymin = 0.5, ymax = 1,
          scale = 1, clip = "inherit", halign = 0.5, valign = 0.5
        )
      ),
      list(mapping = NULL, geom_params = list(
        xmin = 0, xmax = 1,
        ymin = 0, ymax = 0.5, scale = 1, clip = "inherit", halign = 0.5,
        valign = 0.5
      ))
    ),
    coordinates = list(
      limits = list(x = c(
        0,
        1
      ), y = c(0, 1))
    ),
    labels = list()
  )

  expect_lists_equal_with_tolerance(expected_plot, actual_plot)
})
