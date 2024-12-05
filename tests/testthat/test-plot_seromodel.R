library(testthat)
library(ggplot2)
library(purrr)

# NOTE: If the variable `expected_plot` needs to be updated because of changes in the plots
# You can use dput(extract_plot_data(plot)) to obtain the current structure of the plot
skip_on_cran()

# Common data ----
data(veev2012)
serosurvey <- veev2012

# Fit models for testing
suppressWarnings(
  seromodel_constant <- fit_seromodel(
    serosurvey = serosurvey,
    iter = 20
  )
)

suppressWarnings(
  seromodel_constant_serorev <- fit_seromodel(
    serosurvey = serosurvey,
    is_seroreversion = TRUE,
    seroreversion_prior = sf_normal(0, 1e-4),
    iter = 200
  )
)

suppressWarnings(
  seromodel_age <- fit_seromodel(
    serosurvey = serosurvey,
    model_type = "age",
    foi_index = get_foi_index(serosurvey, group_size = 20, model_type = "age"),
    iter = 400
  )
)

suppressWarnings(
  seromodel_time <- fit_seromodel(
    serosurvey = serosurvey,
    model_type = "time",
    foi_index = get_foi_index(serosurvey, group_size = 20, model_type = "time"),
    iter = 400
  )
)

create_prepared_serosurvey <- function(actual_serosurvey) {
  prepared_serosurvey <- serofoi:::prepare_serosurvey_for_plotting(
    add_age_group_to_serosurvey(actual_serosurvey)
  )
  return(prepared_serosurvey)
}

# Test plot_serosurvey ----
test_that("plot_serosurvey creates a ggplot with correct structure", {
  skip_on_cran()

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
    coordinates = list(limits = list(
      x = c(min(serosurvey$age_min), max(serosurvey$age_max)),
      y = c(min(prepared_serosurvey$seroprev_lower), max(prepared_serosurvey$seroprev_upper))
    )),
    labels = list(
      x = "Age", y = "Seroprevalence", ymin = "seroprev_lower",
      ymax = "seroprev_upper", size = "n_sample"
    )
  )

  expect_lists_equal_with_tolerance(expected_plot, actual_plot)
})

test_that("plot_serosurvey creates a binned ggplot with correct structure", {
  bin_step <- 10

  expect_warning(
    actual_plot <- extract_plot_data(
      plot_serosurvey(serosurvey, bin_serosurvey = TRUE, bin_step = bin_step)),
    "The last age interval will be truncated"
  )

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

# Test plot_summary ----
test_that("plot_summary creates a ggplot with correct structure", {
  seromodel <- seromodel_constant_serorev

  loo_estimate_digits <- 1
  central_estimate_digits <- 2
  seroreversion_digits <- 2
  rhat_digits <- 2
  size_text <- 11

  expect_warning(
    plot <- plot_summary(
      seromodel = seromodel,
      serosurvey = serosurvey,
      loo_estimate_digits = loo_estimate_digits,
      central_estimate_digits = central_estimate_digits,
      rhat_digits = rhat_digits,
      size_text = size_text
    ),
    "Some Pareto k diagnostic values are too high"
  )

  actual_data <- plot$data

  actual_plot <- extract_plot_data(plot)

  expect_warning(
    summary <- summarise_seromodel(
      seromodel = seromodel,
      serosurvey = serosurvey,
      loo_estimate_digits = loo_estimate_digits,
      central_estimate_digits = central_estimate_digits,
      rhat_digits = rhat_digits
    ),
    "Some Pareto k diagnostic values are too high"
  )

  expected_data <- data.frame(
    row = c(7, 6, 5, 4, 3, 2, 1),
    text = c(
      paste0("model_name: ", summary$model_name),
      paste0("elpd_loo: ", summary$elpd_loo),
      paste0("foi: ", summary$foi),
      paste0("foi_rhat: ", summary$foi_rhat),
      paste0("seroreversion_rate: ", summary$seroreversion_rate),
      paste0("seroreversion_rate_rhat: ", summary$seroreversion_rate_rhat),
      paste0("converged: ", summary$converged)
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
  seromodel <- seromodel_constant

  expect_warning(
    plot <- plot_seromodel(
      seromodel = seromodel,
      serosurvey = serosurvey
    ),
    "Some Pareto k diagnostic values are too high"
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




test_that("plot_seromodel with age foi creates a ggplot with correct structure", {
  old_opts <- options()
  seromodel <- seromodel_age

  age_foi_df <- data.frame(
    age = seq(1, 60, 1),
    foi = extract_central_estimates(
      seromodel,
      serosurvey,
      par_name = "foi_expanded"
    ) |> dplyr::pull(median)
  )

  expect_warning(
    plot <- plot_seromodel(
      seromodel = seromodel,
      serosurvey = serosurvey,
      foi_df = age_foi_df
    ),
    "Some Pareto k diagnostic values are too high"
  )

  actual_plot <- extract_plot_data(plot)

  expected_plot <- list(classes = c("gg", "ggplot"), layers = list(
    list(
      mapping = NULL,
      geom_params = list(
        xmin = 0, xmax = 1, ymin = 0.75, ymax = 1,
        scale = 1, clip = "inherit", halign = 0.5, valign = 0.5
      )
    ),
    list(mapping = NULL, geom_params = list(
      xmin = 0, xmax = 1,
      ymin = 0.5, ymax = 0.75, scale = 1, clip = "inherit",
      halign = 0.5, valign = 0.5
    )), list(mapping = NULL, geom_params = list(
      xmin = 0, xmax = 1, ymin = 0.25, ymax = 0.5, scale = 1,
      clip = "inherit", halign = 0.5, valign = 0.5
    )), list(
      mapping = NULL, geom_params = list(
        xmin = 0, xmax = 1,
        ymin = 0, ymax = 0.25, scale = 1, clip = "inherit",
        halign = 0.5, valign = 0.5
      )
    )
  ), coordinates = list(
    limits = list(x = c(0, 1), y = c(0, 1))
  ), labels = list())

  expect_lists_equal_with_tolerance(expected_plot, actual_plot)

  options(old_opts)
})





test_that("plot_seromodel with time foi creates a ggplot with correct structure", {
  seromodel <- seromodel_time

  time_foi_df <- data.frame(
    year = seq(1952, 2011, 1),
    foi = extract_central_estimates(
      seromodel,
      serosurvey,
      par_name = "foi_expanded"
    ) |> dplyr::pull(median)
  )

  expect_warning(
    plot <- plot_seromodel(
      seromodel = seromodel,
      serosurvey = serosurvey,
      foi_df = time_foi_df
    ),
    "Some Pareto k diagnostic values are too high"
  )

  actual_plot <- extract_plot_data(plot)

  expected_plot <- list(classes = c("gg", "ggplot"), layers = list(
    list(
      mapping = NULL,
      geom_params = list(
        xmin = 0, xmax = 1, ymin = 0.75, ymax = 1,
        scale = 1, clip = "inherit", halign = 0.5, valign = 0.5
      )
    ),
    list(mapping = NULL, geom_params = list(
      xmin = 0, xmax = 1,
      ymin = 0.5, ymax = 0.75, scale = 1, clip = "inherit",
      halign = 0.5, valign = 0.5
    )), list(mapping = NULL, geom_params = list(
      xmin = 0, xmax = 1, ymin = 0.25, ymax = 0.5, scale = 1,
      clip = "inherit", halign = 0.5, valign = 0.5
    )), list(
      mapping = NULL, geom_params = list(
        xmin = 0, xmax = 1,
        ymin = 0, ymax = 0.25, scale = 1, clip = "inherit",
        halign = 0.5, valign = 0.5
      )
    )
  ), coordinates = list(
    limits = list(x = c(0, 1), y = c(0, 1))
  ), labels = list())

  expect_lists_equal_with_tolerance(expected_plot, actual_plot)
})





# Test plot_seroprevalence_estimates ----
test_that("plot_seroprevalence_estimates creates a ggplot with correct structure", {
  seromodel <- seromodel_constant

  plot <- plot_seroprevalence_estimates(
    seromodel = seromodel,
    serosurvey = serosurvey
  )
  actual_plot <- extract_plot_data(plot)

  expected_plot <- list(
    classes = c("gg", "ggplot"), layers = list(
      list(mapping = list(
        ymin = "seroprev_lower", ymax = "seroprev_upper"
      ), geom_params = list(
        na.rm = FALSE, orientation = NA, width = 0.1
      )), list(mapping = list(
        y = "seroprev", size = "n_sample"
      ), geom_params = list(na.rm = FALSE)),
      list(mapping = list(x = "age", y = "median"), geom_params = list(
        na.rm = FALSE, orientation = NA
      )), list(mapping = list(
        x = "age", ymin = "lower", ymax = "upper"
      ), geom_params = list(
        na.rm = FALSE, orientation = NA, outline.type = "both"
      ))
    ),
    coordinates = list(limits = list(x = c(0, max(serosurvey$age_max)), y = NULL)),
    labels = list(
      x = "Age", y = "Seroprevalence", ymin = "seroprev_lower",
      ymax = "seroprev_upper", size = "n_sample"
    )
  )

  expect_lists_equal_with_tolerance(expected_plot, actual_plot)
})

# Test plot_foi_estimates ----
test_that("plot_foi_estimates creates a ggplot with correct structure", {
  seromodel <- seromodel_age
  plot <- plot_foi_estimates(
    seromodel = seromodel,
    serosurvey = serosurvey,
    foi_max = NULL,
    alpha = 0.05
  )

  foi_central_estimates <- extract_central_estimates(
    seromodel = seromodel,
    serosurvey = serosurvey,
    alpha = 0.05,
    par_name = "foi_expanded"
  )
  foi_max <- max(foi_central_estimates$upper)

  actual_plot <- extract_plot_data(plot)

  expected_plot <- list(
    classes = c("gg", "ggplot"), layers = list(list(mapping = list(
      ymin = "lower", ymax = "upper"
    ), geom_params = list(
      na.rm = FALSE,
      orientation = NA, outline.type = "both"
    )), list(mapping = list(
      y = "median"
    ), geom_params = list(na.rm = FALSE, orientation = NA))),
    coordinates = list(limits = list(x = NULL, y = c(0, foi_max))), labels = list(
      x = "Age", y = "Force-of-Infection", ymin = "lower",
      ymax = "upper"
    )
  )

  expect_lists_equal_with_tolerance(expected_plot, actual_plot)
})

# Test plot_rhats ----
test_that("plot_rhats creates a ggplot with correct structure", {
  seromodel <- seromodel_age

  plot <- plot_rhats(
    seromodel = seromodel,
    serosurvey = serosurvey,
    par_name = "foi_expanded"
  )
  actual_plot <- extract_plot_data(plot)

  rhats <- bayesplot::rhat(seromodel, "foi_expanded")


  expected_plot <- list(classes = c("gg", "ggplot"), layers = list(
    list(mapping = list(
      yintercept = "yintercept"
    ), geom_params = list(na.rm = FALSE)),
    list(mapping = list(y = "rhat"), geom_params = list(
      na.rm = FALSE,
      orientation = NA
    )), list(
      mapping = list(y = "rhat"),
      geom_params = list(na.rm = FALSE)
    )
  ), coordinates = list(
    limits = list(
      x = NULL,
      y = c(
        min(1.0, min(rhats)),
        max(1.02, max(rhats))
      )
    )
  ), labels = list(
    y = "Convergence (r-hats)", x = "Age", yintercept = "yintercept"
  ))

  expect_lists_equal_with_tolerance(expected_plot, actual_plot)
})
