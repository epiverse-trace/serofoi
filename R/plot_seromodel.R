prepare_serosurvey_for_plotting <- function(
  serosurvey,
  alpha = 0.05
  ) {

  serosurvey <- serosurvey %>%
    add_age_group_to_serosurvey() %>%
    cbind(
      Hmisc::binconf(
        serosurvey$n_seropositive,
        serosurvey$sample_size,
        alpha = alpha,
        method = "exact",
        return.df = TRUE
      )
    ) %>%
    dplyr::rename(
      seroprev = "PointEst",
      seroprev_lower = "Lower",
      seroprev_upper = "Upper"
    ) %>%
    dplyr::arrange(.data$age_group) %>%
    dplyr::relocate(age_group)
}

plot_serosurvey <- function(
    serosurvey,
    size_text = 11
    ) {
  serosurvey <- validate_serosurvey(serosurvey = serosurvey) %>%
    prepare_serosurvey_for_plotting()

  min_prev <- min(serosurvey$seroprev_lower)
  max_prev <- max(serosurvey$seroprev_upper)

  seroprev_plot <- ggplot2::ggplot(
    data = serosurvey,
    ggplot2::aes(x = .data$age_group)
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = .data$seroprev_lower,
        ymax = .data$seroprev_upper,
      ),
      width = 0.1
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        y = .data$seroprev,
        size = .data$sample_size
      ),
      fill = "#7a0177", colour = "black", shape = 21
    ) +
    ggplot2::coord_cartesian(
      xlim = c(min(serosurvey$age_min), max(serosurvey$age_max)),
      ylim = c(min_prev, max_prev)
    ) +
    ggplot2::theme_bw(size_text) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::ylab("Seroprevalence") +
    ggplot2::xlab("Age")

  return(seroprev_plot)
}

extract_central_estimates <- function(
  seromodel,
  serosurvey,
  alpha = 0.05,
  par_name = "foi_vector"
) {
  samples <- rstan::extract(seromodel, par_name)[[1]]
  central_estimates <- data.frame(
    median = apply(samples, 2, quantile, 0.5),
    lower = apply(samples, 2, quantile, alpha),
    upper = apply(samples, 2, quantile, 1 - alpha)
  )

  return(central_estimates)
}

plot_seroprevalence_estimates <- function(
  seromodel,
  serosurvey,
  alpha = 0.05,
  ...
) {
  checkmate::assert_class(seromodel, "stanfit", null.ok = TRUE)

  seroprevalence_central_estimates <- extract_central_estimates(
    seromodel = seromodel,
    serosurvey = serosurvey,
    alpha = alpha,
    par_name = "prob_infected_expanded"
  ) %>%
    mutate(age = seq(1, max(serosurvey$age_max)))

  seroprevalence_plot <- plot_serosurvey(
    serosurvey = serosurvey,
    ...
    ) +
    ggplot2::geom_line(
      data = seroprevalence_central_estimates,
      ggplot2::aes(x = age, y = median),
      colour = "#7a0177"
    ) +
    ggplot2::geom_ribbon(
      data = seroprevalence_central_estimates,
      ggplot2::aes(x = age, ymin = lower, ymax = upper),
      fill = "#c994c7", alpha = 0.5
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0, max(serosurvey$age_max))
    )

  return(seroprevalence_plot)
}

plot_foi_estimates <- function(
  seromodel,
  serosurvey,
  alpha = 0.05,
  foi_df = NULL,
  size_text = 11,
  foi_max = NULL
) {
  # TODO: Add checks for foi_df (size, colnames, etc.)
  checkmate::assert_class(seromodel, "stanfit", null.ok = TRUE)

  foi_central_estimates <- extract_central_estimates(
    seromodel = seromodel,
    serosurvey = serosurvey,
    alpha = alpha,
    par_name = "foi_expanded"
  )

  if (is.null(foi_max))
    foi_max <- max(foi_central_estimates$upper)

  if (startsWith(seromodel@model_name, "age")) {
    xlab <- "Age"
    ages <- 1:max(serosurvey$age_max)
    foi_central_estimates <- mutate(
      foi_central_estimates,
      age = ages
    )
    if(!is.null(foi_df)) {
      foi_central_estimates <- foi_central_estimates %>%
        left_join(foi_df, by = "age")
    }
    foi_plot <- ggplot2::ggplot(
      data = foi_central_estimates, ggplot2::aes(x = age)
    )
  } else if (startsWith(seromodel@model_name, "time")) {
    xlab <- "Year"
    ages <- rev(1:max(serosurvey$age_max))
    years <- unique(serosurvey$tsur) - ages
    foi_central_estimates <- mutate(
      foi_central_estimates,
      year = years
    )
    if(!is.null(foi_df)) {
      foi_central_estimates <- foi_central_estimates %>%
        left_join(foi_df, by = "year")
    }
    foi_plot <- ggplot2::ggplot(
      data = foi_central_estimates, ggplot2::aes(x = year)
    )
  }

  foi_plot <- foi_plot +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = .data$lower,
        ymax = .data$upper
      ),
      fill = "#41b6c4",
      alpha = 0.5
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$median),
      colour = "#253494"
    ) +
    ggplot2::theme_bw(size_text) +
    ggplot2::coord_cartesian(ylim = c(0, foi_max)) +
    ggplot2::ylab("Force-of-Infection") +
    ggplot2::xlab(xlab)

  if (!is.null(foi_df)) {
    foi_plot <- foi_plot +
      ggplot2::geom_line(
        ggplot2::aes(y = foi),
        colour = "#b30909"
      )
  }

  return(foi_plot)
}
