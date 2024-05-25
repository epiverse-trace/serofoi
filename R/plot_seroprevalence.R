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

plot_seroprevalence_estimates <- function(
  seromodel,
  serosurvey,
  alpha = 0.05,
  ...
) {
  checkmate::assert_class(seromodel, "stanfit", null.ok = TRUE)

  seroprevalence_samples <- rstan::extract(
    seromodel, "prob_infected_expanded"
    )[[1]]
  seroprevalence_estimates <- data.frame(
    age = seq(1, max(serosurvey$age_max)),
    median = apply(seroprevalence_samples, 2, quantile, 0.5),
    lower = apply(seroprevalence_samples, 2, quantile, alpha),
    upper = apply(seroprevalence_samples, 2, quantile, 1 - alpha)
  )

  seroprevalence_plot <- plot_serosurvey(
    serosurvey = serosurvey,
    ...
    ) +
    ggplot2::geom_line(
      data = seroprevalence_estimates,
      ggplot2::aes(x = age, y = median),
      colour = "#7a0177"
    ) +
    ggplot2::geom_ribbon(
      data = seroprevalence_estimates,
      ggplot2::aes(x = age, ymin = lower, ymax = upper),
      fill = "#c994c7", alpha = 0.5
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0, max(serosurvey$age_max))
    )

  return(seroprevalence_plot)
}
