#' Prepares serosurvey for plotting
#'
#' Adds seroprevalence values with corresponing binomial confidence interval
#' @inheritParams fit_seromodel
#' @param alpha 1 - alpha indicates the confidence level to be used
#' @return serosurvey with additional columns:
#' \describe{
#'  \item{seroprev}{Seroprevalence computed as the proportion of positive
#'                  cases `n_seropositive` in the number of samples
#'                  `sample_size` for each age group}
#'  \item{seroprev_lower}{Lower limit of the binomial confidence interval
#'                        of `seroprev`}
#'  \item{seroprev_upper}{Upper limit of the binomial confidence interval
#'                        of `seroprev`}
#' }
#' @export
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

#' Plots seroprevalence from the given serosurvey
#'
#' @inheritParams fit_seromodel
#' @param size_text Size of text for plotting (`base_size` in
#' [ggplot2][ggplot2::theme_bw])
#' @return ggplot object with seroprevalence plot
#' @export
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

#' Extracts central estimates from stan_fit object for specified parameter
#'
#' @param seromodel stan_fit object obtained from sampling a model
#' with [fit_seromode]
#' @inheritParams fit_seromodel
#' @param alpha 1 - alpha indicates the credibility level to be used
#' @param par_name String specifying the parameter to be extracted
#' from `seromodel`
#' @returns A dataframe with the following columns
#' \describe{
#'  \item{`median`}{Median of the samples computed as the 0.5 quantile}
#'  \item{`lower`}{Lower quantile `alpha`}
#'  \item{`upper`}{Upper quantile `1 - alpha`}
#' }
#' @export
extract_central_estimates <- function(
  seromodel,
  serosurvey,
  alpha = 0.05,
  par_name = "foi_vector"
) {
  samples <- rstan::extract(seromodel, par_name)[[1]] %>%
    as.matrix() #to deal with 1-time estimates
  central_estimates <- data.frame(
    median = apply(samples, 2, quantile, 0.5),
    lower = apply(samples, 2, quantile, alpha),
    upper = apply(samples, 2, quantile, 1 - alpha)
  )

  return(central_estimates)
}

#' Plot seroprevalence estimates on top of the serosurvey
#'
#' @inheritParams extract_central_estimates
#' @inheritParams plot_serosurvey
#' @returns ggplot object with seroprevalence estimates and serisurvey plots
#' @export
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

#' Plots force-of-infection central estimates
#'
#' @inheritParams extract_central_estimates
#' @inheritParams fit_seromodel
#' @param foi_df Dataframe with columns
#' \describe{
#'  \item{`year`/`age`}{Year/Age (depending on the model)}
#'  \item{`foi`}{Force-of-infection values by year/age}
#' }
#' @inheritParams plot_serosurvey
#' @param foi_max Max force-of-infection value for plotting
#' @return ggplot object with estimated force-of-infection
#' @export
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

  model_name <- seromodel@model_name
  stopifnot(
    "seromodel@name should start with either 'age' or 'time'" =
    startsWith(model_name, "age") | startsWith(model_name, "time")
  )

  foi_central_estimates <- extract_central_estimates(
    seromodel = seromodel,
    serosurvey = serosurvey,
    alpha = alpha,
    par_name = "foi_expanded"
  )

  if (is.null(foi_max))
    foi_max <- max(foi_central_estimates$upper)

  if (startsWith(model_name, "age")) {
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
  } else if (startsWith(model_name, "time")) {
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

#' Plot r-hats convergence criteria for the specified model
#'
#' @inheritParams extract_central_estimates
#' @inheritParams plot_serosurvey
#' @return ggplot object showing the r-hats of the model to be compared with the
#' convergence criteria (horizontal dashed line)
plot_rhats <- function(
  seromodel,
  serosurvey,
  par_name = "foi_expanded",
  size_text = 11
) {
  checkmate::assert_class(seromodel, "stanfit", null.ok = TRUE)

  model_name <- seromodel@model_name
  stopifnot(
    "seromodel@name should start with either 'age' or 'time'" =
    startsWith(model_name, "age") | startsWith(model_name, "time")
  )

  rhats <- bayesplot::rhat(seromodel, par_name)

  if (startsWith(model_name, "age")) {
    xlab <- "Age"
    ages <- 1:max(serosurvey$age_max)
    rhats_df <- data.frame(
      age = ages,
      rhat = rhats
    )

    rhats_plot <- ggplot2::ggplot(
      data = rhats_df, ggplot2::aes(x = age)
    )
  } else if (startsWith(model_name, "time")) {
    xlab <- "Year"
    ages <- rev(1:max(serosurvey$age_max))
    years <- unique(serosurvey$tsur) - ages
    rhats_df <- data.frame(
      year = years,
      rhat = rhats
    )

    rhats_plot <- ggplot2::ggplot(
      data = rhats_df, ggplot2::aes(x = year)
    )
  }

  rhats_plot <- rhats_plot +
    ggplot2::geom_hline(
      yintercept = 1.01,
      linetype = 'dashed'
    ) +
    ggplot2::geom_line(ggplot2::aes(y = rhat)) +
    ggplot2::geom_point(ggplot2::aes(y = rhat)) +
    ggplot2::coord_cartesian(
      ylim = c(
        min(1.0, min(rhats_df$rhat)),
        max(1.02, max(rhats_df$rhat))
      )
    ) +
    ggplot2::theme_bw(size_text) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab("Convergence (r-hats)")

    return(rhats_plot)
}

#' Summarise specified model
#'
#' @inheritParams extract_central_estimates
#' @param elpd_digits Number of elpd digits to show
#' @param foi_digits Number of foi digits to show
#' @param seroreversion_digits Number of seroreversion rate digits to show
#' @return A list showing
#' \describe{
#'  \item{`model_name`}{Name of the model}
#'  \item{`elpd`}{elpd and its standard deviation}
#'  \item{`foi`}{Estimated foi with credible interval (for 'constant' model)}
#'  \item{`foi_rhat`}{foi rhat value (for 'constant' model)}
#'  \item{`seroreversion_rate`}{Estimated seroreversion rate}
#'  \item{`rhat_seroreversion_rate`}{Seroreversion rate rhat value}
#' }
#' @export
summarise_seromodel <- function(
  seromodel,
  serosurvey,
  alpha = 0.05,
  elpd_digits = 1,
  foi_digits = 2,
  seroreversion_digits = 2,
  rhat_digits = 2
) {
  checkmate::assert_class(seromodel, "stanfit", null.ok = TRUE)

  loo_fit <- loo::loo(
    seromodel,
    pars = c(parameter_name = "log_likelihood")
  )
  elp <- loo_fit$estimates["elpd_loo", ] %>% round(elpd_digits)

  summary_list <- list(
    model_name = seromodel@model_name,
    elpd = paste0(elpd[1], "(se=", elpd[2], ")")
  )

  model_name <- seromodel@model_name
  if (startsWith(model_name, "constant")) {
    foi_central_estimates <- extract_central_estimates(
      seromodel = seromodel,
      serosurvey = serosurvey,
      alpha = alpha,
      par_name = "foi"
    ) %>%
    signif(foi_digits)

    foi_summary <- paste0(
      foi_central_estimates$median,
      "(", 1 - alpha, "% CI, ",
      foi_central_estimates$lower, "-",
      foi_central_estimates$upper, ")"
    )
    rhat_foi <- bayesplot::rhat(seromodel, "foi") %>%
      signif(rhat_digits)
    summary_list <- append(
      summary_list,
      list(
        foi = foi_summary,
        rhat_foi = rhat_foi
      )
    )
  }

  if (!endsWith(model_name, "no_seroreversion")) {
    seroreversion_rate_central_estimates <- extract_central_estimates(
      seromodel = seromodel,
      serosurvey = serosurvey,
      alpha = alpha,
      par_name = "seroreversion_rate"
    ) %>%
      signif(seroreversion_digits)

    seroreversion_rate_summary <- paste0(
      seroreversion_rate_central_estimates$median,
      "(", 1 - alpha, "% CI, ",
      seroreversion_rate_central_estimates$lower, "-",
      seroreversion_rate_central_estimates$upper, ")"
    )
    rhat_seroreversion_rate <- bayesplot::rhat(
      seromodel,
      "seroreversion_rate"
      ) %>%
      signif(rhat_digits)
    summary_list <- append(
      summary_list,
      list(
        seroreversion_rate = seroreversion_rate_summary,
        rhat_seroreversion_rate = rhat_seroreversion_rate
      )
    )
  }
  return(summary_list)
}

#' Plots model summary
#'
#' @inheritParams summarise_seromodel
#' @return ggplot object with a summary of the specified model
#' @export
plot_summary <- function(
  seromodel,
  serosurvey,
  ...
) {
  checkmate::assert_class(seromodel, "stanfit", null.ok = TRUE)

  summary_table <- summarise_seromodel(
    seromodel = seromodel,
    serosurvey = serosurvey,
    ...
    ) %>%
    t() #convert summary to table

  summary_df <- data.frame(
    row = NCOL(summary_table):1,
    text = paste0(colnames(summary_table), ": ", summary_table[1, ])
  )

  summary_plot <- ggplot2::ggplot(
    summary_df,
    ggplot2::aes(x = 1, y = row)) +
    ggplot2::scale_y_continuous(
      limits = c(0, nrow(summary_df) + 1),
      breaks = NULL
    ) +
    ggplot2::theme_void() +
    ggplot2::geom_text(
      ggplot2::aes(label = text),
      fontface = "bold"
    )

  return(summary_plot)
}

