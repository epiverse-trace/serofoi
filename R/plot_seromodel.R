#' Prepares serosurvey for plotting
#'
#' Adds seroprevalence values with corresponding binomial confidence interval
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
prepare_serosurvey_for_plotting <- function( #nolint
  serosurvey,
  alpha = 0.05
  ) {

  serosurvey <- serosurvey %>%
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
    dplyr::relocate(.data$age_group)
}

#' Construct age-group variable from age column
#'
#' Generates age intervals of length step in the interval spanned by
#' `age_min` and `age_max` in a serosurvey.
#' In cases where `max(age_max)%%(step+1)!=0`, the last age interval is
#' truncated and will have a different length than the others.
#' @inheritParams plot_serosurvey
#' @param  step step used to split the age interval
#' @return Serosurvey with addition factor variable grouping `age_intervals`.
#'  The interval is taken as closed to the right and to the left.
get_age_intervals <- function(serosurvey, step) {
  age_min <- min(serosurvey$age_min)
  age_max <- max(serosurvey$age_max)

  checkmate::assert_int(age_min, lower = 0)
  checkmate::assert_int(age_max, lower = age_min)
  checkmate::assert_int(step, lower = 2, upper = age_max)

  limits_low <- as.integer(
    seq(
      age_min, age_max,
      by = step
    )
  )

  if ((age_max - age_min) %% step != 0) {
    warn_msg <- "(age_min - age_max) is not an integer multiple of step.
    The last age interval will be truncated to "
    warn_msg <- paste0(
      warn_msg, "[", limits_low[length(limits_low)], ",", age_max, "]"
    )
    warning(warn_msg)
  }

  # prepare breaks
  lim_breaks <- c(limits_low, age_max)

  # define age groups closed to the left and closed to the right
  survey_features <- data.frame(
    age_min = limits_low,
    age_max = limits_low + step - 1
  ) %>% add_age_bins()

  serosurvey$age_interval <- cut(
      x = serosurvey$age_group,
      breaks = lim_breaks,
      include.lowest = TRUE, right = FALSE,
      labels = survey_features$group
    )

  return(serosurvey)
}

#' Plots seroprevalence from the given serosurvey
#'
#' @inheritParams fit_seromodel
#' @param size_text Size of text for plotting (`base_size` in
#' [ggplot2][ggplot2::theme_bw])
#' @param bin_serosurvey If `TRUE`, `serodata` is binned by means of
#'   `prepare_bin_serosurvey`. Otherwise, age groups are kept as originally input.
#' @param bin_step Integer specifying the age groups bin size to be used when
#' `bin_serosurvey` is set to `TRUE`.
#' @return ggplot object with seroprevalence plot
#' @export
plot_serosurvey <- function(
    serosurvey,
    size_text = 11,
    bin_serosurvey = FALSE,
    bin_step = 5
    ) {
  serosurvey <- validate_serosurvey(serosurvey = serosurvey) %>%
    add_age_group_to_serosurvey()

  if (bin_serosurvey) {
    age_min <- min(serosurvey$age_min)
    age_max <- max(serosurvey$age_max)
    checkmate::assert_int(bin_step, lower = 2, upper = age_max)

    serosurvey <- get_age_intervals(
      serosurvey = serosurvey,
      step = bin_step
    )

    serosurvey <- serosurvey %>%
      dplyr::group_by(.data$age_interval) %>%
      dplyr::summarise(
        sample_size = sum(.data$sample_size),
        n_seropositive = sum(.data$n_seropositive)
      ) %>%
      dplyr::mutate(
        survey_year = unique(serosurvey$survey_year),
        age_min = as.integer(gsub("[[]|\\,.*", "\\1", .data$age_interval)) + 1,
        age_max = as.integer(gsub(".*\\,|[]]", "\\1", .data$age_interval))
      ) %>%
      add_age_group_to_serosurvey()
  }

  serosurvey <- prepare_serosurvey_for_plotting(serosurvey)

  min_prev <- min(serosurvey$seroprev_lower)
  max_prev <- max(serosurvey$seroprev_upper)

  seroprev_plot <- ggplot2::ggplot(
    data = serosurvey,
    ggplot2::aes(x = .data$age_group)
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = .data$seroprev_lower,
        ymax = .data$seroprev_upper
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
#' @returns ggplot object with seroprevalence estimates and serosurveys plots
#' @export
plot_seroprevalence_estimates <- function(
  seromodel,
  serosurvey,
  alpha = 0.05,
  size_text = 11,
  bin_serosurvey = FALSE,
  bin_step = 5
) {
  checkmate::assert_class(seromodel, "stanfit", null.ok = TRUE)

  seroprevalence_central_estimates <- data.frame( #nolint
    median = 0.0,
    lower = 0.0,
    upper = 0.0,
    age = 0
  ) %>%
  rbind(
    extract_central_estimates(
      seromodel = seromodel,
      serosurvey = serosurvey,
      alpha = alpha,
      par_name = "prob_infected_expanded"
  ) %>%
    mutate(age = seq(1, max(serosurvey$age_max)))
  )

  seroprevalence_plot <- plot_serosurvey(
    serosurvey = serosurvey,
    size_text = size_text,
    bin_serosurvey = bin_serosurvey,
    bin_step = bin_step
    ) +
    ggplot2::geom_line(
      data = seroprevalence_central_estimates,
      ggplot2::aes(x = .data$age, y = median),
      colour = "#7a0177"
    ) +
    ggplot2::geom_ribbon(
      data = seroprevalence_central_estimates,
      ggplot2::aes(x = .data$age, ymin = .data$lower, ymax = .data$upper),
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
  foi_max = NULL,
  size_text = 11
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
    if (!is.null(foi_df)) {
      foi_central_estimates <- foi_central_estimates %>%
        left_join(foi_df, by = "age")
    }
    foi_plot <- ggplot2::ggplot(
      data = foi_central_estimates, ggplot2::aes(x = .data$age)
    )
  } else if (startsWith(model_name, "time")) {
    checkmate::assert_names(names(serosurvey), must.include = "survey_year")
    xlab <- "Year"
    ages <- rev(1:max(serosurvey$age_max))
    years <- unique(serosurvey$survey_year) - ages
    foi_central_estimates <- mutate(
      foi_central_estimates,
      year = years
    )
    if (!is.null(foi_df)) {
      foi_central_estimates <- foi_central_estimates %>%
        left_join(foi_df, by = "year")
    }
    foi_plot <- ggplot2::ggplot(
      data = foi_central_estimates, ggplot2::aes(x = .data$year)
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
        ggplot2::aes(y = .data$foi),
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
      data = rhats_df, ggplot2::aes(x = .data$age)
    )
  } else if (startsWith(model_name, "time")) {
    checkmate::assert_names(names(serosurvey), must.include = "survey_year")
    xlab <- "Year"
    ages <- rev(1:max(serosurvey$age_max))
    years <- unique(serosurvey$survey_year) - ages
    rhats_df <- data.frame(
      year = years,
      rhat = rhats
    )

    rhats_plot <- ggplot2::ggplot(
      data = rhats_df, ggplot2::aes(x = .data$year)
    )
  }

  rhats_plot <- rhats_plot +
    ggplot2::geom_hline(
      yintercept = 1.01,
      linetype = "dashed"
    ) +
    ggplot2::geom_line(ggplot2::aes(y = .data$rhat)) +
    ggplot2::geom_point(ggplot2::aes(y = .data$rhat)) +
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

#' Plots model summary
#'
#' @inheritParams summarise_seromodel
#' @inheritParams plot_serosurvey
#' @return ggplot object with a summary of the specified model
#' @export
plot_summary <- function(
  seromodel,
  serosurvey,
  loo_estimate_digits = 1,
  central_estimate_digits = 2,
  rhat_digits = 2,
  size_text = 11
) {
  checkmate::assert_class(seromodel, "stanfit", null.ok = TRUE)

  summary_table <- summarise_seromodel(
    seromodel = seromodel,
    serosurvey = serosurvey,
    loo_estimate_digits = loo_estimate_digits,
    central_estimate_digits = central_estimate_digits,
    rhat_digits = rhat_digits
    ) %>%
    t() #convert summary to table

  summary_df <- data.frame(
    row = rev(seq_len(NCOL(summary_table))),
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
      fontface = "bold",
      size = size_text / 2.5
    )

  return(summary_plot)
}

#' Visualise results of the provided model
#'
#' @inheritParams plot_summary
#' @inheritParams plot_seroprevalence_estimates
#' @inheritParams plot_foi_estimates
#' @inheritParams plot_rhats
#' @return seromodel summary plot
#' @export
plot_seromodel <- function(
  seromodel,
  serosurvey,
  alpha = 0.05,
  bin_serosurvey = FALSE,
  bin_step = 5,
  foi_df = NULL,
  foi_max = NULL,
  loo_estimate_digits = 1,
  central_estimate_digits = 2,
  seroreversion_digits = 2,
  rhat_digits = 2,
  size_text = 11
) {
  checkmate::assert_class(seromodel, "stanfit", null.ok = TRUE)

  summary_plot <- plot_summary(
    seromodel = seromodel,
    serosurvey = serosurvey,
    loo_estimate_digits = loo_estimate_digits,
    central_estimate_digits = central_estimate_digits,
    rhat_digits = rhat_digits,
    size_text = size_text
  )

  seroprev_plot <- plot_seroprevalence_estimates(
    seromodel = seromodel,
    serosurvey = serosurvey,
    alpha = alpha,
    size_text = size_text,
    bin_serosurvey = bin_serosurvey,
    bin_step = bin_step
  )

  plot_list <- list(
    summary_plot,
    seroprev_plot
  )

  model_name <- seromodel@model_name
  if (!startsWith(model_name, "constant")) {
    foi_plot <- plot_foi_estimates(
      seromodel,
      serosurvey,
      alpha = alpha,
      foi_df = foi_df,
      foi_max = foi_max,
      size_text
    )

    rhats_plot <- plot_rhats(
      seromodel = seromodel,
      serosurvey = serosurvey,
      size_text = size_text
    )

    plot_list <- c(
      plot_list,
      list(foi_plot, rhats_plot)
    )
  }

  seromodel_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = 1)
  return(seromodel_plot)
}
