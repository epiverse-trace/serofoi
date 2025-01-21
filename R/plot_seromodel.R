#' Prepares serosurvey for plotting
#'
#' Adds seroprevalence values with corresponding binomial confidence interval
#' @inheritParams fit_seromodel
#' @param alpha 1 - alpha indicates the confidence level to be used
#' @return serosurvey with additional columns:
#' \describe{
#'  \item{seroprev}{Seroprevalence computed as the proportion of positive
#'                  cases `n_seropositive` in the number of samples
#'                  `n_sample` for each age group}
#'  \item{seroprev_lower}{Lower limit of the binomial confidence interval
#'                        of `seroprev`}
#'  \item{seroprev_upper}{Upper limit of the binomial confidence interval
#'                        of `seroprev`}
#' }
prepare_serosurvey_for_plotting <- function( #nolint
  serosurvey,
  alpha = 0.05
  ) {
  # The binomial confidence interval calculation is based on:
  # https://forum.posit.co/t/apply-binomial-test-for-each-row-in-a-data-table/32112/2 #nolint
  serosurvey$seroprev <- serosurvey$n_seropositive / serosurvey$n_sample
  serosurvey <- dplyr::mutate(
    serosurvey,
    seroprev = n_seropositive / n_sample,
    binconf = purrr::pmap(
      .l = serosurvey,
      .f = purrr::lift_vd(..f = function(dat)
      {
        ci <- stats::binom.test(
          x = dat["n_seropositive"],
          n = dat["n_sample"],
          p = dat["seroprev"],
          conf.level = 1 - alpha)$conf.int
        names(x = ci) <- c("seroprev_lower", "seroprev_upper")
        return(ci)
      })
    )
  ) |>
  tidyr::unnest_wider(.data$binconf) |>
    dplyr::arrange(.data$age_group) |>
    dplyr::relocate(!!dplyr::sym("age_group"))

  return(serosurvey)
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
  survey_features <- add_age_bins(
    data.frame(
      age_min = limits_low,
      age_max = limits_low + step - 1
    )
  )

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
#'   `prepare_bin_serosurvey`.
#'   Otherwise, age groups are kept as originally input.
#' @param bin_step Integer specifying the age groups bin size to be used when
#' `bin_serosurvey` is set to `TRUE`.
#' @return ggplot object with seroprevalence plot
#' @examples
#' # Chikungunya example serosurvey
#' data(chik2015)
#' plot_serosurvey(chik2015)
#'
#' # VEEV example serosurvey
#' data(veev2012)
#' plot_serosurvey(veev2012)
#'
#' # Chagas disease example serosurvey
#' data(chagas2012)
#' plot_serosurvey(chagas2012, bin_serosurvey = TRUE)
#' @export
plot_serosurvey <- function(
    serosurvey,
    size_text = 11,
    bin_serosurvey = FALSE,
    bin_step = 5
    ) {
  serosurvey <- add_age_group_to_serosurvey(validate_serosurvey(serosurvey))

  if (bin_serosurvey) {
    age_max <- max(serosurvey$age_max)
    checkmate::assert_int(bin_step, lower = 2, upper = age_max)

    serosurvey <- get_age_intervals(
      serosurvey = serosurvey,
      step = bin_step
    )

    serosurvey <- dplyr::group_by(serosurvey, .data$age_interval) |>
      dplyr::summarise(
        n_sample = sum(.data$n_sample),
        n_seropositive = sum(.data$n_seropositive)
      ) |>
      dplyr::mutate(
        age_min = as.integer(gsub("[[]|\\,.*", "\\1", .data$age_interval)) + 1,
        age_max = as.integer(gsub(".*\\,|[]]", "\\1", .data$age_interval))
      ) |>
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
        size = .data$n_sample
      ),
      fill = "#7a0177", colour = "black", shape = 21
    ) +
    ggplot2::coord_cartesian(
      xlim = c(min(serosurvey$age_min), max(serosurvey$age_max)),
      ylim = c(min_prev, max_prev),
      default = TRUE
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
#' @examples
#' data(veev2012)
#' seromodel <- fit_seromodel(veev2012, iter = 100)
#' central_estimates <- extract_central_estimates(
#'   seromodel,
#'   veev2012,
#'   par_name = "foi"
#' )
#' @export
extract_central_estimates <- function(
  seromodel,
  serosurvey,
  alpha = 0.05,
  par_name = "foi_vector"
) {
  samples <- as.matrix( #to deal with 1-time estimates
    rstan::extract(seromodel, par_name)[[1]]
  )
  central_estimates <- data.frame(
    median = apply(samples, 2, stats::quantile, 0.5),
    lower = apply(samples, 2, stats::quantile, alpha),
    upper = apply(samples, 2, stats::quantile, 1 - alpha)
  )

  return(central_estimates)
}

#' Plot seroprevalence estimates on top of the serosurvey
#'
#' @inheritParams extract_central_estimates
#' @inheritParams plot_serosurvey
#' @returns ggplot object with seroprevalence estimates and serosurveys plots
#' @examples
#' data(veev2012)
#' seromodel <- fit_seromodel(veev2012, iter = 100)
#' plot_seroprevalence_estimates(seromodel, veev2012)
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
  ) |>
  rbind(
    extract_central_estimates(
      seromodel = seromodel,
      serosurvey = serosurvey,
      alpha = alpha,
      par_name = "prob_infected_expanded"
  ) |>
    dplyr::mutate(age = seq(1, max(serosurvey$age_max)))
  )

  seroprevalence_plot <- plot_serosurvey(
    serosurvey = serosurvey,
    size_text = size_text,
    bin_serosurvey = bin_serosurvey,
    bin_step = bin_step
    ) +
    ggplot2::geom_line(
      data = seroprevalence_central_estimates,
      ggplot2::aes(x = .data$age, y = .data$median),
      colour = "#7a0177"
    ) +
    ggplot2::geom_ribbon(
      data = seroprevalence_central_estimates,
      ggplot2::aes(x = .data$age, ymin = .data$lower, ymax = .data$upper),
      fill = "#c994c7", alpha = 0.5
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0, max(serosurvey$age_max)),
      default = TRUE
    )

  return(seroprevalence_plot)
}

#' Plots force-of-infection central estimates
#'
#' @inheritParams extract_central_estimates
#' @inheritParams fit_seromodel
#' @inheritParams plot_rhats
#' @param foi_df Dataframe with columns
#' \describe{
#'  \item{`year`/`age`}{Year/Age (depending on the model)}
#'  \item{`foi`}{Force-of-infection values by year/age}
#' }
#' @inheritParams plot_serosurvey
#' @param foi_max Max force-of-infection value for plotting
#' @return ggplot object with estimated force-of-infection
#' @examples
#' data(chagas2012)
#' seromodel <- fit_seromodel(
#'   serosurvey = chagas2012,
#'   model_type = "time",
#'   foi_index = data.frame(
#'     year = 1935:2011,
#'     foi_index = c(rep(1, 46), rep(2, 31))
#'   ),
#'   iter = 100,
#'   chains = 2
#' )
#' plot_foi_estimates(seromodel, chagas2012)
#' @export
plot_foi_estimates <- function(
  seromodel,
  serosurvey,
  alpha = 0.05,
  foi_df = NULL,
  foi_max = NULL,
  size_text = 11,
  plot_constant = FALSE,
  x_axis = NA
) {
  # TODO: Add checks for foi_df (size, colnames, etc.)
  checkmate::assert_class(seromodel, "stanfit", null.ok = TRUE)

  model_name <- seromodel@model_name
  stopifnot(
    "seromodel@name should start with either 'age', 'time' or 'constant' " =
      startsWith(model_name, "age") | startsWith(model_name, "time") |
      startsWith(model_name, "constant")
  )

  error_msg_x_axis <- paste0(
    "x_axis specification as either 'age' or 'time' when plotting ",
    "constant FOI estimates is required to avoid ambiguity"
  )
  validate_plot_constant(
    plot_constant = plot_constant,
    x_axis = x_axis,
    model_name = model_name,
    error_msg_x_axis = error_msg_x_axis
  )

  foi_central_estimates <- extract_central_estimates(
    seromodel = seromodel,
    serosurvey = serosurvey,
    alpha = alpha,
    par_name = "foi_expanded"
  )

  if (is.null(foi_max))
    foi_max <- max(foi_central_estimates$upper)

  if (startsWith(model_name, "age") || (plot_constant && (x_axis == "age"))) {
    xlab <- "Age"
    ages <- 1:max(serosurvey$age_max)
    foi_central_estimates <- dplyr::mutate(
      foi_central_estimates,
      age = ages
    )

    if (!is.null(foi_df)) {
      foi_central_estimates <- dplyr::left_join(
        foi_central_estimates, foi_df,
        by = "age"
      )
    }

    foi_plot <- ggplot2::ggplot(
      data = foi_central_estimates, ggplot2::aes(x = .data$age)
    )

  } else if (
    startsWith(model_name, "time") ||
    (plot_constant && x_axis == "time")
    ) {
    checkmate::assert_names(names(serosurvey), must.include = "survey_year")
    xlab <- "Year"
    ages <- rev(1:max(serosurvey$age_max))
    years <- unique(serosurvey$survey_year) - ages
    foi_central_estimates <- dplyr::mutate(
      foi_central_estimates,
      year = years
    )
    if (!is.null(foi_df)) {
      foi_central_estimates <- dplyr::left_join(
        foi_central_estimates, foi_df,
        by = "year"
      )
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
    ggplot2::coord_cartesian(ylim = c(0, foi_max), default = TRUE) +
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
#' @inheritParams plot_serosurvey
#' @inheritParams plot_summary
#' @param x_axis either `"time"` or `"age"`. Specifies time axis values
#' label for constant model additional plots. Only relevant when
#'and `seromodel@model_name == "constant"`
#' @return ggplot object showing the r-hats of the model to be compared with the
#' convergence criteria (horizontal dashed line)
#' @examples
#' data(chagas2012)
#' seromodel <- fit_seromodel(
#'   serosurvey = chagas2012,
#'   model_type = "time",
#'   foi_index = data.frame(
#'     year = 1935:2011,
#'     foi_index = c(rep(1, 46), rep(2, 31))
#'   ),
#'   iter = 100,
#'   chains = 2
#' )
#' plot_rhats(seromodel, chagas2012)
#' @export
plot_rhats <- function(
  seromodel,
  serosurvey,
  size_text = 11,
  plot_constant = FALSE,
  x_axis = NA
) {
  checkmate::assert_class(seromodel, "stanfit", null.ok = TRUE)

  model_name <- seromodel@model_name
  stopifnot(
    "seromodel@name should start with either 'age', 'time' or 'constant' " =
    startsWith(model_name, "age") | startsWith(model_name, "time") |
    startsWith(model_name, "constant")
  )

  error_msg_x_axis <- paste0(
    "x_axis specification as either 'age' or 'time' when plotting rhat ",
    "estimates of constant models is required to avoid ambiguity"
  )
  validate_plot_constant(
    plot_constant = plot_constant,
    x_axis = x_axis,
    model_name = model_name,
    error_msg_x_axis = error_msg_x_axis
  )

  rhats <- bayesplot::rhat(seromodel, "foi_expanded")

  if (startsWith(model_name, "age") || (plot_constant && (x_axis == "age"))) {
    xlab <- "Age"
    ages <- 1:max(serosurvey$age_max)
    rhats_df <- data.frame(
      age = ages,
      rhat = rhats
    )

    rhats_plot <- ggplot2::ggplot(
      data = rhats_df, ggplot2::aes(x = .data$age)
    )
  } else if (
    startsWith(model_name, "time") ||
    (plot_constant && x_axis == "time")
    ) {
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
      ),
      default = TRUE
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
#' @param plot_constant boolean specifying whether to plot single FOI estimate
#' and its corresponding rhat value instead of showing this information in the
#' summary. Only relevant when `seromodel@model_name == "constant"`)
#' @return ggplot object with a summary of the specified model
#' @examples
#' data(veev2012)
#' seromodel <- fit_seromodel(veev2012, iter = 100)
#' plot_summary(seromodel, veev2012)
#' @export
plot_summary <- function(
  seromodel,
  serosurvey,
  loo_estimate_digits = 1,
  central_estimate_digits = 2,
  rhat_digits = 2,
  size_text = 11,
  plot_constant = FALSE
) {
  checkmate::assert_class(seromodel, "stanfit", null.ok = TRUE)

  summary_list <- summarise_seromodel(
    seromodel = seromodel,
    serosurvey = serosurvey,
    loo_estimate_digits = loo_estimate_digits,
    central_estimate_digits = central_estimate_digits,
    rhat_digits = rhat_digits
  )

  if (plot_constant) {
    if (startsWith(seromodel@model_name, "constant")) {
      drop <- c("foi", "foi_rhat")
      summary_list <- summary_list[!(names(summary_list) %in% drop)]
    } else {
      error_msg <- paste0(
        "plot_constant is only relevant when ",
        "`seromodel@model_name == 'constant'`"
      )
      stop(error_msg)
    }
  }

  summary_df <- data.frame(
    row = rev(seq_len(NCOL(t(summary_list)))),
    text = paste0(colnames(t(summary_list)), ": ", summary_list)
  )

  summary_plot <- ggplot2::ggplot(
    summary_df,
    ggplot2::aes(x = 1, y = .data$row)) +
    ggplot2::scale_y_continuous(
      limits = c(0, nrow(summary_df) + 1),
      breaks = NULL
    ) +
    ggplot2::theme_void() +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$text),
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
#' @param seroreversion_digits Number of seroreversion rate digits
#' @return seromodel summary plot
#' @examples
#' data(veev2012)
#' seromodel <- fit_seromodel(veev2012, iter = 100)
#' plot_seromodel(seromodel, veev2012)
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
  size_text = 11,
  plot_constant = FALSE,
  x_axis = NA
) {
  checkmate::assert_class(seromodel, "stanfit", null.ok = TRUE)

  model_name <- seromodel@model_name
  error_msg_x_axis <- paste0(
    "x_axis specification as either 'age' or 'time' when plotting ",
    "FOI and rhat estimates of constant models is required to avoid ambiguity"
  )
  validate_plot_constant(
    plot_constant = plot_constant,
    x_axis = x_axis,
    model_name = model_name,
    error_msg_x_axis = error_msg_x_axis
  )

  summary_plot <- plot_summary(
    seromodel = seromodel,
    serosurvey = serosurvey,
    loo_estimate_digits = loo_estimate_digits,
    central_estimate_digits = central_estimate_digits,
    rhat_digits = rhat_digits,
    size_text = size_text,
    plot_constant = plot_constant
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

  # This condition (!p | q) is equivalent to !(p & !q)
  # This is to preserve the default behavior for single FOI estimation for
  # constant models
  if (!startsWith(model_name, "constant") || plot_constant) {
    foi_plot <- plot_foi_estimates(
      seromodel,
      serosurvey,
      alpha = alpha,
      foi_df = foi_df,
      foi_max = foi_max,
      size_text,
      plot_constant = plot_constant,
      x_axis = x_axis
    )

    rhats_plot <- plot_rhats(
      seromodel = seromodel,
      serosurvey = serosurvey,
      size_text = size_text,
      plot_constant = plot_constant,
      x_axis = x_axis
    )

    plot_list <- c(
      plot_list,
      list(foi_plot, rhats_plot)
    )
  }

  seromodel_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = 1)
  return(seromodel_plot)
}
