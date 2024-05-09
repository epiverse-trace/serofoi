#' Generate seropositivity plot from a raw serological
#' survey dataset
#'
#' @inheritParams prepare_serodata
#' @inheritParams get_prev_expanded
#' @param size_text Text size use in the theme of the graph returned by the
#'   function.
#' @return A ggplot object containing the seropositivity-vs-age graph of the raw
#'   data of a given seroprevalence survey with its corresponding binomial
#'   confidence interval.
#' @examples
#' data(chagas2012)
#' serodata <- prepare_serodata(chagas2012)
#' plot_seroprev(serodata, size_text = 15)
#' @export
plot_seroprev <- function(serodata,
                          size_text = 6,
                          bin_data = TRUE,
                          bin_step = 5) {
  serodata <- validate_prepared_serodata(serodata = serodata)
  if (bin_data) {
    if (any(serodata$age_max - serodata$age_min > 2)) {
      warn_msg <- paste0(
        "Make sure `serodata` is already grouped by age. ",
        "Skipping binning in seroprevalence plotting."
      )
      warning(warn_msg)
      bin_data <- FALSE
      } else {
        checkmate::assert_int(bin_step, lower = 2)
      }
    }

  if (bin_data) {
    serodata <- prepare_bin_data(
      serodata,
      bin_step = bin_step
      )
    }

  min_prev <- min(serodata$prev_obs_lower)
  max_prev <- max(serodata$prev_obs_upper)

  seroprev_plot <- ggplot2::ggplot(
    data = serodata,
    ggplot2::aes(x = .data$age_mean_f)
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = .data$prev_obs_lower,
        ymax = .data$prev_obs_upper
        ),
        width = 0.1
      ) +
    ggplot2::geom_point(
      ggplot2::aes(
        y = .data$prev_obs,
        size = .data$total
        ),
      fill = "#7a0177", colour = "black", shape = 21
  ) +
    ggplot2::coord_cartesian(
      xlim = c(min(serodata$age_min), max(serodata$age_max)),
      ylim = c(min_prev, max_prev)
      ) +
    ggplot2::theme_bw(size_text)  +
    ggplot2::theme(legend.position = "none") +
    ggplot2::ylab("seropositivity") +
    ggplot2::xlab("age")

  return(seroprev_plot)
}

#' Generate seropositivity plot corresponding to the specified
#' fitted serological model
#'
#' This function generates a seropositivity plot of the specified serological
#' model object. This includes the original data grouped by age as well as the
#' obtained fitting from the model implementation. Age is located on the x axis
#' and seropositivity on the y axis with its corresponding confidence interval.
#' @inheritParams get_foi_central_estimates
#' @inheritParams run_seromodel
#' @inheritParams get_prev_expanded
#' @param size_text Text size of the graph returned by the function.
#' @return A ggplot object containing the seropositivity-vs-age graph including
#'   the data, the fitted model and their corresponding confidence intervals.
#' @examples
#' data(chagas2012)
#' serodata <- prepare_serodata(chagas2012)
#' seromodel_object <- fit_seromodel(
#'   serodata = serodata,
#'   foi_model = "constant",
#'   iter = 1000
#' )
#' plot_seroprev_fitted(seromodel_object,
#'   serodata = serodata,
#'   size_text = 15
#' )
#' @export
plot_seroprev_fitted <- function(seromodel_object,
                                 serodata,
                                 size_text = 6,
                                 bin_data = TRUE,
                                 bin_step = 5,
                                 alpha = 0.05
                                 ) {
  checkmate::assert_class(seromodel_object, "stanfit", null.ok = TRUE)
  validate_prepared_serodata(serodata)

  foi <- rstan::extract(seromodel_object, "foi", inc_warmup = FALSE)[[1]]

  prev_expanded <- get_prev_expanded(
    foi,
    serodata = serodata,
    alpha = alpha,
    bin_data = bin_data,
    bin_step = bin_step
  )
  prev_plot <-
    ggplot2::ggplot(prev_expanded, ggplot2::aes(x = .data$age)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = .data$predicted_prev_lower,
        ymax = .data$predicted_prev_upper
      ),
      fill = "#c994c7"
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$predicted_prev),
      colour = "#7a0177"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$prev_obs_lower, ymax = .data$prev_obs_upper),
      width = 0.1
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = .data$prev_obs, size = .data$total),
      fill = "#7a0177",
      colour = "black",
      shape = 21
    ) +
    ggplot2::theme_bw(size_text) +
    ggplot2::coord_cartesian(
      xlim = c(min(serodata$age_min), max(serodata$age_max)),
      ylim = c(
        min(prev_expanded$prev_obs_lower, prev_expanded$predicted_prev_lower),
        max(prev_expanded$prev_obs_upper, prev_expanded$predicted_prev_upper)
        )
    ) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::ylab("seropositivity") +
    ggplot2::xlab("age")

  return(prev_plot)
}

# TODO Complete @param documentation

#' Generate force-of-infection plot corresponding to the
#' specified fitted serological model
#'
#' This function generates a force-of-infection plot from the results obtained
#' by fitting a serological model. This includes the corresponding binomial
#' confidence interval. The x axis corresponds to the decades covered by the
#' survey the y axis to the force-of-infection.
#' @inheritParams get_foi_central_estimates
#' @param size_text Text size use in the theme of the graph returned by the
#'   function.
#' @param max_lambda TBD
#' @param foi_sim TBD
#' @return A ggplot2 object containing the force-of-infection vs time including
#'   the corresponding confidence interval.
#' @examples
#' data(chagas2012)
#' serodata <- prepare_serodata(chagas2012)
#' seromodel_object <- fit_seromodel(
#'   serodata = serodata,
#'   foi_model = "constant",
#'   iter = 1000
#' )
#' cohort_ages <- get_cohort_ages(serodata)
#' plot_foi(
#'   seromodel_object = seromodel_object,
#'   cohort_ages = cohort_ages,
#'   size_text = 15
#' )
#' @export
plot_foi <- function(seromodel_object,
                     cohort_ages,
                     max_lambda = NA,
                     size_text = 25,
                     foi_sim = NULL) {
  checkmate::assert_class(seromodel_object, "stanfit", null.ok = TRUE)
  #-------- This bit is to get the actual length of the foi data
  foi_data <- get_foi_central_estimates(
    seromodel_object = seromodel_object,
    cohort_ages = cohort_ages
  )

  foi_plot <-
    ggplot2::ggplot(foi_data, ggplot2::aes(x = .data$year)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = .data$lower,
        ymax = .data$upper
      ),
      fill = "#41b6c4",
      alpha = 0.5
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$medianv),
      colour = "#253494",
      size = size_text / 8
    ) +
    ggplot2::theme_bw(size_text) +
    ggplot2::coord_cartesian(ylim = c(0, max_lambda)) +
    ggplot2::ylab("Force-of-Infection") +
    ggplot2::xlab("year")

  if (!is.null(foi_sim)) {
    if (nrow(foi_data) != length(foi_sim)) {
      warn_msg <- paste0(
        "`foi_sim` has different length than `exposure_years`. ",
        "Dropping last elements of `foi_sim`"
      )
      warning(warn_msg)
      remove_x_values <- length(foi_sim) - nrow(foi_data)
      foi_sim_data <- data.frame(
        year = foi_data$year,
        foi_sim = foi_sim[-(1:remove_x_values)]
      )
      foi_plot <- foi_plot +
        ggplot2::geom_line(
          data = foi_sim_data, ggplot2::aes(y = foi_sim),
          colour = "#b30909",
          size = size_text / 8
        )
    } else {
      foi_sim_data <- data.frame(
        year = foi_data$year,
        foi_sim = foi_sim
      )
      foi_plot <- foi_plot +
        ggplot2::geom_line(
          data = foi_sim_data, ggplot2::aes(y = foi_sim),
          colour = "#b30909",
          size = size_text / 8
        )
    }
  }

  return(foi_plot)
}

#' Generate plot of the R-hat estimates for the specified fitted
#' serological model
#'
#' This function generates a plot of the R-hat estimates obtained for a
#' specified fitted serological model `seromodel_object`. The x axis corresponds
#' to the decades covered by the survey and the y axis to the value of the
#' rhats. All rhats must be smaller than 1 to ensure convergence (for further
#' details check [rhat][bayesplot::rhat]).
#' @inheritParams get_foi_central_estimates
#' @param size_text Text size use in the theme of the graph returned by the
#'   function.
#' @return The rhats-convergence plot of the selected model.
#' @examples
#' data(chagas2012)
#' serodata <- prepare_serodata(chagas2012)
#' seromodel_object <- fit_seromodel(
#'   serodata = serodata,
#'   foi_model = "constant",
#'   iter = 1000
#' )
#' cohort_ages <- get_cohort_ages(serodata = serodata)
#' plot_rhats(seromodel_object,
#'   cohort_ages = cohort_ages,
#'   size_text = 15
#' )
#' @export
plot_rhats <- function(seromodel_object,
                       cohort_ages,
                       size_text = 25) {
  checkmate::assert_class(seromodel_object, "stanfit", null.ok = TRUE)

  rhats <- get_table_rhats(
    seromodel_object = seromodel_object,
    cohort_ages = cohort_ages
  )

  rhats_plot <-
    ggplot2::ggplot(rhats, ggplot2::aes(.data$year, .data$rhat)) +
    ggplot2::geom_line(colour = "purple") +
    ggplot2::geom_point() +
    ggplot2::coord_cartesian(
      ylim = c(
        min(1.0, min(rhats$rhat)),
        max(1.02, max(rhats$rhat))
        )
      ) +
    ggplot2::geom_hline(
      yintercept = 1.01,
      colour = "blue",
      size = size_text / 12
    ) +
    ggplot2::theme_bw(size_text) +
    ggplot2::ylab("Convergence (R^)")

  return(rhats_plot)
}

#' Generate vertical arrangement of plots showing a summary of a
#' model, the estimated seroprevalence, the force-of-infection fit and the R-hat
#' estimates plots.
#'
#' @inheritParams get_foi_central_estimates
#' @inheritParams run_seromodel
#' @inheritParams get_prev_expanded
#' @param size_text Text size use in the theme of the graph returned by the
#'   function.
#' @param max_lambda TBD
#' @param foi_sim TBD
#' @return A ggplot object with a vertical arrange containing the
#'   seropositivity, force of infection, and convergence plots.
#' @examples
#' data(chagas2012)
#' serodata <- prepare_serodata(chagas2012)
#' seromodel_object <- fit_seromodel(
#'   serodata = serodata,
#'   foi_model = "constant",
#'   iter = 1000
#' )
#' plot_seromodel(seromodel_object,
#'   serodata = serodata,
#'   size_text = 15
#' )
#' @export
plot_seromodel <- function(seromodel_object,
                           serodata,
                           alpha = 0.05,
                           max_lambda = NA,
                           size_text = 25,
                           bin_data = TRUE,
                           bin_step = 5,
                           foi_sim = NULL) {
  checkmate::assert_class(seromodel_object, "stanfit", null.ok = TRUE)
  serodata <- validate_serodata(serodata)

  cohort_ages <- get_cohort_ages(serodata = serodata)

  prev_plot <- plot_seroprev_fitted(
    seromodel_object = seromodel_object,
    serodata = serodata,
    alpha = alpha,
    size_text = size_text,
    bin_data = bin_data,
    bin_step = bin_step
  )

  foi_plot <- plot_foi(
    seromodel_object = seromodel_object,
    cohort_ages = cohort_ages,
    max_lambda = max_lambda,
    size_text = size_text,
    foi_sim = foi_sim
  )

  rhats_plot <- plot_rhats(
    seromodel_object = seromodel_object,
    cohort_ages = cohort_ages,
    size_text = size_text
  )

  model_summary <- extract_seromodel_summary(
    seromodel_object = seromodel_object,
    serodata = serodata
  )
  summary_table <- t(
    dplyr::select(
      model_summary,
      c("foi_model", "dataset", "elpd", "se", "converged")
    )
  )
  summary_plot <-
    plot_info_table(summary_table, size_text = size_text)

  plot_arrange <- cowplot::plot_grid(
    summary_plot,
    prev_plot,
    foi_plot,
    rhats_plot,
    ncol = 1,
    nrow = 4,
    rel_heights = c(0.5, 1, 1, 1)
  )
  return(plot_arrange)
}

#' Generate plot summarizing a given table
#'
#' @param info_table Table with the information to be summarised
#' @param size_text Text size of the graph returned by the function
#' @return ggplot object summarizing the information in `info_table`
#' @examples
#' serodata <- prepare_serodata(chagas2012)
#' seromodel_object <- fit_seromodel(
#'   serodata = serodata,
#'   foi_model = "constant",
#'   iter = 1000
#' )
#' seromodel_summary <- extract_seromodel_summary(
#'   seromodel_object = seromodel_object,
#'   serodata = serodata
#' )
#' info_table <- t(seromodel_summary)
#' plot_info_table(info_table, size_text = 15)
#' @export
plot_info_table <- function(info_table, size_text) {
  dato <- data.frame(
    y = NROW(info_table):seq_len(1),
    text = paste0(rownames(info_table), ": ", info_table[, 1])
  )
  p <- ggplot2::ggplot(dato, ggplot2::aes(x = 1, y = .data$y)) +
    ggplot2::scale_y_continuous(
      limits = c(0, NROW(info_table) + 1),
      breaks = NULL
      ) +
    ggplot2::theme_void() +
    ggplot2::geom_text(ggplot2::aes(label = text),
      size = size_text / 2.5,
      fontface = "bold"
    )

  return(p)
}
