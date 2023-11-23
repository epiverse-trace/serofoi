#' Function that generates the sero-positivity plot from a raw serological
#' survey dataset
#'
#' @inheritParams prepare_serodata
#' @param size_text Text size use in the theme of the graph returned by the
#'   function.
#' @return A ggplot object containing the seropositivity-vs-age graph of the raw
#'   data of a given seroprevalence survey with its corresponging binomial
#'   confidence interval.
#' @examples
#' data(chagas2012)
#' serodata <- prepare_serodata(chagas2012)
#' plot_seroprev(serodata, size_text = 15)
#' @export
plot_seroprev <- function(serodata,
                          size_text = 6) {
  xx <- prepare_bin_data(serodata)
  seroprev_plot <-
    ggplot2::ggplot(data = xx) +
    ggplot2::geom_errorbar(
      ggplot2::aes(age, ymin = p_obs_bin_l, ymax = p_obs_bin_u),
      width = 0.1
    ) +
    ggplot2::geom_point(
      ggplot2::aes(age, p_obs_bin, size = xx$bin_size),
      fill = "#7a0177", colour = "black", shape = 21
    ) +
    ggplot2::theme_bw(size_text) +
    ggplot2::coord_cartesian(xlim = c(0, 60), ylim = c(0, 1)) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::ylab("Sero-positivity") +
    ggplot2::xlab("Age")

  return(seroprev_plot)
}

#' Function that generates a seropositivity plot corresponding to the specified
#' fitted serological model
#'
#' This function generates a seropositivity plot of the specified serological
#' model object. This includes the original data grouped by age as well as the
#' obtained fitting from the model implementation. Age is located on the x axis
#' and seropositivity on the y axis with its corresponding confidence interval.
#' @inheritParams get_foi_central_estimates
#' @inheritParams run_seromodel
#' @param size_text Text size of the graph returned by the function.
#' @return A ggplot object containing the seropositivity-vs-age graph including
#'   the data, the fitted model and their corresponding confindence intervals.
#' @examples
#' data(chagas2012)
#' serodata <- prepare_serodata(chagas2012)
#' seromodel_object <- run_seromodel(
#'   serodata = serodata,
#'   foi_model = "constant",
#'   n_iters = 1000
#' )
#' plot_seroprev_fitted(seromodel_object,
#'   serodata = serodata,
#'   size_text = 15
#' )
#' @export
plot_seroprev_fitted <- function(seromodel_object,
                                 serodata,
                                 size_text = 6) {
  if (!is.character(seromodel_object)) {
    if (class(seromodel_object@sim$samples) != "NULL") {
      foi <- rstan::extract(seromodel_object, "foi", inc_warmup = FALSE)[[1]]
      prev_expanded <- get_prev_expanded(
        foi,
        serodata = serodata,
        bin_data = TRUE
      )
      prev_plot <-
        ggplot2::ggplot(prev_expanded, ggplot2::aes(x = age)) +
        ggplot2::geom_ribbon(
          ggplot2::aes(
            ymin = predicted_prev_lower,
            ymax = predicted_prev_upper
          ),
          fill = "#c994c7"
        ) +
        ggplot2::geom_line(
          ggplot2::aes(y = predicted_prev),
          colour = "#7a0177"
        ) +
        ggplot2::geom_errorbar(
          ggplot2::aes(ymin = p_obs_bin_l, ymax = p_obs_bin_u),
          width = 0.1
        ) +
        ggplot2::geom_point(
          ggplot2::aes(y = p_obs_bin, size = bin_size),
          fill = "#7a0177",
          colour = "black",
          shape = 21
        ) +
        ggplot2::theme_bw(size_text) +
        ggplot2::coord_cartesian(xlim = c(0, 60), ylim = c(0, 1)) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::ylab("Sero-positivity") +
        ggplot2::xlab("Age")
    }
  } else {
    message("model did not run")
    print_warning <- "errors"
    df <- data.frame()

    prev_plot <- ggplot2::ggplot(df) +
      ggplot2::geom_point() +
      ggplot2::xlim(0, 10) +
      ggplot2::ylim(0, 10) +
      ggplot2::annotate("text",
        x = 4,
        y = 5,
        label = print_warning
      ) +
      ggplot2::theme_bw(25) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank()
      ) +
      ggplot2::ylab(" ") +
      ggplot2::xlab(" ")
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10))
  }

  return(prev_plot)
}

# TODO Complete @param documentation

#' Function that generates a Force-of-Infection plot corresponding to the
#' specified fitted serological model
#'
#' This function generates a Force-of-Infection plot from the results obtained
#' by fitting a serological model. This includes the corresponding binomial
#' confidence interval. The x axis corresponds to the decades covered by the
#' survey the y axis to the Force-of-Infection.
#' @inheritParams get_foi_central_estimates
#' @param size_text Text size use in the theme of the graph returned by the
#'   function.
#' @param max_lambda TBD
#' @param foi_sim TBD
#' @return A ggplot2 object containing the Force-of-infection-vs-time including
#'   the corresponding confidence interval.
#' @examples
#' data(chagas2012)
#' serodata <- prepare_serodata(chagas2012)
#' seromodel_object <- run_seromodel(
#'   serodata = serodata,
#'   foi_model = "constant",
#'   n_iters = 1000
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
  if (!is.character(seromodel_object)) {
    if (class(seromodel_object@sim$samples) != "NULL") {
      foi <- rstan::extract(seromodel_object,
        "foi",
        inc_warmup = FALSE
      )[[1]]
      #-------- This bit is to get the actual length of the foi data
      foi_data <- get_foi_central_estimates(
        seromodel_object = seromodel_object,
        cohort_ages = cohort_ages
      )

      #--------
      foi_data$medianv[1] <- NA
      foi_data$lower[1] <- NA
      foi_data$upper[1] <- NA

      foi_plot <-
        ggplot2::ggplot(foi_data) +
        ggplot2::geom_ribbon(
          ggplot2::aes(
            x = year,
            ymin = lower,
            ymax = upper
          ),
          fill = "#41b6c4",
          alpha = 0.5
        ) +
        ggplot2::geom_line(ggplot2::aes(x = year, y = medianv),
          colour = "#253494",
          size = size_text / 8
        ) +
        ggplot2::theme_bw(size_text) +
        ggplot2::coord_cartesian(ylim = c(0, max_lambda)) +
        ggplot2::ylab("Force-of-Infection") +
        ggplot2::xlab("Year")
      # TODO Add warning for foi_sim of different length than exposure years
      if (!is.null(foi_sim)) {
        if (nrow(foi_data) != length(foi_sim)) {
          remove_x_values <- length(foi_sim) - nrow(foi_data)
          foi_sim_data <- data.frame(
            year = foi_data$year,
            foi_sim = foi_sim[-(1:remove_x_values)]
          )
          foi_plot <- foi_plot +
            ggplot2::geom_line(
              data = foi_sim_data, ggplot2::aes(x = year, y = foi_sim),
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
              data = foi_sim_data, ggplot2::aes(x = year, y = foi_sim),
              colour = "#b30909",
              size = size_text / 8
            )
        }
      }
    }
  } else {
    message("model did not run")
    print_warning <- "errors"
    df <- data.frame()

    foi_plot <- ggplot2::ggplot(df) +
      ggplot2::geom_point() +
      ggplot2::xlim(0, 10) +
      ggplot2::ylim(0, 10) +
      ggplot2::annotate("text",
        x = 4,
        y = 5,
        label = print_warning
      ) +
      ggplot2::theme_bw(25) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank()
      ) +
      ggplot2::ylab(" ") +
      ggplot2::xlab(" ")
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10))
  }

  return(foi_plot)
}

#' Function that generates a plot of the R-hat estimates of the specified fitted
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
#' seromodel_object <- run_seromodel(
#'   serodata = serodata,
#'   foi_model = "constant",
#'   n_iters = 1000
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
  if (!is.character(seromodel_object)) {
    if (class(seromodel_object@sim$samples) != "NULL") {
      rhats <- get_table_rhats(
        seromodel_object = seromodel_object,
        cohort_ages = cohort_ages
      )

      rhats_plot <-
        ggplot2::ggplot(rhats, ggplot2::aes(year, rhat)) +
        ggplot2::geom_line(colour = "purple") +
        ggplot2::geom_point() +
        ggplot2::coord_cartesian(ylim = c(0.7, 2)) +
        ggplot2::geom_hline(
          yintercept = 1.1,
          colour = "blue",
          size = size_text / 12
        ) +
        ggplot2::theme_bw(size_text) +
        ggplot2::ylab("Convergence (R^)")
    }
  } else {
    message("model did not run")
    print_warning <- "errors"
    df <- data.frame()

    rhats_plot <- ggplot2::ggplot(df) +
      ggplot2::geom_point() +
      ggplot2::xlim(0, 10) +
      ggplot2::ylim(0, 10) +
      ggplot2::annotate("text",
        x = 4,
        y = 5,
        label = print_warning
      ) +
      ggplot2::theme_bw(25) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank()
      ) +
      ggplot2::ylab(" ") +
      ggplot2::xlab(" ")
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10))
  }

  return(rhats_plot)
}

# TODO Complete @param documentation

#' Function that generates a vertical arrange of plots showing a summary of a
#' model, the estimated seroprevalence, the Force-of-Infection fit and the R-hat
#' estimates plots.
#'
#' @inheritParams get_foi_central_estimates
#' @inheritParams run_seromodel
#' @param size_text Text size use in the theme of the graph returned by the
#'   function.
#' @param max_lambda TBD
#' @param foi_sim TBD
#' @return A ggplot object with a vertical arrange containing the
#'   seropositivity, force of infection, and convergence plots.
#' @examples
#' data(chagas2012)
#' serodata <- prepare_serodata(chagas2012)
#' seromodel_object <- run_seromodel(
#'   serodata = serodata,
#'   foi_model = "constant",
#'   n_iters = 1000
#' )
#' plot_seromodel(seromodel_object,
#'   serodata = serodata,
#'   size_text = 15
#' )
#' @export
plot_seromodel <- function(seromodel_object,
                           serodata,
                           max_lambda = NA,
                           size_text = 25,
                           foi_sim = NULL) {
  if (!is.character(seromodel_object)) {
    if (class(seromodel_object@sim$samples) != "NULL") {
      cohort_ages <- get_cohort_ages(serodata = serodata)

      prev_plot <- plot_seroprev_fitted(
        seromodel_object = seromodel_object,
        serodata = serodata,
        size_text = size_text
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
    }
  } else {
    message("model did not run")
    print_warning <- "errors"
    df <- data.frame()

    g0 <- ggplot2::ggplot(df) +
      ggplot2::geom_point() +
      ggplot2::xlim(0, 10) +
      ggplot2::ylim(0, 10) +
      ggplot2::annotate("text",
        x = 4,
        y = 5,
        label = print_warning
      ) +
      ggplot2::theme_bw(25) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank()
      ) +
      ggplot2::ylab(" ") +
      ggplot2::xlab(" ")
    g1 <- g0
    # TODO: This
    g0 <- g0 + ggplot2::labs(subtitle = seromodel_object$model_name) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 10))

    plot_arrange <-
      cowplot::plot_grid(g0, g1, g1, g1, g1, ncol = 1, nrow = 5)
  }

  return(plot_arrange)
}

# TODO Improve documentation of @return.
# TODO Give more details about the generated plot
#' Function that generates a plot for a given table
#'
#' @param info the information that will be contained in the table
#' @param size_text Text size of the graph returned by the function
#' @return p the plot for the given table
#' @examples
#' serodata <- prepare_serodata(chagas2012)
#' seromodel_object <- run_seromodel(
#'   serodata = serodata,
#'   foi_model = "constant",
#'   n_iters = 1000
#' )
#' seromodel_summary <- extract_seromodel_summary(
#'   seromodel_object = seromodel_object,
#'   serodata = serodata
#' )
#' info <- t(seromodel_summary)
#' plot_info_table(info, size_text = 15)
#' @export
plot_info_table <- function(info, size_text) {
  dato <- data.frame(
    y = NROW(info):seq_len(1),
    text = paste0(rownames(info), ": ", info[, 1])
  )
  p <- ggplot2::ggplot(dato, ggplot2::aes(x = 1, y = y)) +
    ggplot2::scale_y_continuous(limits = c(0, NROW(info) + 1), breaks = NULL) +
    ggplot2::theme_void() +
    ggplot2::geom_text(ggplot2::aes(label = text),
      size = size_text / 2.5,
      fontface = "bold"
    )

  return(p)
}
