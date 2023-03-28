#' Function that generates the sero-positivity plot from a raw serological survey dataset.
#'
#' @param serodata A data frame containing the data from a seroprevalence survey.
#' This data frame must contain the following columns:
#' \tabular{ll}{
#' \code{survey} \tab survey Label of the current survey \cr \tab \cr
#' \code{total} \tab Number of samples for each age group\cr \tab \cr
#' \code{counts} \tab Number of positive samples for each age group\cr \tab \cr
#' \code{age_min} \tab age_min \cr \tab \cr
#' \code{age_max} \tab age_max \cr \tab \cr
#' \code{year_init} \tab year_init \cr \tab \cr
#' \code{year_end} \tab year_end \cr \tab \cr
#' \code{tsur} \tab Year in which the survey took place \cr \tab \cr
#' \code{country} \tab The country where the survey took place \cr \tab \cr
#' \code{test} \tab The type of test taken \cr \tab \cr
#' \code{antibody} \tab antibody \cr \tab \cr
#' }
#' @param size_text Text size use in the theme of the graph returned by the function.
#' @return A ggplot object containing the seropositivity-vs-age graph of the raw data of a given seroprevalence survey with its corresponging binomial confidence interval.
#' @examples
#' \dontrun{
#'  data_test <- prepare_serodata(serodata)
#'  seromodel_object <- run_seromodel(
#'  serodata = data_test,
#'  seromodel_name = "constant_foi_bi",
#'  n_iters = 1000
#')
#' plot_seroprev(seromodel_object, size_text = 15)
#' }
#' @export
plot_seroprev <- function(serodata,
                          size_text = 6) {
  xx <- prepare_bin_data(serodata)
  seroprev_plot <-
    ggplot2::ggplot(data = xx) +
    ggplot2::geom_errorbar(ggplot2::aes(age, ymin = p_obs_bin_l, ymax = p_obs_bin_u), width = 0.1) +
    ggplot2::geom_point(ggplot2::aes(age, p_obs_bin, size = xx$bin_size), fill = "#7a0177", colour = "black", shape = 21) +
    ggplot2::theme_bw(size_text) +
    ggplot2::coord_cartesian(xlim = c(0, 60), ylim = c(0,1)) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::ylab("Sero-positivity") + ggplot2::xlab("Age")

  return(seroprev_plot)

}

#' Function that generates a seropositivity plot corresponding to the specified fitted serological model.
#'
#' This function generates a seropositivity plot of the specified serological model object. This includes the original data grouped by age
#' as well as the obtained fitting from the model implementation. Age is located on the x axis and seropositivity on the y axis with its 
#' corresponding confidence interval.
#' @param seromodel_object Object containing the results of fitting a model by means of \link{run_seromodel}.
#' @param size_text Text size of the graph returned by the function.
#' @return A ggplot object containing the seropositivity-vs-age graph including the data, the fitted model and their corresponding confindence intervals.
#' @examples
#' \dontrun{
#' data("serodata")
#' data_test <- prepare_serodata(serodata)
#' seromodel_object <- run_seromodel(serodata = data_test,
#'                                   seromodel_name = "constant_foi_bi",
#'                                   n_iters = 1000)
#' plot_seroprev_fitted(seromodel_object, size_text = 15)
#' }
#' @export
plot_seroprev_fitted <- function(seromodel_object,
                                 size_text = 6) {

  if (is.character(seromodel_object$fit) == FALSE)  {
    if  (class(seromodel_object$fit@sim$samples)  != "NULL" ) {

      foi <- rstan::extract(seromodel_object$fit, "foi", inc_warmup = FALSE)[[1]]
      prev_expanded <- get_prev_expanded(foi, serodata = seromodel_object$serodata)
      prev_plot <-
        ggplot2::ggplot(prev_expanded) +
        ggplot2::geom_ribbon(
          ggplot2::aes(
            x = age,
            ymin = predicted_prev_lower,
            ymax = predicted_prev_upper
          ),
          fill = "#c994c7"
        ) +
        ggplot2::geom_line(ggplot2::aes(x = age, y = predicted_prev), colour = "#7a0177") +
        ggplot2::geom_errorbar(ggplot2::aes(age, ymin = p_obs_bin_l, ymax = p_obs_bin_u),
                               width = 0.1) +
        ggplot2::geom_point(
          ggplot2::aes(age, p_obs_bin, size = bin_size),
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
    print("model did not run")
    print_warning <- "errors"
    df <- data.frame()

    prev_plot <- ggplot2::ggplot(df) +
      ggplot2::geom_point() +
      ggplot2::xlim(0, 10) +
      ggplot2::ylim(0, 10) +
      ggplot2::annotate("text",
                        x = 4,
                        y = 5,
                        label = print_warning) +
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

#' Function that generates a Force-of-Infection plot corresponding to the specified fitted serological model.
#'
#' This function generates a Force-of-Infection plot from the results obtained by fitting a serological model.
#' This includes the corresponding binomial confidence interval. 
#' The x axis corresponds to the decades covered by the survey the y axis to the Force-of-Infection.
#' @param seromodel_object Object containing the results of fitting a model by means of \link{run_seromodel}.
#' @param size_text Text size use in the theme of the graph returned by the function.
#' @return A ggplot2 object containing the Force-of-infection-vs-time including the corresponding confidence interval.
#' @examples
#' \dontrun{
#'    data_test <- prepare_serodata(serodata)
#'    seromodel_object <- run_seromodel(
#'    serodata = data_test,
#'    seromodel_name = "constant_foi_bi",
#'    n_iters = 1000
#' )
#' plot_foi(seromodel_object, size_text = 15)
#' }
#' @export
plot_foi <- function(seromodel_object,
                     lambda_sim = NA,
                     max_lambda = NA,
                     size_text = 25) {
  if (is.character(seromodel_object$fit) == FALSE) {
    if (class(seromodel_object$fit@sim$samples) != "NULL") {
      foi <- rstan::extract(seromodel_object$fit,
                            "foi",
                            inc_warmup = FALSE)[[1]]

      #-------- This bit is to get the actual length of the foi data
      foi_data <- seromodel_object$foi_cent_est

      if (!is.na(lambda_sim)) {
        lambda_mod_length <- NROW(foi_data)
        lambda_sim_length <- length(lambda_sim)

        if (lambda_mod_length < lambda_sim_length) {
          remove_x_values <- lambda_sim_length - lambda_mod_length
          lambda_sim <- lambda_sim[-c(1:remove_x_values)]
        }

        foi_data$simulated <- lambda_sim
      }

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
                           size = size_text / 8) +
        ggplot2::theme_bw(size_text) +
        ggplot2::coord_cartesian(ylim = c(0, max_lambda)) +
        ggplot2::ylab("Force-of-Infection") +
        ggplot2::xlab("Year")
    }
  } else {
    print("model did not run")
    print_warning <- "errors"
    df <- data.frame()

    foi_plot <- ggplot2::ggplot(df) +
      ggplot2::geom_point() +
      ggplot2::xlim(0, 10) +
      ggplot2::ylim(0, 10) +
      ggplot2::annotate("text",
                        x = 4,
                        y = 5,
                        label = print_warning) +
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

#' Function that generates a plot of the R-hat estimates of the specified fitted serological model.
#'
#' This function generates a plot of the R-hat estimates obtained for a specified fitted serological model \code{seromodel_object}. 
#' The x axis corresponds to the decades covered by the survey and the y axis to the value of the rhats. 
#' All rhats must be smaller than 1 to ensure convergence (for further details check \link[bayesplot]{rhat}).
#' @param seromodel_object Object containing the results of fitting a model by means of \link{run_seromodel}.
#' @param size_text Text size use in the theme of the graph returned by the function.
#' @return The rhats-convergence plot of the selected model.
#' @examples
#' \dontrun{
#' data("serodata")
#' data_test <- prepare_serodata(serodata)
#' seromodel_object <- run_seromodel(
#'  serodata = data_test,
#'  seromodel_name = "constant_foi_bi",
#'  n_iters = 1000
#')
#' plot_rhats(seromodel_object, 
#'            size_text = 15)
#' }
#' @export
plot_rhats <- function(seromodel_object,
                       size_text = 25) {
  if (is.character(seromodel_object$fit) == FALSE) {
    if (class(seromodel_object$fit@sim$samples) != "NULL") {
      rhats <- get_table_rhats(seromodel_object)

      rhats_plot <-
        ggplot2::ggplot(rhats, ggplot2::aes(year, rhat)) +
        ggplot2::geom_line(colour = "purple") +
        ggplot2::geom_point() +
        ggplot2::coord_cartesian(ylim = c(0.7, 2)) +
        ggplot2::geom_hline(yintercept = 1.1,
                            colour = "blue",
                            size = size_text / 12) +
        ggplot2::theme_bw(size_text) +
        ggplot2::ylab("Convergence (R^)")
    }
  } else {
    print("model did not run")
    print_warning <- "errors"
    df <- data.frame()

    rhats_plot <- ggplot2::ggplot(df) +
      ggplot2::geom_point() +
      ggplot2::xlim(0, 10) +
      ggplot2::ylim(0, 10) +
      ggplot2::annotate("text",
                        x = 4,
                        y = 5,
                        label = print_warning) +
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


#' Function that generates a vertical arrange of plots showing a summary of a given implemented model, as well as its corresponding
#' seroprevalence, Force-of-Infection and R-hat estimates plots.
#'
#' @param seromodel_object Object containing the results of fitting a model by means of \link{run_seromodel}.
#' @param size_text Text size use in the theme of the graph returned by the function.
#' @return A ggplot object with a vertical arrange containing the seropositivity, force of infection, and convergence plots.
#' @examples
#' \dontrun{
#' data_test <- prepare_serodata(serodata)
#' seromodel_object <- run_seromodel(
#'  serodata = data_test,
#'  seromodel_name = "constant_foi_bi",
#'  n_iters = 1000
#')
#' plot_seromodel(seromodel_object, size_text = 15)
#' }
#' @export
plot_seromodel <- function(seromodel_object,
                       lambda_sim = NA,
                       max_lambda = NA,
                       size_text = 25) {
  if (is.character(seromodel_object$fit) == FALSE) {
    if (class(seromodel_object$fit@sim$samples) != "NULL") {
      prev_plot <- plot_seroprev_fitted(seromodel_object = seromodel_object,
                                 size_text = size_text)

      foi_plot <- plot_foi(
        seromodel_object = seromodel_object,
        lambda_sim = lambda_sim,
        max_lambda = max_lambda,
        size_text = size_text
      )

      rhats_plot <- plot_rhats(seromodel_object = seromodel_object,
                               size_text = size_text)

      summary_table <- t(
        dplyr::select(seromodel_object$model_summary, 
        c('seromodel_name', 'dataset', 'elpd', 'se', 'converged')))
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
    print("model did not run")
    print_warning <- "errors"
    df <- data.frame()

    g0 <- ggplot2::ggplot(df) +
      ggplot2::geom_point() +
      ggplot2::xlim(0, 10) +
      ggplot2::ylim(0, 10) +
      ggplot2::annotate("text",
                        x = 4,
                        y = 5,
                        label = print_warning) +
      ggplot2::theme_bw(25) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank()
      ) +
      ggplot2::ylab(" ") +
      ggplot2::xlab(" ")
    g1 <- g0
    g0 <- g0 + ggplot2::labs(subtitle = seromodel_object$model) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 10))

    plot_arrange <-
      cowplot::plot_grid(g0, g1, g1, g1, g1, ncol = 1, nrow = 5)
  }

  return(plot_arrange)
}


#' Auxiliary function that generates a plot of a given table
#'
#' @param info the information that will be contained in the table
#' @param size_text Text size of the graph returned by the function
#' @return p, a variable that will be used in the \link{visualisation} module
#' @examples
#' \dontrun{
#'  data_test <- prepare_serodata(serodata)
#'  seromodel_object <- run_seromodel(
#'  serodata = data_test,
#'  seromodel_name = "constant_foi_bi",
#'  n_iters = 1000
#')
#' info = t(seromodel_object$model_summary)
#' plot_info_table (info, size_text = 15)
#' }
#' @export
plot_info_table <- function(info, size_text) {
  dato <- data.frame(y = NROW(info):seq_len(1),
                     text = paste0(rownames(info), ": ", info[, 1]))
  p <- ggplot2::ggplot(dato, ggplot2::aes(x = 1, y = y)) +
    ggplot2::scale_y_continuous(limits = c(0, NROW(info) + 1), breaks = NULL) +
    ggplot2::theme_void() +
    ggplot2::geom_text(ggplot2::aes(label = text),
                       size = size_text / 2.5,
                       fontface = "bold")

  return(p)
}