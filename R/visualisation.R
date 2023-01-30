#' Generate sero-positivity plot from raw data
#'
#' Function that generates the sero positivity plot
#' @param data_test Object that the prepared data_test
#' @param xlabel Label of axis x
#' @param ylabel Label of axis y
#' @return The sero-positivity plot
#' @export
plot_seroprev <- function(data_test, size_text = 6) {

  # OJO!! Aquí está pendiente agregar la función que genera el binned prevalence (revisar plot_seroprev_fitted)

      seroprev_plot <-
        ggplot2::ggplot(data = data_test) +
        ggplot2::geom_errorbar(ggplot2::aes(age_mean_f, ymin = prev_obs_lower, ymax = prev_obs_upper), width = 0.1) +
        ggplot2::geom_point(ggplot2::aes(age_mean_f, prev_obs), fill = "#7a0177", colour = "black", shape = 21) +
        ggplot2::theme_bw(size_text) +
        ggplot2::coord_cartesian(xlim = c(0, 60), ylim = c(0,1)) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::ylab("Sero-positivity") + ggplot2::xlab("Age")



  return(seroprev_plot)

}





#' Generate sero-positivity plot with fitted model
#'
#' Function that generates the sero positivity plot
#' @param model_object Object that the run_model function returns with the results of the fit
#' @param model Refers to the model selected
#' @param xlabel Label of axis x
#' @param ylabel Label of axis y
#' @return The sero-positivity plot
#' @export
plot_seroprev_fitted <- function(model_object,
                          size_text = 6) {

  if (is.character(model_object$fit) == FALSE)  {
    if  (class(model_object$fit@sim$samples)  != "NULL" ) {

      foi <- rstan::extract(model_object$fit, "foi", inc_warmup = FALSE)[[1]]
      prev_expanded <- get_prev_expanded(foi, model_data = model_object$model_data)

      prev_plot <-
        ggplot2::ggplot(prev_expanded) +
        ggplot2::geom_ribbon(ggplot2::aes(x = age, ymin = predicted_prev_lower, ymax = predicted_prev_upper), fill = "#c994c7") +
        ggplot2::geom_line(ggplot2::aes(x = age, y = predicted_prev), colour = "#7a0177") +
        ggplot2::geom_errorbar(ggplot2::aes(age, ymin = p_obs_bin_l, ymax = p_obs_bin_u), width = 0.1) +
        ggplot2::geom_point(ggplot2::aes(age, p_obs_bin, size = bin_size), fill = "#7a0177", colour = "black", shape = 21) +
        ggplot2::theme_bw(size_text) +
        ggplot2::coord_cartesian(xlim = c(0, 60), ylim = c(0,1)) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::ylab("Sero-positivity") + ggplot2::xlab("Age")

    } } else

    {
      print("model did not run")
      print_warning <- "errors"
      df <- data.frame()

      prev_plot <-  ggplot2::ggplot(df) + ggplot2::geom_point() + ggplot2::xlim(0, 10) + ggplot2::ylim(0, 10) +
                    ggplot2::annotate("text", x = 4, y = 5, label = print_warning) +
                    ggplot2::theme_bw(25) +
                    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank()) +
                    ggplot2::ylab(" ") + ggplot2::xlab(" ")
                    ggplot2::theme(plot.title = ggplot2::element_text(size = 10))

    }

  return(prev_plot)

}

#' Generate Force-of-Infection Plot
#'
#' Function that generates the force of infection plot
#' @param model_object Object that the run_model function returns with the results of the fit
#' @param model Refers to model selected
#' @param xlabel Label of axis x
#' @param ylabel Label of axis y
#' @return Force of infection plot
#' @export
plot_foi <- function(model_object,
                     lambda_sim = NA,
                     max_lambda = NA,
                     size_text = 25) {

  if (is.character(model_object$fit) == FALSE)  {
    if  (class(model_object$fit@sim$samples)  != "NULL" ) {

      foi <- rstan::extract(model_object$fit,
                            "foi",
                            inc_warmup = FALSE)[[1]]

      #-------- This bit is to get the actual length of the foi data
      foi_data <- model_object$foi_cent_est

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
        ggplot2::geom_ribbon(ggplot2::aes(x = year, ymin = lower, ymax = upper), fill = "#41b6c4", alpha = 0.5) +
        ggplot2::geom_line(ggplot2::aes(x = year, y = medianv), colour = "#253494", size = size_text/8) +
        ggplot2::theme_bw(size_text) +
        ggplot2::coord_cartesian(ylim = c(0, max_lambda)) +
        ggplot2::ylab("Force-of-Infection") + ggplot2::xlab("Year")
      }

    } else

    {
      print("model did not run")
      print_warning <- "errors"
      df <- data.frame()

      foi_plot <- ggplot2::ggplot(df) + ggplot2::geom_point() + ggplot2::xlim(0, 10) + ggplot2::ylim(0, 10) +
                  ggplot2::annotate("text", x = 4, y = 5, label = print_warning) +
                  ggplot2::theme_bw(25) +
                  ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank()) +
                  ggplot2::ylab(" ") + ggplot2::xlab(" ")
                  ggplot2::theme(plot.title = ggplot2::element_text(size = 10))

    }

  return(foi_plot)

}

#' Generate Rhats-Convergence Plot
#'
#' Function that generates the convergence graph of a model
#' @param model_object Object that the run_model function returns with the results of the fit
#' @param model Refers to model selected
#' @param xlabel Label of axis x
#' @param ylabel Label of axis y
#' @return The rhats-convergence plot of the selected model
#' @export
plot_rhats <- function(model_object,
                       size_text = 25) {

  if (is.character(model_object$fit) == FALSE)  {
    if  (class(model_object$fit@sim$samples)  != "NULL" ) {

      rhats <- get_table_rhats(model_object)

      rhats_plot <- ggplot2::ggplot(rhats, ggplot2::aes(year, rhat)) +
                    ggplot2::geom_line(colour = "purple") +
                    ggplot2::geom_point() +
                    ggplot2::coord_cartesian(ylim = c(0.7, 2)) +
                    ggplot2::geom_hline(yintercept = 1.1, colour = "blue", size = size_text/12) +
                    ggplot2::theme_bw(size_text) +
                    ggplot2::ylab("Convergence (R^)")

    } } else

    {
      print("model did not run")
      print_warning <- "errors"
      df <- data.frame()

      rhats_plot <- ggplot2::ggplot(df) + ggplot2::geom_point() + ggplot2::xlim(0, 10) + ggplot2::ylim(0, 10) +
        ggplot2::annotate("text", x = 4, y = 5, label = print_warning) +
        ggplot2::theme_bw(25) +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank()) +
        ggplot2::ylab(" ") + ggplot2::xlab(" ")
      ggplot2::theme(plot.title = ggplot2::element_text(size = 10))

    }

  return(rhats_plot)

}


#' Generate a vertical arrange of plots summarizing the results of the model implementation
#'
#' Function that generates the combined graph
#' @param model_object Object that the run_model function returns with the results of the fit
#' @param model Refers to model selected
#' @param xlabel Label of axis x
#' @param ylabel Label of axis y
#' @return The combined plots
#' @export
plot_model <- function(model_object,
                       lambda_sim = NA,
                       max_lambda = NA,
                       size_text = 25) {

  if (is.character(model_object$fit) == FALSE)  {
    if  (class(model_object$fit@sim$samples)  != "NULL" ) {
      prev_plot <- plot_seroprev(model_object = model_object,
                                 size_text = size_text)

      foi_plot <- plot_foi(model_object = model_object,
                          lambda_sim = lambda_sim,
                          max_lambda = max_lambda,
                          size_text = size_text)

      rhats_plot <- plot_rhats(model_object = model_object,
                              size_text = size_text)

      summary_plot <- plot_info_table(t(model_object$model_summary), size_text = size_text)

      plot_arrange <- gridExtra::grid.arrange(summary_plot,
                                              prev_plot,
                                              foi_plot,
                                              rhats_plot, nrow = 4,
                                              heights = c(1.5, 1, 1, 1))

    } } else

    {
      print("model did not run")
      print_warning <- "errors"
      df <- data.frame()

      g0 <- ggplot2::ggplot(df) + ggplot2::geom_point() + ggplot2::xlim(0, 10) + ggplot2::ylim(0, 10) +
        ggplot2::annotate("text", x = 4, y = 5, label = print_warning) +
        ggplot2::theme_bw(25) +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank()) +
        ggplot2::ylab(" ") + ggplot2::xlab(" ")
      g1 <- g0
      g0 <- g0 + ggplot2::labs(subtitle = model_object$model) +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 10))

      plot_arrange <- gridExtra::grid.arrange(g0, g1, g1, g1, g1, nrow = 5)

    }

  return(plot_arrange)

}


#' Plot Info Table
#'
#' Function that generates the information table
#' @param model_data refers to data of the each model
#' @param info the information that will be contained in the table
#' @param size_text text size
#' @return The previous expanded graphic
#' @export
plot_info_table <- function(info, size_text){

  dato <- data.frame(
    y = NROW(info):1,
    text = paste0(rownames(info), ": ",info[,1])
  )
  p <- ggplot2::ggplot(dato, ggplot2::aes(x = 1, y = y)) +
    ggplot2::scale_y_continuous(limits = c(0, NROW(info) + 1), breaks = NULL) +
    # scale_x_continuous(breaks=NULL) +
    ggplot2::theme_void() +
    ggplot2::geom_text(ggplot2::aes(label = text), size = size_text/2.5, fontface = "bold")

  return(p)
}
