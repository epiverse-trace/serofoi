#' Generate Combined Plots
#'
#' Función que genera la gráfica combinada
#' Function that generates the combined graph
#' @param res
#' @param data Data
#' @param xlabel Label of axis x
#' @param ylabel Label of axis y
#' @return The combined plots
#' @export
generate_combined_plots <- function(res, dat, lambda_sim = NA, max_lambda = NA, size_text = 25) {

  if (is.character(res$fit) == FALSE)  {
    if  (class(res$fit@sim$samples)  != "NULL" ) {

      summary_mod <- extract_summary_mod(res, dat)
      foi <- rstan::extract(res$fit, 'foi', inc_warmup = FALSE)[[1]]


      prev_expanded <- get_prev_expanded(foi, dat)
      plot_prev <-
        ggplot2::ggplot(prev_expanded) +
        ggplot2::geom_ribbon(ggplot2::aes(x = age, ymin = predicted_prev_lower, ymax = predicted_prev_upper), fill = '#c994c7') +
        ggplot2::geom_line(ggplot2::aes(x = age, y = predicted_prev), colour = '#7a0177') +
        ggplot2::geom_errorbar(ggplot2::aes(age, ymin = p_obs_bin_l, ymax = p_obs_bin_u), width = 0.1) +
        ggplot2::geom_point(ggplot2::aes(age, p_obs_bin, size = bin_size), fill = '#7a0177', colour = 'black', shape = 21) +
        ggplot2::theme_bw(size_text) +
        ggplot2::coord_cartesian(xlim = c(0, 60), ylim = c(0,1)) +
        ggplot2::theme(legend.position = 'none') +
        ggplot2::ylab('Sero-positivity') + ggplot2::xlab("Age")



      #-------- This bit is to get the actual length of the foi data
      foi_dat <- res$foi_cent_est
      if (!is.na(lambda_sim)) {
        lambda_mod_length <- NROW(foi_dat)
        lambda_sim_length <- length(lambda_sim)

        if (lambda_mod_length < lambda_sim_length) {
          remove_x_values <- lambda_sim_length - lambda_mod_length
          lambda_sim <- lambda_sim[-c(1:remove_x_values)]
        }

        foi_dat$simulated <- lambda_sim

      }

      #--------
      foi_dat$medianv[1] <- NA
      foi_dat$lower[1] <- NA
      foi_dat$upper[1] <- NA


      plot_foi <-
        ggplot2::ggplot(foi_dat) +
        ggplot2::geom_ribbon(ggplot2::aes(x = year, ymin = lower, ymax = upper), fill = '#41b6c4', alpha = 0.5) +
        ggplot2::geom_line(ggplot2::aes(x = year, y = medianv), colour = '#253494', size = size_text/8) +
        ggplot2::theme_bw(size_text) +
        ggplot2::coord_cartesian(ylim = c(0, max_lambda)) +
        ggplot2::ylab('Force-of-Infection') + ggplot2::xlab("Year")



      if (!is.na(lambda_sim)) {
        lambda_plot <- plot_foi + ggplot2::geom_line(ggplot2::aes(x = year, y = simulated), colour = 'red', size = 1.5)
      }

      rhats <- get_table_rhats(res)

      plot_rhats <- ggplot2::ggplot(rhats, ggplot2::aes(year, rhat)) +
        ggplot2::geom_line(colour = 'purple') +
        ggplot2::geom_point() +
        ggplot2::coord_cartesian(ylim = c(0.7, 2)) +
        ggplot2::geom_hline(yintercept = 1.1, colour = 'blue', size = size_text/12) +
        ggplot2::theme_bw(size_text) +
        ggplot2::ylab('Convergence (R^)')


      plot_data <- plot_info_table(t(summary_mod), size_text = size_text)


      plot_arrange <- grid.arrange(plot_data,
                                   plot_prev,
                                   plot_foi,
                                   plot_rhats, nrow = 4,
                                   heights = c(1.5, 1, 1, 1))

      res_plots <- list(plots = list(plot_data = plot_data,
                                      plot_prev = plot_prev,
                                      plot_foi  = plot_foi,
                                      plot_rhats = plot_rhats),
                         summary_mod     = summary_mod,
                         rhats           = rhats,
                         prev_expanded = prev_expanded)



    } } else

    {
      print('model did not run')
      print_warning <- 'errors'
      df <- data.frame()

      g0 <- ggplot2::ggplot(df) + ggplot2::geom_point() + ggplot2::xlim(0, 10) + ggplot2::ylim(0, 10) +
        ggplot2::annotate("text", x = 4, y = 5, label = print_warning) +
        ggplot2::theme_bw(25) +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank()) +
        ggplot2::ylab('') + ggplot2::xlab('')
      g1 <- g0
      g0 <- g0 + ggplot2::labs(subtitle = res$model) +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 10))


      plots <- grid.arrange(g0, g1, g1, g1, g1, nrow = 5)
      res_p <- list(plots = plots,
                     loo_fit = c('Not available'),
                     summary_mod     = c('Not available'),
                     prev_expanded = c('Not available')
      )

      res_plots <-  list(plots = plots,
                          summary_mod  =   "model did not run")


    }



  return(res_plots)


}

#' Get Model Comparison Plot
#'
#' Función que obtiene la gráfica de comparación de modelos
#' Function that obtains the model comparison plot
#' @param res_comp
#' @param xlabel Label of axis x
#' @param ylabel Label of axis y
#' @return The model comparison graph

#========================= PLOT comparison
get_model_comparison_plot <- function(res_comp) {
  model_comp  <- res_comp$model_comp
  best_model  <- as.character(res_comp$best_model_data1$model)
  best_modelP <- as.numeric(res_comp$best_model_data1$pvalue)

  emptyp <-  ggplot2::ggplot(data = data.frame()) +
    ggplot2::geom_point() +
    ggplot2::xlim(0,1) +  ggplot2::ylim(0,1) +  ggplot2::theme_void() +
    ggplot2::annotate('text',  x = .5, y = .6, label = best_model, size = 14) +
    ggplot2::annotate('text',  x = .5, y = .55, label = '(best model)', size = 13)


  infot <- dplyr::filter(model_comp, converged == 'Yes') %>% dplyr::select(model, difference, diff_se, pvalue, best) %>%
    dplyr::mutate(diff = round(difference, 2),
           diff_se = round(diff_se, 2),
           pvalue = round(pvalue, 4))

  blank <- data.frame(x = 1:10, y = 1:100)
  table_pars <-
    ggplot2::ggplot(blank, ggplot2::aes(x, y)) +
    ggplot2::geom_blank() + ggplot2::ylab('') + ggplot2::xlab('') +
    ggplot2::annotation_custom(gridExtra::tableGrob(d = infot,
                                theme = gridExtra::ttheme_default(base_size = 20))) +
    ggplot2::theme_void(30)

  pf <- cowplot::plot_grid(emptyp, table_pars, nrow = 2 )

  return(pf)

}

#' Get Vertical Plot Arrange per Model
#'
#' Función que genera el grafico en un arreglo vertical por modelo
#' Function that generates the graph in a vertical arrange per model
#' @param PPC
#' @param xlabel Label of axis x
#' @param ylabel Label of axis y
#' @return The vertical plot arrange per model
#' @export

vertical_plot_arrange_per_model <- function(PPC){

  pp <- gridExtra::grid.arrange(PPC$plots$plot_data,
                     PPC$plots$plot_prev,
                     PPC$plots$plot_foi,
                     PPC$plots$plot_rhats,
                     nrow = 4,
                     heights = c(1.5, 1, 1, 1))
  return(pp)

}

#' Generate Info Table Plot
#'
#' Función que genera la tabla de información
#' Function that generates the information table
#' @param data Data
#' @param info
#' @param size_text
#' @return The previous expanded graphic
#' @export

plot_info_table <- function(info, size_text){

  dato <- data.frame(
    y = NROW(info):1,
    text = paste0(rownames(info), ': ',info[,1])
  )
  p <- ggplot2::ggplot(dato, ggplot2::aes(x = 1, y = y)) +
    ggplot2::scale_y_continuous(limits = c(0, NROW(info) + 1), breaks = NULL) +
    # scale_x_continuous(breaks=NULL) +
    ggplot2::theme_void() +
    ggplot2::geom_text(ggplot2::aes(label = text), size = size_text/2.5, fontface = 'bold')

  return(p)
}