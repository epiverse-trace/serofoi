#' Run Model 0
#'
#' Función que corre el modelo 0
#' Function that runs model 0
#' @param my_dir
#' @param survey
#' @param data dat0
#' @param n_iters
#' @param n_warmup
#' @param model model_0
#' @param name_model name_model_0
#' @return results of model 0
#' @export
run_model_0 <- function(my_dir,
                        survey,
                        dat0,
                        n_iters,
                        n_warmup,
                        model_0, name_model_0) {

  t0 <- Sys.time()
  my_dir0 <- paste0(config::get("test_files_path"), my_dir)


  data <- dplyr::filter(dat0, .data$survey == survey) %>% dplyr::arrange(.data$age_mean_f) %>%
    dplyr::mutate(birth_year = .data$tsur - .data$age_mean_f)

  model_0 <- fit_model(model = model_0, data,
                       m_name = name_model_0, n_iters = n_iters); print(paste0(survey, "finished ------ model_0"))


  #  ---- Plotting
  foi_model <- rstan::extract(model_0$fit, "foi", inc_warmup = FALSE)[[1]]
  max_lambda <-  (as.numeric(quantile(foi_model, 0.95))) * 1.3

  plots_model_0    <- generate_combined_plots(result = model_0, data = data, lambda_sim = NA, max_lambda) ; print(paste0(survey, "finished ------ plots_model_0"))


  plots_model_0$summary_model$difference <- 0; plots_model_0$summary_model$diff_se <- 1;

  model_comparison <- rbind(plots_model_0$summary_model)


  model_0$prev_expanded <- plots_model_0$prev_expanded


  # ---------------------
  name_fitting_plots  <- paste0(my_dir0, "/fitting_plots/", survey, ".png")
  name_fitting_result  <- paste0(my_dir0, "/fitting_results/",survey, ".RDS")

  # browser()
  result_comp <- compare_and_save_best_model(survey = survey,
                                             model_comparison = model_comparison,
                                             name_fitting_result = name_fitting_result,
                                             model_0 = model_0)


  model_comp    <- result_comp$model_comp
  model_comp_plot <- get_model_comparison_plot(result_comp)


  parrange_0 <- vertical_plot_arrange_per_model(plots_model_0)


  grDevices::png(name_fitting_plots, width = 2000, height = 2000)

  grid.arrange(parrange_0,
               nrow = 1)
  dev.off()

  print(paste("end", survey))


  t5 <- Sys.time()
  time_taken <- t5 - t0
  print(time_taken)



  result_survey <- list(data = data,
                        time_taken = time_taken,
                        model_0 = model_0,
                        model_comp = model_comp,
                        plots_model_0 = plots_model_0)



  saveRDS(result_survey, name_fitting_result)


}

#' Dir Results
#'
#' Función guarda los resultados en una dirección
#' Function saves the results in a place
#' @param name_dir
#' @return dir
#' @export
dir_results <- function(name_dir)
{

  my_dir <- paste0("tests/", name_dir)
  dir_plots <- paste0(my_dir, "/fitting_plots")
  dir_posterior  <- paste0(my_dir, "/fitting_results")

  if (dir.exists(my_dir) == FALSE) {
    dir.create(my_dir)
  }

  if (dir.exists(dir_plots) == FALSE) {
    dir.create(dir_plots)
  }

  if (dir.exists(dir_posterior) == FALSE) {
    dir.create(dir_posterior)
  }

  print(paste0("NOTE: My results will be stored at:_________/", my_dir))

}


#. --- Function for cleaning seroprevalence data
# conf <- data.frame(Hmisc::binconf(dat$counts, dat$total,method="exact"))
# dat  <- cbind(dat, conf) %>%
#   rename (prev_obs = PointEst,
#           prev_obs_lower = Lower,
#           prev_obs_upper = Upper
#   )



# For other models

#,
#model_1, name_model_1,
#model_2, name_model_2


#model_1 <- fit_model(model = model_1, data,
#m_name = name_model_1, n_iters = n_iters); print(paste0(survey, "finished ------ model_1"))

#model_2 <- fit_model_log(model = model_2, data,
#m_name = name_model_2, n_iters = n_iters); print(paste0(survey, "finished ------ model_2"))

#Plotting

#plots_model_1    <- generate_combined_plots(result = model_1, data = data, lambda_sim = NA, max_lambda) ; print(paste0(survey, "finished ------ plots_model_1"))
#plots_model_2    <- generate_combined_plots(result = model_2, data = data, lambda_sim = NA, max_lambda) ; print(paste0(survey, "finished ------ plots_model_2"))

# ---- Statistical Checking
#dif_m0_m1 <- loo::loo_compare(model_0$loo_fit, model_1$loo_fit)
#dif_m0_m2 <- loo::loo_compare(model_0$loo_fit, model_2$loo_fit)

#plots_model_1$summary_model$difference <- dif_m0_m1[1];   plots_model_1$summary_model$diff_se <- dif_m0_m1[2];
#plots_model_2$summary_model$difference <- dif_m0_m2[1];   plots_model_2$summary_model$diff_se <- dif_m0_m2[2];

#,
#plots_model_1$summary_model, plots_model_2$summary_model

#model_1$prev_expanded <- plots_model_1$prev_expanded
#model_2$prev_expanded <- plots_model_2$prev_expanded

#,
#model_1 = model_1,
#model_2 = model_2

#parrange_1 <- vertical_plot_arrange_per_model(plots_model_1)
#parrange_2 <- vertical_plot_arrange_per_model(plots_model_2)

#parrange_1,
#parrange_2,

#model_1 = model_1,
#model_2 = model_2,

#,
#plots_model_1 = plots_model_1,
#plots_model_2 = plots_model_2

