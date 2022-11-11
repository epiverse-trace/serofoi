run_save_models <- function(my_dir,
                          suv,
                          dat0,
                          n_iters,
                          n_warmup,
                          Model0, NameModel0,
                          Model1, NameModel1,
                          Model2, NameModel2) {

  t0 <- Sys.time()
  my_dir0 <- paste0(config::get("test_files_path"), my_dir)


  dat <- filter(dat0, survey == suv) %>% arrange(age_mean_f) %>%
    mutate(birth_year = tsur - age_mean_f)

  mod_0 <- fit_model(model = Model0, dat,
                     m_name = NameModel0, n_iters = n_iters); print(paste0(suv, ' finished ------ Model_0'))

  mod_1 <- fit_model(model = Model1, dat,
                     m_name = NameModel1, n_iters = n_iters); print(paste0(suv, ' finished ------ Model_1'))

  mod_2 <- fit_model_log(model = Model2, dat,
                         m_name = NameModel2, n_iters = n_iters); print(paste0(suv, ' finished ------ Model_2'))


  #  ---- Plotting
  foi_mod <- rstan::extract(mod_2$fit, 'foi', inc_warmup = FALSE)[[1]]
  max_lambda <-  (as.numeric(quantile(foi_mod, 0.95))) * 1.3

  PlotsM0    <- generate_combined_plots(res = mod_0, dat = dat, lambda_sim = NA, max_lambda) ; print(paste0(suv, ' finished ------ Plots_M0'))
  PlotsM1    <- generate_combined_plots(res = mod_1, dat = dat, lambda_sim = NA, max_lambda) ; print(paste0(suv, ' finished ------ Plots_M1'))
  PlotsM2    <- generate_combined_plots(res = mod_2, dat = dat, lambda_sim = NA, max_lambda) ; print(paste0(suv, ' finished ------ Plots_M2'))

  # ---- Statistical Checking
  dif_m0_m1 <- loo_compare (mod_0$loo_fit, mod_1$loo_fit)
  dif_m0_m2 <- loo_compare (mod_0$loo_fit, mod_2$loo_fit)

  PlotsM0$summary_mod$difference <- 0; PlotsM0$summary_mod$diff_se <- 1;
  PlotsM1$summary_mod$difference <- dif_m0_m1[1];   PlotsM1$summary_mod$diff_se <- dif_m0_m1[2];
  PlotsM2$summary_mod$difference <- dif_m0_m2[1];   PlotsM2$summary_mod$diff_se <- dif_m0_m2[2];

  model_comparison <- rbind(PlotsM0$summary_mod,
                            PlotsM1$summary_mod, PlotsM2$summary_mod)


  mod_0$prev_expanded <- PlotsM0$prev_expanded
  mod_1$prev_expanded <- PlotsM1$prev_expanded
  mod_2$prev_expanded <- PlotsM2$prev_expanded

  # ---------------------
  name_fitting_plots  <- paste0(my_dir0, '/fitting_plots/', suv, '.png')
  name_fitting_res  <- paste0(my_dir0, '/fitting_results/',suv, '.RDS')

  # browser()
  res_comp <- compare_and_save_best_model(survey = suv,
                                          model_comparison = model_comparison,
                                          name_fitting_res = name_fitting_res,
                                          mod_0 = mod_0,
                                          mod_1 = mod_1,
                                          mod_2 = mod_2)

  model_comp    <- res_comp$model_comp
  mod_comp_plot <- get_model_comparison_plot(res_comp)


  parrange_0 <- vertical_plot_arrange_per_model(PlotsM0)
  parrange_1 <- vertical_plot_arrange_per_model(PlotsM1)
  parrange_2 <- vertical_plot_arrange_per_model(PlotsM2)




  png(name_fitting_plots, width = 2000, height = 2000)

  grid.arrange(parrange_0,
               parrange_1,
               parrange_2,
               nrow = 1)
  dev.off()

  print(paste('end', suv))


  t5 <- Sys.time()
  time_taken <- t5 - t0
  print(time_taken)



  res_survey <- list(dat = dat,
                     time_taken = time_taken,
                     mod_0 = mod_0,
                     mod_1 = mod_1,
                     mod_2 = mod_2,
                     model_comp = model_comp,
                     PlotsM0 = PlotsM0,
                     PlotsM1 = PlotsM1,
                     PlotsM2 = PlotsM2)



  saveRDS(res_survey, name_fitting_res)


}

dir_results <- function(name_dir)
{

  my_dir <- paste0('tests/', name_dir)
  dir_plots <- paste0(my_dir, '/fitting_plots')
  dir_posterior  <- paste0(my_dir, '/fitting_results')

  if (dir.exists(my_dir) == FALSE) {
    dir.create(my_dir)
  }

  if (dir.exists(dir_plots) == FALSE) {
    dir.create(dir_plots)
  }

  if (dir.exists(dir_posterior) == FALSE) {
    dir.create(dir_posterior)
  }

  print(paste0("NOTE: My results will be sortored at:_________/", my_dir))

}
