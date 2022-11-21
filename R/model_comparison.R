#' Fit Model Log
#'
#' Función que ajusta el modelo logarítmico a los datos
#' Function that fits the logarithmic model to the data
#' @param data dat
#' @param model
#' @param m_name
#' @param n_iters
#' @param n_thin
#' @param delta
#' @param mtreed
#' @param mDecades
#' @return res
#' @export
compare_and_save_best_model <- function(survey,
                                         model_comparison,
                                         name_fitting_res,
                                         mod_0, mod_1, mod_2)
{


  # --------------- Model comparison
  model_comp <- model_comparison %>% dplyr::select(-.data$performance)
  model_comp <- model_comp %>% dplyr::mutate(pvalue = 1 - stats::pnorm(.data$difference/.data$diff_se,0,1))
  model_comp$pvalue[is.nan(model_comp$pvalue)] <- 1
  model_comp$pvalue <- model_comp$pvalue * stats::runif(NROW(model_comp), min = 1, max = 1.0001)# I make this just to ensure I get different values

  model_comp$better <- NA
  model_comp$better[model_comp$difference > 0] <- 'Yes'
  model_comp$better[model_comp$difference <= 0] <- 'No'

  model_comp$better[model_comp$model == 'Constant'] <- "-"
  model_comp$pvalue[model_comp$model == 'Constant'] <- 0


  model_comp$converged[model_comp$elpd == -1.000e+10] <- 'No'


  ds_one <- dplyr::filter(model_comp, .data$converged == 'Yes')
  print(paste0('number of converged models = ', NROW(ds_one)))


  elps_order <-  rev(sort(ds_one$elpd)) #[1:3]
  best <- dplyr::filter(model_comp, .data$elpd %in% elps_order) %>% dplyr::arrange(-.data$elpd)# This is to make sure I keep only three
  best_model1 <- as.character(best$model[1])
  best_model2 <- as.character(best$model[2])
  best_model3 <- as.character(best$model[3])



  model_comp$best <- NA
  model_comp$best[model_comp$model == best_model1] <- 1
  model_comp$best[model_comp$model == best_model2] <- 2
  model_comp$best[model_comp$model == best_model3] <- 3

  model_comp <- model_comp %>% dplyr::arrange(best)
  model_comp$pvalue <- round(model_comp$pvalue, 6)

  # --------------- Best model
  best_model_data1 <- dplyr::filter(model_comp, best == 1) # Here I choose the maximun difference rather than the lowest p value
  best_model_data2 <- dplyr::filter(model_comp, best == 2) # Here I choose the maximun difference rather than the lowest p value
  best_model_data3 <- dplyr::filter(model_comp, best == 3) # Here I choose the maximun difference rather than the lowest p value

  RealYexpo  <-  mod_0$RealYexpo
  best_model_1 <- as.character(best_model_data1$model)
  best_model_2 <- as.character(best_model_data2$model)
  best_model_3 <- as.character(best_model_data3$model)

  if (best_model_1 == mod_0$model) {
    res_file_1 <- mod_0
  }

  if (best_model_1 ==  mod_1$model) {
    res_file_1 <- mod_1
  }

  if (best_model_1 ==  mod_2$model) {
    res_file_1 <- mod_2
  }

  # ------- Best 2
  if (best_model_2  == mod_0$model) {
    res_file_2 <- mod_0
  }

  if (best_model_2 ==  mod_1$model) {
    res_file_2 <- mod_1
  }

  if (best_model_2 ==  mod_2$model) {
    res_file_2 <- mod_2
  }

  # ------- Best 3
  if (best_model_3  == mod_0$model) {
    res_file_3 <- mod_0
  }

  if (best_model_3 ==  mod_1$model) {
    res_file_3 <- mod_1
  }

  if (best_model_3 ==  mod_2$model) {
    res_file_3 <- mod_2
  }

  # --- save_best_model
  extract_and_save(res_file_1, res_file_2, res_file_3,
                   best_model_1, best_model_2, best_model_3,
                   name_fitting_res,
                   survey, RealYexpo)

  res_comp <- list(best_model_data1 = best_model_data1,
                    best_model_data2 = best_model_data2,
                    model_comp = model_comp)
  return(res_comp)

}

#' Extract and Save
#'
#' Función que extrae y guarda resultados
#' Function that extracts and saves results
#' @param survey
#' @param res res_file_1
#' @param res res_file_2
#' @param res res_file_3
#' @param RealYexpo
#' @param best_model best_model_1
#' @param best_model best_model_2
#' @param best_model best_model_3
#' @return extracted and saved results
#' @export
extract_and_save <- function(res_file_1, res_file_2, res_file_3,
                             best_model_1, best_model_2, best_model_3,
                             name_file,
                             survey, RealYexpo) {

  foi_0 <- rstan::extract(res_file_1$fit, 'foi', inc_warmup = FALSE)[[1]]
  foi_1 <- rstan::extract(res_file_2$fit, 'foi', inc_warmup = FALSE)[[1]]
  foi_2 <- rstan::extract(res_file_3$fit, 'foi', inc_warmup = FALSE)[[1]]


  foi_cent_est1 <- data.frame(year  = RealYexpo,
                              lower = apply(foi_0, 2, function(x) quantile(x, 0.1)),
                              upper = apply(foi_0, 2, function(x) quantile(x, 0.9)),
                              median = apply(foi_0, 2, function(x) quantile(x, 0.5))) %>%
    dplyr::mutate(best = 'best1', name_model = best_model_1)

  foi_cent_est2 <- data.frame(year  = RealYexpo,
                              lower = apply(foi_1, 2, function(x) quantile(x, 0.1)),
                              upper = apply(foi_1, 2, function(x) quantile(x, 0.9)),
                              median = apply(foi_1, 2, function(x) quantile(x, 0.5))) %>%
    dplyr::mutate(best = 'best2', name_model = best_model_2)

  foi_cent_est3 <- data.frame(year  = RealYexpo,
                              lower = apply(foi_1, 2, function(x) quantile(x, 0.1)),
                              upper = apply(foi_1, 2, function(x) quantile(x, 0.9)),
                              median = apply(foi_1, 2, function(x) quantile(x, 0.5))) %>%
    dplyr::mutate(best = 'best3', name_model = best_model_3)


  foi_cent_est <- rbind(foi_cent_est1, foi_cent_est2, foi_cent_est3)


  foi_0_post_1000s <- dplyr::sample_n(as.data.frame(foi_0), size = 1000) %>% dplyr::mutate(best = 'best1', name_model = best_model_1)
  foi_1_post_1000s <- dplyr::sample_n(as.data.frame(foi_1), size = 1000) %>% dplyr::mutate(best = 'best2', name_model = best_model_2)
  foi_2_post_1000s <- dplyr::sample_n(as.data.frame(foi_2), size = 1000) %>% dplyr::mutate(best = 'best3', name_model = best_model_3)

  foi_post_1000s <- rbind(foi_0_post_1000s, foi_1_post_1000s, foi_2_post_1000s)

  colnames(foi_post_1000s)[1:length(RealYexpo)] <- RealYexpo


  prev1 <- res_file_1$prev_expanded %>% dplyr::mutate(best = 'best1', name_model = best_model_1)
  prev2 <- res_file_2$prev_expanded %>% dplyr::mutate(best = 'best2', name_model = best_model_2)
  prev3 <- res_file_3$prev_expanded %>% dplyr::mutate(best = 'best3', name_model = best_model_3)


  prevalence_expanded <- rbind(prev1, prev2, prev3)


  fres <- list(dataset = survey,
               foi_cent_est = foi_cent_est,
               foi_post_1000s = foi_post_1000s,
               prevalence    = prevalence_expanded)

  saveRDS(fres, name_file)

}

#' Extract Summary Mod
#'
#' Función que hace un resumen de los modelos
#' Function that summarizes the models
#' @param res
#' @param data data
#' @return summary of the models
#' @export
extract_summary_mod <- function(res, dat) {

  model_name <- res$model
  #------- Loo estimates

  loo_fit <- res$loo_fit
  if (sum(is.na(loo_fit)) < 1)
  {
    lll <- as.numeric((round(loo_fit$estimates[1,],2)))} else
    {
      lll <- c(-1e10, 0)
    }



  summary_mod <- data.frame(model = res$model,
                            dataset = dat$survey[1],
                            country = dat$country[1],
                            year    = dat$tsur[1],
                            test    = dat$test[1],
                            antibody = dat$antibody[1],
                            n_sample = sum(dat$total),
                            n_agec  = length(dat$age_mean_f),
                            n_iter  = res$n_iters,
                            performance = "_____",
                            elpd = lll[1],
                            se = lll[2],
                            converged = NA
  )

  rhats <- get_table_rhats(res)
  if (any(rhats$rhat > 1.1 ) == FALSE) {
    summary_mod$converged = "Yes"  }


  return(summary_mod)
}

#' Get Table Rhats
#'
#' Función que hace la tabla de los rhats
#' Function that makes the rhats table
#' @param res
#' @return rhats table
#' @export
get_table_rhats <- function(res) {

  rhats <- bayesplot::rhat(res$fit, "foi")

  if (any(is.nan(rhats))) {
    rhats[which(is.nan(rhats))] <- 0}

  res_rhats <- data.frame(year = res$RealYexpo, rhat = rhats)
  res_rhats$rhat[res_rhats$rhat == 0] <- NA # This is because I'm not estimating these foi values

  return(res_rhats)
}
