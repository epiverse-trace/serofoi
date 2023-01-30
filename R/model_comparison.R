#' Get Table Rhats
#'
#' Function that makes the rhats table
#' @param model_object model_object
#' @return rhats table
#' @export
get_table_rhats <- function(model_object) {

  rhats <- bayesplot::rhat(model_object$fit, "foi")

  if (any(is.nan(rhats))) {
    rhats[which(is.nan(rhats))] <- 0}
  model_rhats <- data.frame(year = model_object$real_yexpo, rhat = rhats)
  model_rhats$rhat[model_rhats$rhat == 0] <- NA # This is because I'm not estimating these foi values

  return(model_rhats)
}


#' Get Model Table Comparison
#' Provides a table with statistics for comparison between models and selection
#' @param model_objects_list model_objects to compare
#' @return comparison table
#' @export
get_comparison_table <- function(model_objects_list) {


  dif_m0_m1 <- loo::loo_compare(model_objects_list$m0.loo_fit,
                                 model_objects_list$m1.loo_fit)

  dif_m0_m2 <- loo::loo_compare(model_objects_list$m0.loo_fit,
                                 model_objects_list$m2.loo_fit)

  # Aquí pendiente revisar <diference> desde la función summary_model
  # No estoy segura que este parámetro venga bien desde allá ni tampoco que esté bien acá

  # model_comp$better <- NA
  # model_comp$better[model_comp$difference > 0] <- 'Yes'
  # model_comp$better[model_comp$difference <= 0] <-'No'
  # model_comp$better[model_comp$model == 'constant_foi_bi'] <- "-"


  model_objects_list$m0.model_summary$difference <- 0
  model_objects_list$m0.model_summary$diff_se <- 1

  model_objects_list$m1.model_summary$difference <- dif_m0_m1[1]
  model_objects_list$m1.model_summary$diff_se <- dif_m0_m1[2]

  model_objects_list$m2.model_summary$difference <- dif_m0_m2[1]
  model_objects_list$m2.model_summary$diff_se <- dif_m0_m2[2]

  model_comp <- rbind(model_objects_list$m0.model_summary,
                      model_objects_list$m1.model_summary,
                      model_objects_list$m2.model_summary)

  model_comp$converged[model_comp$elpd == -1.000e+10] <- 'No' #
  ds_one <- dplyr::filter(model_comp, converged == 'Yes')
  print(paste0('number of converged models = ', NROW(ds_one)))

  # Ordering the best model based on elpd values
  elps_order <-  rev(sort(ds_one$elpd))
  best <- dplyr::filter(model_comp, elpd %in% elps_order) %>% dplyr::arrange(-.data$elpd)# This is to make sure I keep only three
  best_model1 <- as.character(best$model[1])
  best_model2 <- as.character(best$model[2])
  best_model3 <- as.character(best$model[3])

  model_comp$best_elpd <- NA
  model_comp$best_elpd[model_comp$model == best_model1] <- 1
  model_comp$best_elpd[model_comp$model == best_model2] <- 2
  model_comp$best_elpd[model_comp$model == best_model3] <- 3
  model_comp <- model_comp %>% dplyr::arrange(.data$best_elpd)

  # Estimating p-values to check the difference between the models m0 and other models is actually important
  model_comp <- model_comp %>% dplyr::mutate(pvalue = 1 - stats::pnorm(difference/diff_se,0,1))
  model_comp$pvalue[is.nan(model_comp$pvalue)] <- 1
  model_comp$pvalue <- model_comp$pvalue * stats::runif(NROW(model_comp), min = 1, max = 1.0001)# I make this just to ensure I get different values
  model_comp$pvalue[model_comp$model == 'constant_foi_bi'] <- 0
  model_comp$pvalue <- round(model_comp$pvalue, 6)

  return(model_comp)
}
