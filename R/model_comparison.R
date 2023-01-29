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
#' Provides a table with statistics for comparison between models
#' @param model_objects_list model_objects to compare within a list
#' @return comparison table
#' @export
get_comparison_table <- function(model_objects_list) {

  dif_m0_m1 <- loo::loo_compare (model_objects_list$m0.loo_fit,
                            model_objects_list$m1.loo_fit)

  dif_m0_m2 <- loo::loo_compare (model_objects_list$m0.loo_fit,
                                 model_objects_list$m2.loo_fit)

  model_objects_list$m0.model_summary$difference <- 0
  model_objects_list$m0.model_summary$diff_se <- 1

  model_objects_list$m1.model_summary$difference <- dif_m0_m1[1]
  model_objects_list$m1.model_summary$diff_se <- dif_m0_m1[2]

  model_objects_list$m2.model_summary$difference <- dif_m0_m2[1]
  model_objects_list$m2.model_summary$diff_se <- dif_m0_m2[2]

  comparison_table <- rbind(model_objects_list$m0.model_summary,
                            model_objects_list$m1.model_summary,
                            model_objects_list$m2.model_summary)


  return(comparison_table)
}


#' Select best model
#' Frovides an authomated selection of best model based on statistic parameters
#' @param model_objects_list model_objects to compare
#' @return compasiron table
#' @export
get_comparison_table <- function(model_objects_list) {

  dif_m0_m1 <- loo::loo_compare (model_objects_list$m0.loo_fit,
                                 model_objects_list$m1.loo_fit)

  dif_m0_m2 <- loo::loo_compare (model_objects_list$m0.loo_fit,
                                 model_objects_list$m2.loo_fit)

  model_objects_list$m0.model_summary$difference <- 0
  model_objects_list$m0.model_summary$diff_se <- 1

  model_objects_list$m1.model_summary$difference <- dif_m0_m1[1]
  model_objects_list$m1.model_summary$diff_se <- dif_m0_m1[2]

  model_objects_list$m2.model_summary$difference <- dif_m0_m2[1]
  model_objects_list$m2.model_summary$diff_se <- dif_m0_m2[2]

  comparison_table <- rbind(model_objects_list$m0.model_summary,
                            model_objects_list$m1.model_summary,
                            model_objects_list$m2.model_summary)


  return(comparison_table)
}


