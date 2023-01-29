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


#' Get Model Comparison
#'
#' Function that provides a statisctica comparison between models
#' @param model_objects_list model_objects to compare
#' @return compasiron table
#' @export
get_comparison_table <- function(model_objects_list) {
  mod1 <-model_objects_list[1]
  dif_m0_m1 <- loo_compare (mod_0$loo_fit, mod_1$loo_fit)


  return(comparison_table)
}
