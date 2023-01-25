#' Get Table Rhats
#'
#' Función que hace la tabla de los rhats
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
