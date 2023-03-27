#' Get Table Rhats
#'
#' Function that makes the rhats table
#' @param model_object model_object
#' @return rhats table
#' @examples
#' \dontrun{
#' serodata <- prepare_serodata(serodata = serodata, alpha = 0.05)
#' model_object <- run_seromodel(
#' serodata = serodata, seromodel_name = "constant_foi_bi")
#' get_table_rhats (model_object)
#' }
#' @export
get_table_rhats <- function(model_object) {
  rhats <- bayesplot::rhat(model_object$fit, "foi")

  if (any(is.nan(rhats))) {
    rhats[which(is.nan(rhats))] <- 0
  }
  model_rhats <- data.frame(year = model_object$exposure_years, rhat = rhats)
  model_rhats$rhat[model_rhats$rhat == 0] <- NA

  return(model_rhats)
}