#' Get Table Rhats
#'
#' Function that makes the rhats table
#' @param model_object model_object
#' @return rhats table
#' @examples
#' \dontrun{
#' seroprev_data <- prepare_seroprev_data(seroprev_data = serodata, alpha = 0.05)
#' model_object <- run_seroprev_model(
#' seroprev_data = seroprev_data, seroprev_model_name = "constant_foi_bi")
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