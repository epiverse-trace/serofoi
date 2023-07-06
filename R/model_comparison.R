#' Method for extracting a dataframe containing the R-hat estimates for a given
#' serological model
#'
#' This method relies in the function [rhat][bayesplot::rhat] to extract the
#' R-hat estimates of the serological model object `seromodel_object` and
#' returns a table a dataframe with the estimates for each year of birth.
#' @inheritParams get_foi_central_estimates
#' @return rhats table
#' @examples
#' \dontrun{
#' data(chagas2012)
#' data_test <- prepare_serodata(serodata = chagas2012)
#' model_constant <- run_seromodel(serodata = data_test,
#'                                 foi_model = "constant",
#'                                 n_iters = 1500)
#' get_table_rhats(model_object = model_constant)
#' }
#' @export
get_table_rhats <- function(seromodel_object,
                            cohort_ages) {
  rhats <- bayesplot::rhat(seromodel_object, "foi")

  if (any(is.nan(rhats))) {
    rhats[which(is.nan(rhats))] <- NA
  }
  model_rhats <- data.frame(year = seromodel_object$exposure_years, rhat = rhats)
  return(model_rhats)
}
