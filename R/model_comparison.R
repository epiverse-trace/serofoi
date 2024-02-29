#' Build dataframe containing the R-hat estimates for a given
#' serological model
#'
#' This function relies on [rhat][bayesplot::rhat] to extract the
#' R-hat estimates of the serological model object `seromodel_object` and
#' returns a table a dataframe with the estimates for each year of birth.
#' @inheritParams get_foi_central_estimates
#' @return rhats table
#' @examples
#' data(chagas2012)
#' serodata <- prepare_serodata(serodata = chagas2012)
#' model_constant <- run_seromodel(
#'   serodata = serodata,
#'   foi_model = "constant",
#'   iter = 1500
#' )
#' cohort_ages <- get_cohort_ages(serodata)
#' get_table_rhats(
#'   seromodel_object = model_constant,
#'   cohort_ages = cohort_ages
#' )
#' @export
get_table_rhats <- function(seromodel_object,
                            cohort_ages) {
  rhats <- bayesplot::rhat(seromodel_object, "foi")

  if (any(is.nan(rhats))) {
    rhats[which(is.nan(rhats))] <- 0
  }
  model_rhats <- data.frame(year = cohort_ages$birth_year, rhat = rhats)
  model_rhats$rhat[model_rhats$rhat == 0] <- NA

  return(model_rhats)
}
