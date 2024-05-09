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
#' model_constant <- fit_seromodel(
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
  checkmate::assert_class(seromodel_object, "stanfit")

  rhats <- bayesplot::rhat(seromodel_object, "foi")

  if (any(is.nan(rhats))) {
    warn_msg <- paste0(
      length(which(is.nan(rhats))),
      " rhat values are `nan`, ",
      "indicating the model may not have run correctly for those times.\n",
      "Setting those rhat values to `NA`."
    )
    warning(warn_msg)
    rhats[which(is.nan(rhats))] <- NA
  }
  model_rhats <- data.frame(
    year = cohort_ages$birth_year,
    rhat = rhats
    )

  return(model_rhats)
}
