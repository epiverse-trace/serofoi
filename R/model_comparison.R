#' Method for extracting a dataframe containing the R-hat estimates for a given serological model
#' 
#' This method relies in the function \link[bayesplot]{rhat} to extract the R-hat estimates of the serological model object 
#' \code{seromodel_object} and returns a table a dataframe with the estimates for each year of birth. 
#' @param seromodel_object seromodel_object
#' @return rhats table
#' @examples
#' data(chagas2012)
#' serodata <- prepare_serodata(serodata = chagas2012)
#' model_constant <- run_seromodel(serodata = serodata,
#'                                 foi_model = "constant",
#'                                 n_iters = 1500)
#' get_table_rhats(seromodel_object = model_constant,
#'                 serodata = serodata)
#' @export
get_table_rhats <- function(seromodel_object, serodata) {
  rhats <- bayesplot::rhat(seromodel_object, "foi")

  if (any(is.nan(rhats))) {
    rhats[which(is.nan(rhats))] <- 0
  }
  cohort_ages <- get_cohort_ages(serodata)
  model_rhats <- data.frame(year = cohort_ages$birth_year, rhat = rhats)
  model_rhats$rhat[model_rhats$rhat == 0] <- NA

  return(model_rhats)
}