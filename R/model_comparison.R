#' Method for extracting a dataframe containing the R-hat estimates for a given serological model
#' 
#' This method relies in the function \link[bayesplot]{rhat} to extract the R-hat estimates of the serological model object 
#' \code{seromodel_object} and returns a table a dataframe with the estimates for each year of birth. 
#' @param seromodel_object seromodel_object
#' @return rhats table
#' @examples
#' \dontrun{
#' data("serodata")
#' data_test <- prepare_serodata(serodata = serodata)
#' model_constant <- run_seromodel(serodata = data_test, 
#'                                 foi_model = "constant", 
#'                                 n_iters = 1500)
#' get_table_rhats(model_constant, data_test)
#' }
#' @export
get_table_rhats <- function(seromodel_object,
                            serodata) {
  rhats <- bayesplot::rhat(seromodel_object, "foi")
  exposure_years <- (min(serodata$birth_year):serodata$tsur[1])[-1]
  if (any(is.nan(rhats))) {
    rhats[which(is.nan(rhats))] <- 0
  }
  model_rhats <- data.frame(year = exposure_years, rhat = rhats)
  model_rhats$rhat[model_rhats$rhat == 0] <- NA

  return(model_rhats)
}