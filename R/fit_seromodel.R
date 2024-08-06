#' Adds age group marker to serosurvey
#'
#' @inheritParams fit_seromodel
add_age_group_to_serosurvey <- function(serosurvey) {
  if (!any(colnames(serosurvey) == "age_group")) {
    serosurvey <- serosurvey %>%
      dplyr::mutate(
        age_group = floor((.data$age_min + .data$age_max) / 2)
      )
  } else {
    message("Using `age_group` already present in serosurvey")
  }
  return(serosurvey)
}

#' Runs specified stan model for the force-of-infection
#'
#' @param serosurvey
#' \describe{
#'   \item{`tsur`}{Year in which the survey took place}
#'   \item{`age_min`}{Floor value of the average between age_min and age_max}
#'   \item{`age_max`}{The size of the sample}
#'   \item{`sample_size`}{Number of samples for each age group}
#'   \item{`n_seropositive`}{Number of positive samples for each age group}
#' }
#' @param model_type Type of the model. Either "constant", "age" or "time"
#' @param foi_prior Force-of-infection distribution specified by means of
#'  the helper functions. Currently available options are:
#' \describe{
#'  \item{`[sf_normal]`}
#'  \item{`[sf_uniform]`}
#'  \item{`[sf_none]`}
#' }
#' @param foi_index Integer vector specifying the age-groups for which
#' force-of-infection values will be estimated
#' @param is_seroreversion Boolean specifying whether to include
#' seroreversion rate estimation in the model
#' @param seroreversion_prior seroreversion distribution specified by means of
#'  the helper functions. Currently available options are:
#' \describe{
#'  \item{`[sf_normal]`}
#'  \item{`[sf_uniform]`}
#'  \item{`[sf_none]`}
#' }
#' @param ... Additional parameters for [rstan][rstan::sampling]
#' @returns stan_fit object with force-of-infection and seroreversion
#' (when applicable) samples
#' @export
fit_seromodel <- function(
  serosurvey,
  model_type = "constant",
  foi_prior = sf_normal(),
  foi_index = NULL,
  is_seroreversion = FALSE,
  seroreversion_prior = sf_uniform(),
  ...
) {
  serosurvey <- serosurvey %>%
    validate_serosurvey() %>%
    add_age_group_to_serosurvey()

  stan_data <- build_stan_data(
    serosurvey = serosurvey,
    model_type = model_type,
    foi_prior = foi_prior,
    foi_index = foi_index,
    is_seroreversion = is_seroreversion,
    seroreversion_prior = seroreversion_prior
  )

  if (is_seroreversion)
    model_name <- paste0(model_type, "_seroreversion")
  else
    model_name <- paste0(model_type, "_no_seroreversion")

  # model <- stan_models[[model_name]]
  model <- stanmodels[[model_name]]
  seromodel <- rstan::sampling(
    model,
    data = stan_data,
    ...
  )
  seromodel@model_name <- model_name
  return(seromodel)
}
