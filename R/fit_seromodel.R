#' Adds age group marker to serosurvey
#'
#' @inheritParams fit_seromodel
add_age_group_to_serosurvey <- function(serosurvey) {
  if (any(colnames(serosurvey) == "age_group")) {
    message("Using `age_group` already present in serosurvey")
  } else {
    serosurvey <- serosurvey %>%
      dplyr::mutate(
        age_group = floor((.data$age_min + .data$age_max) / 2)
      )
    }
  return(serosurvey)
}

#' Sets initialization function for sampling
#'
#' @inheritParams fit_seromodel
#' @param foi_init Initialization function for sampling. If null, default is chosen
#' depending on the foi-scale of the model
#' @export
set_foi_init <- function(
  foi_init,
  is_log_foi,
  foi_index
) {

  # Set default behavior for initialization
  if (is.null(foi_init)) {
    config_file <- system.file("extdata", "config.yml", package = "serofoi")
    init_default <- config::get(file = config_file, "priors")$defaults$init

    if (is_log_foi) {
      foi_init <- function() {
        list(log_foi_vector = rep(log(init_default), max(foi_index)))
      }
    } else {
      foi_init <- function() {
        list(foi_vector = rep(init_default, max(foi_index)))
      }
    }
  }

  checkmate::assert_class(foi_init, "function")
  checkmate::assert_double(unlist(foi_init()[[1]]), len = max(foi_index))

  return(foi_init)
}

#' Runs specified stan model for the force-of-infection
#'
#' @param serosurvey
#' \describe{
#'   \item{`survey_year`}{Year in which the survey took place
#'         (only needed to plot time models)}
#'   \item{`age_min`}{Floor value of the average between age_min and age_max}
#'   \item{`age_max`}{The size of the sample}
#'   \item{`n_sample`}{Number of samples for each age group}
#'   \item{`n_seropositive`}{Number of positive samples for each age group}
#' }
#' @param model_type Type of the model. Either "constant", "age" or "time"
#' @param is_log_foi Boolean to set logarithmic scale in the FOI
#' @param foi_prior Force-of-infection distribution specified by means of
#'  the helper functions. Currently available options are:
#' \describe{
#'  \item{[sf_normal]}{Function to set normal distribution priors}
#'  \item{[sf_uniform]}{Function to set uniform distribution priors}
#' }
#' @param foi_sigma_rw Prior distribution for the standard deviation of the
#' force-of-infection. Currently available options are:
#' \describe{
#'  \item{[sf_normal]}{Function to set normal distribution prior.
#'                     Available for time models in the log-scale}
#'  \item{[sf_cauchy]}{Function to set Cauchy distribution prior.
#'                     Available for time models in regular scale.}
#' }
#' @param foi_index Integer vector specifying the age-groups for which
#' force-of-infection values will be estimated. It can be specified by
#' means of [get_foi_index]
#' @inheritParams set_foi_init
#' @param is_seroreversion Boolean specifying whether to include
#' seroreversion rate estimation in the model
#' @param seroreversion_prior seroreversion distribution specified by means of
#'  the helper functions. Currently available options are:
#' \describe{
#'  \item{[sf_normal]}{Function to set normal distribution priors}
#'  \item{[sf_uniform]}{Function to set uniform distribution priors}
#'  \item{[sf_none]}{Function to set no prior distribution}
#' }
#' @param ... Additional parameters for [rstan][rstan::sampling]
#' @returns stan_fit object with force-of-infection and seroreversion
#' (when applicable) samples
#' @examples
#' data(veev2012)
#' seromodel <- fit_seromodel(
#' serosurvey = veev2012,
#'   model_type = "time",
#'   foi_index = get_foi_index(veev2012, group_size = 30, model_type = "time")
#' )
#' @export
fit_seromodel <- function(
  serosurvey,
  model_type = "constant",
  is_log_foi = FALSE,
  foi_prior = sf_normal(),
  foi_sigma_rw = sf_none(),
  foi_index = NULL,
  foi_init = NULL,
  is_seroreversion = FALSE,
  seroreversion_prior = sf_normal(),
  ...
) {
  serosurvey <- serosurvey %>%
    validate_serosurvey() %>%
    add_age_group_to_serosurvey()

  stopifnot(
    "model_type must be either 'constant', 'time' or 'age'" =
    model_type %in% c("constant", "time", "age")
  )

  stan_data <- build_stan_data(
    serosurvey = serosurvey,
    model_type = model_type,
    foi_prior = foi_prior,
    foi_index = foi_index,
    is_log_foi = is_log_foi,
    foi_sigma_rw = foi_sigma_rw,
    is_seroreversion = is_seroreversion,
    seroreversion_prior = seroreversion_prior
  )

  foi_init <- set_foi_init(
    foi_init = foi_init,
    is_log_foi = is_log_foi,
    foi_index = stan_data$foi_index
  )

  # Assigning right name to the model based on user specifications
  model_name <- model_type
  if (is_log_foi) {
    model_name <- paste0(model_name, "_log")
  }
  if (is_seroreversion)
    model_name <- paste0(model_name, "_seroreversion")
  else
    model_name <- paste0(model_name, "_no_seroreversion")

  # Compile or load Stan model
  model <- stanmodels[[model_name]]

  seromodel <- rstan::sampling(
    model,
    data = stan_data,
    init = foi_init,
    ...
  )
  seromodel@model_name <- model_name
  return(seromodel)
}
