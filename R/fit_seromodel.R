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
  model <- rstan::stan_model(paste0("inst/stan/", model_name, ".stan"))

  seromodel <- rstan::sampling(
    model,
    data = stan_data,
    ...
  )

  return(seromodel)
}
