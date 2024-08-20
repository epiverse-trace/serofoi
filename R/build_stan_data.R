#' Sets normal distribution parameters for sampling
#'
#' @param mean Mean of the normal distribution
#' @param sd Standard deviation of the normal distribution
#' @return List with specified statistics and name of the model
#' @export
sf_normal <- function(mean = 0, sd = 1) {
  # Restricting normal inputs to be non-negative
  if (mean < 0 || sd <= 0) {
    stop(
      "Normal distribution only accepts",
      " `mean>=0` and `sd>0` for mean and standard deviation")
  }

  return(list(mean = mean, sd = sd, name = "normal"))
}

#' Sets uniform distribution parameters for sampling
#'
#' @param min Minimum value of the random variable of the uniform distribution
#' @param max Maximum value of the random variable of the uniform distribution
#' @return List with specified statistics and name of the model
#' @export
sf_uniform <- function(min = 0, max = 10) {
  # Restricting uniform inputs to be non-negative
  if (min < 0 || (min >= max)) {
    stop(
      "Uniform distribution only accepts",
      " 0<=min<max"
    )
  }
  if (min >= max) {
    message("Uniform distribution only accepts min < max")
  }

  return(list(min = min, max = max, name = "uniform"))
}

#' Sets Cauchy distribution parameters for sampling
#'
#' @param scale Scale
#'  of the normal distribution
#' @param sd Standard deviation of the normal distribution
#' @return List with specified statistics and name of the model
#' @export
sf_cauchy <- function(location = 0, scale = 1) {
  # Restricting normal inputs to be non-negative
  if (location < 0 || scale < 0) {
    stop(
      "Normal distribution only accepts",
      " `location>=0` and `scale=>0` for mean and standard deviation")
  }

  return(list(location = location, scale = scale, name = "cauchy"))
}

#' Sets empty distribution
sf_none <- function() {
  return(list(name = "none"))
}

#' Generates force-of-infection indexes for heterogeneous age groups
#'
#' The max value of the force-of-infection indexes correspond to
#' the number of foi values to be estimated when sampling.
#' @inheritParams fit_seromodel
#' @param group_size Age groups size
#' @return Integer vector with the indexes numerating each year/age
#' (depending on the model).
#' @export
get_foi_index <- function(
  serosurvey,
  group_size
  ) {
    checkmate::assert_int(
      group_size,
      lower = 1,
      upper = max(serosurvey$age_max)
      )

    foi_index <- unlist(
      purrr::map(
        seq(
          1,
          max(serosurvey$age_max) / group_size,
          1),
        rep,
        times = group_size
      )
    )

    foi_index <- c(
      foi_index,
      rep(
        max(foi_index),
        max(serosurvey$age_max) - length(foi_index)
      )
    )

  return(foi_index)
}

#' Set stan data defaults for sampling
#'
#' @param stan_data List to be passed to [rstan][rstan::sampling]
#' @inheritParams fit_seromodel
#' @return List with default values of stan data for sampling
set_stan_data_defaults <- function(
    stan_data,
    is_log_foi = FALSE,
    is_seroreversion = FALSE
) {
  config_file <- system.file("extdata", "config.yml", package = "serofoi")
  prior_default <- config::get(file = config_file, "priors")$defaults

  foi_defaults <- list(
    foi_prior_index = prior_default$index,
    foi_min = prior_default$min,
    foi_max = prior_default$max,
    foi_mean = prior_default$mean,
    foi_sd = prior_default$sd
  )
  # Add sigma defaults depending on scale
  if (is_log_foi) {
    # Normal distribution defaults
    foi_defaults <- c(
      foi_defaults,
      list(
        foi_sigma_rw_loc = prior_default$mean,
        foi_sigma_rw_sc = prior_default$sd
      )
    )
  } else {
    # Cauchy distribution defaults
    foi_defaults <- c(
      foi_defaults,
      list(
        foi_sigma_rw_loc = prior_default$location,
        foi_sigma_rw_sc = prior_default$scale
      )
    )
  }
  stan_data <- c(
    stan_data,
    foi_defaults
  )

  if (is_seroreversion) {
    seroreversion_defaults <- list(
      seroreversion_prior_index = prior_default$index,
      seroreversion_min = prior_default$min,
      seroreversion_max = prior_default$max,
      seroreversion_mean = prior_default$mean,
      seroreversion_sd = prior_default$sd
    )
    stan_data <- c(
      stan_data,
      seroreversion_defaults
    )
  }

  return(stan_data)
}

#' Builds stan data for sampling depending on the selected model
#'
#' @inheritParams fit_seromodel
#' @return List with necessary data for sampling the specified model
#' @export
build_stan_data <- function(
    serosurvey,
    model_type = "constant",
    foi_prior = sf_uniform(),
    foi_index = NULL,
    is_log_foi = FALSE,
    foi_sigma_rw = sf_none(),
    is_seroreversion = FALSE,
    seroreversion_prior = sf_none()
) {

  stan_data <- list(
    n_observations = nrow(serosurvey),
    age_max = max(serosurvey$age_max),
    ages = seq(1, max(serosurvey$age_max), 1),
    n_seropositive = serosurvey$n_seropositive,
    sample_size = serosurvey$sample_size,
    age_groups = serosurvey$age_group
  ) %>%
    set_stan_data_defaults(
      is_log_foi = is_log_foi,
      is_seroreversion = is_seroreversion
      )

  if (is.null(foi_index)) {
    foi_index_default <- get_foi_index(serosurvey = serosurvey, group_size = 1)
    stan_data <- c(
      stan_data,
      list(foi_index = foi_index_default)
    )
  } else {
    # TODO: check that foi_index is the right size
    stan_data <- c(
      stan_data,
      list(foi_index = foi_index)
    )
  }
  config_file <- system.file("extdata", "config.yml", package = "serofoi")
  prior_index <- config::get(file = config_file, "priors")$indexes

  if (foi_prior$name == "uniform") {
    stan_data$foi_prior_index <- prior_index[["uniform"]]
    stan_data$foi_min <- foi_prior$min
    stan_data$foi_max <- foi_prior$max
  } else if (foi_prior$name == "normal") {
    stan_data$foi_prior_index <- prior_index[["normal"]]
    stan_data$foi_mean <- foi_prior$mean
    stan_data$foi_sd <- foi_prior$sd
  }

  if (foi_sigma_rw$name == "cauchy") {
    stan_data$foi_sigma_rw_loc <- foi_sigma_rw$location
    stan_data$foi_sigma_rw_sc <- foi_sigma_rw$scale
  }

  if (is_seroreversion) {
    if (seroreversion_prior$name == "none") {
      stop("seroreversion_prior not specified")
    } else if (seroreversion_prior$name == "uniform") {
      stan_data$seroreversion_prior_index <- prior_index[["uniform"]]
      stan_data$seroreversion_min <- seroreversion_prior$min
      stan_data$seroreversion_max <- seroreversion_prior$max
    } else if (seroreversion_prior$name == "normal") {
      stan_data$seroreversion_prior_index <- prior_index[["normal"]]
      stan_data$seroreversion_mean <- seroreversion_prior$mean
      stan_data$seroreversion_sd <- seroreversion_prior$sd
    }
  }

  return(stan_data)
}
