sf_normal <- function(mean = 0, sd = 1) {
  # Restricting normal inputs to be non-negative
  if(mean < 0 | sd <= 0) {
    msg <- paste0(
      "Normal distribution here only accepts",
      " non-negative values for mean and standard deviation"
    )
    message(msg)
    stop()
  }

  return(list(mean = mean, sd = sd, name = "normal"))
}

sf_uniform <- function(min = 0, max = 10) {
  # Restricting uniform inputs to be non-negative
  if (min < 0 | max < 0) {
    msg <- paste0(
      "Uniform distribution here only accepts",
      " non-negative values for min and max"
    )
    message(msg)
    stop()
  }
  if (min >= max) {
    message("Uniform distribution only accepts min < max")
  }

  return(list(min = min, max = max, name = "uniform"))
}

sf_none <- function() {
  return(list(name = "none"))
}

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

    foi_index <- append(
      foi_index,
      rep(
        max(foi_index),
        max(serosurvey$age_max) - length(foi_index)
      )
    )

  return(foi_index)
}

set_stan_data_defaults <- function(
    stan_data,
    is_seroreversion = FALSE
) {
  config_file <- "inst/extdata/config.yml"
  prior_default <- config::get(file = config_file, "priors")$defaults

  foi_defaults <- list(
    foi_prior_index = prior_default$index,
    foi_min = prior_default$min,
    foi_max = prior_default$max,
    foi_mean = prior_default$mean,
    foi_sd = prior_default$sd
  )
  stan_data <- append(
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
    stan_data <- append(
      stan_data,
      seroreversion_defaults
    )
  }

  return(stan_data)
}

build_stan_data <- function(
    serosurvey,
    model_type = "constant",
    foi_prior = sf_uniform(),
    foi_index = NULL,
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
    set_stan_data_defaults(is_seroreversion = is_seroreversion)

  if (is.null(foi_index)) {
    foi_index_default <- get_foi_index(serosurvey = serosurvey, group_size = 1)
    stan_data <- append(
      stan_data,
      list(foi_index = foi_index_default)
    )
  } else {
    # TODO: check that foi_index is the right size
    stan_data <- append(
      stan_data,
      list(foi_index = foi_index)
    )
  }

  config_file <- "inst/extdata/config.yml"
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

  if (is_seroreversion) {
    if(seroreversion_prior$name == "none") {
      message("seroreversion_prior not specified")
      stop()
    }
    else if (seroreversion_prior$name == "uniform") {
      stan_data$seroreversion_prior_index <- prior_index[["uniform"]]
      stan_data$seroreversion_min <- seroreversion_prior$min
      stan_data$seroreversion_max <- seroreversion_prior$max
    }
    else if (seroreversion_prior$name == "normal") {
      stan_data$seroreversion_prior_index <- prior_index[["normal"]]
      stan_data$seroreversion_mean <- seroreversion_prior$mean
      stan_data$seroreversion_sd <- seroreversion_prior$sd
    }
  }

  return(stan_data)
}
