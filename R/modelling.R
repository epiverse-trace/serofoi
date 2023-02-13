#' Run model
#' Runs the specified stan model for the force-of-infection
#' @param model_data A data frame containing the data from a seroprevalence survey.
#' This data frame must contain the following columns:
#' \tabular{ll}{
#' \code{survey} \tab survey Label of the current survey \cr \tab \cr
#' \code{total} \tab Number of samples for each age group\cr \tab \cr
#' \code{counts} \tab Number of positive samples for each age group\cr \tab \cr
#' \code{age_min} \tab age_min \cr \tab \cr
#' \code{age_max} \tab age_max \cr \tab \cr
#' \code{year_init} \tab year_init \cr \tab \cr
#' \code{year_end} \tab year_end \cr \tab \cr
#' \code{tsur} \tab Year in which the survey took place \cr \tab \cr
#' \code{country} \tab The country where the survey took place \cr \tab \cr
#' \code{test} \tab The type of test taken \cr \tab \cr
#' \code{antibody} \tab antibody \cr \tab \cr
#' \code{age_mean_f} \tab Floor value of the average between age_min and age_max \cr \tab \cr
#' \code{sample_size} \tab The size of the sample \cr \tab \cr
#' \code{birth_year} \tab The year in which the individuals of each age group were bornt \cr \tab \cr
#' \code{prev_obs} \tab Observed prevalence \cr \tab \cr
#' \code{prev_obs_lower} \tab Lower limit of the confidence interval for the observed prevalence \cr \tab \cr
#' \code{prev_obs_upper} \tab Upper limit of the confidence interval for the observed prevalence \cr \tab \cr
#' }
#' The last six colums can be added to \code{model_data} by means of the function \code{\link{prepare_data}}.
#' @param model_name Name of the selected model. Current version provides three options:
#' \describe{
#' \item{\code{constant_foi_bi}}{Runs a constant model}
#' \item{\code{continuous_foi_normal_bi}}{Runs a normal model}
#' \item{\code{continuous_foi_normal_log}}{Runs a normal logarithmic model}
#' }
#' @param n_iters Number of interations for eah chain including the warmup. \code{iter} in \link[rstan]{sampling}.
#' @param n_thin Positive integer specifying the period for saving samples. \code{thin} in \link[rstan]{sampling}.
#' @param delta Real number between 0 and 1 that represents the target average acceptance probability.
#' Increasing the value of \code{delta} will result in a smaller step size and fewer divergences.
#' For further details refer to the \code{control} parameter in \link[rstan]{sampling} or \href{https://mc-stan.org/rstanarm/reference/adapt_delta.html}{here}.
#' @param m_treed Maximum tree depth for the binary tree used in the NUTS stan sampler. For further details refer to the \code{control} parameter in \link[rstan]{sampling}.
#' @param decades Number of decades covered by the survey data.
#' @return model_object (complementar cuando escriba la documentación de fit_model y fit_model_log)
#' @examples
#' model_data <- preprare_data(mydata)
#' run_model (model_data,
#'            model_name = "constant_foi_bi")
#' @export
run_model <- function(model_data,
                      model_name = "constant_foi_bi",
                      n_iters = 1000,
                      n_thin = 2,
                      delta = 0.90,
                      m_treed = 10,
                      decades = 0) {
  survey <- unique(model_data$survey)
  if (length(survey) > 1) warning("You have more than 1 surveys or survey codes")

  if (model_name == "constant_foi_bi") {
    model_0 <- save_or_load_model(model_name = model_name)
    model_object <- fit_model(model_data = model_data,
                              model_name = model_name,
                              n_iters = n_iters,
                              n_thin = n_thin,
                              delta = delta,
                              m_treed = m_treed,
                              decades = decades); print(paste0("serofoi model ",
                                                               model_name,
                                                               " finished running ------"))
  }
  if (model_name == "continuous_foi_normal_bi") {
    model_1 <- save_or_load_model(model_name = model_name)
    model_object <- fit_model(model_data = model_data,
                              model_name = model_name,
                              n_iters = n_iters,
                              n_thin = n_thin,
                              delta = delta,
                              m_treed = m_treed,
                              decades = decades); print(paste0("serofoi model ",
                                                               model_name,
                                                               " finished running ------"))
  }
  if (model_name == "continuous_foi_normal_log") {
    model_object <- fit_model_log(model_data = model_data,
                                  model_name = model_name,
                                  n_iters = n_iters,
                                  n_thin = n_thin,
                                  delta = delta,
                                  m_treed = m_treed,
                                  decades = decades); print(paste0("serofoi model ",
                                                                   model_name,
                                                                   " finished running ------"))
  }
  print(t(model_object$model_summary))
  return(model_object)
}

#' Save or load model
#' This function determines whether the corresponding .RDS file of the selected model exists or not.
#' In case the .RDS file exists, it is read and returned; otherwise, the object model is created through the \link[rstan]{stan_model} function, saved as an .RDS file and returned as the output of the function.
#' @param Name of the selected model. Current version provides three options:
#' \describe{
#' \item{\code{constant_foi_bi}}{Runs a constant model}
#' \item{\code{continuous_foi_normal_bi}}{Runs a normal model}
#' \item{\code{continuous_foi_normal_log}}{Runs a normal logarithmic model}
#' }
#' @return \code{model}. The rstan model object corresponding to the selected model.
#' @examples
#' save_or_load_model(model_name = "constant_foi_bi")
#' @export
save_or_load_model <- function(model_name = "constant_foi_bi") {
  rds_path <- config::get(model_name)$rds_path
  stan_path <- config::get(model_name)$stan_path

  if (!file.exists(rds_path)) {
    warning(paste0("Model ", model_name, " is being compiled for the first time. This might take some minutes"))
    model <- rstan::stan_model(stan_path)
    saveRDS(model, rds_path)
  } else {
    model <- readRDS(rds_path)
  }
  return(model)
}


#' Fit Model
#'
#' Function that fits the selected model to the data
#' @param model_data A data frame containing the data from a seroprevalence survey. For further details refer to \link{run_model}.
#' @param model_name Name of the selected model. Current version provides three options:
#' \describe{
#' \item{\code{constant_foi_bi}}{Runs a constant model}
#' \item{\code{continuous_foi_normal_bi}}{Runs a normal model}
#' \item{\code{continuous_foi_normal_log}}{Runs a normal logarithmic model}
#' }
#' @param n_iters Number of interations for eah chain including the warmup. \code{iter} in \link[rstan]{sampling}.
#' @param n_thin Positive integer specifying the period for saving samples. \code{thin} in \link[rstan]{sampling}.
#' @param delta Real number between 0 and 1 that represents the target average acceptance probability.
#' Increasing the value of \code{delta} will result in a smaller step size and fewer divergences.
#' For further details refer to the \code{control} parameter in \link[rstan]{sampling} or \href{https://mc-stan.org/rstanarm/reference/adapt_delta.html}{here}.
#' @param m_treed Maximum tree depth for the binary tree used in the NUTS stan sampler. For further details refer to the \code{control} parameter in \link[rstan]{sampling}.
#' @param decades Number of decades covered by the survey data.
#' @return model_object
#' @examples
#' model_data <- prepare_data(mydata)
#' fit_model (model_data,
#'            model_name = "constant_foi_bi",
#'            n_iters = n_iters,
#'            n_thin = n_thin,
#'            delta = delta,
#'            m_treed = m_treed,
#'            decades = decades)
#' @export
fit_model <- function(model_data,
                      model_name,
                      n_iters = 1000,
                      n_thin = 2,
                      delta = 0.90,
                      m_treed = 10,
                      decades = 0) {
  # add a warning because there are exceptions where a minimal amount of iterations need to be run
  model <- save_or_load_model(model_name)

  exposure_years <- get_exposure_years(model_data)
  exposure_years <- exposure_years[-length(exposure_years)]
  real_exposure_years <- (min(model_data$birth_year):model_data$tsur[1])[-1]
  exposure_matrix <- get_exposure_matrix(model_data, exposure_years)
  Nobs <- nrow(model_data)

  stan_data <- list(
    Nobs = Nobs,
    Npos = model_data$counts,
    Ntotal = model_data$total,
    Age = model_data$age_mean_f,
    Ymax = max(exposure_years),
    AgeExpoMatrix = exposure_matrix,
    NDecades = decades
  )

  n_warmup <- floor(n_iters / 2)

  f_init <- function() {
    list(foi = rep(0.01, length(exposure_years)))
  }

  fit <- rstan::sampling(
    model,
    data = stan_data,
    iter = n_iters,
    chains = 4,
    init = f_init,
    warmup = n_warmup,
    verbose = FALSE,
    refresh = 0,
    control = list(adapt_delta = delta,
                   max_treedepth = m_treed),
    seed = "12345",
    thin = n_thin
  )

  if (class(fit@sim$samples) != "NULL") {
    loo_fit <- loo::loo(fit, save_psis = TRUE, "logLikelihood")
    foi <- rstan::extract(fit, "foi", inc_warmup = FALSE)[[1]]

    # generates central estimations
    foi_cent_est <- data.frame(
      year = real_exposure_years,
      lower = apply(foi, 2, function(x) quantile(x, 0.05)),

      upper = apply(foi, 2, function(x) quantile(x, 0.95)),

      medianv = apply(foi, 2, function(x) quantile(x, 0.5))
    )


    # generates a sample of iterations
    if (n_iters >= 2000) {
      foi_post_s <- dplyr::sample_n(as.data.frame(foi), size = 1000)
      colnames(foi_post_s) <- real_exposure_years
    } else {
      foi_post_s <- as.data.frame(foi)
      colnames(foi_post_s) <- real_exposure_years
    }

    model_object <- list(
      fit = fit,
      model_data = model_data,
      stan_data = stan_data,
      real_exposure_years = real_exposure_years,
      exposure_years = exposure_years,
      n_iters = n_iters,
      n_thin = n_thin,
      n_warmup = n_warmup,
      model = model_name,
      delta = delta,
      m_treed = m_treed,
      loo_fit = loo_fit,
      foi_cent_est = foi_cent_est,
      foi_post_s = foi_post_s
    )
    model_object$model_summary <-
      extract_model_summary(model_object)
  } else {
    loo_fit <- c(-1e10, 0)
    model_object <- list(
      fit = "no model",
      model_data = model_data,
      stan_data = stan_data,
      real_exposure_years = real_exposure_years,
      exposure_years = exposure_years,
      n_iters = n_iters,
      n_thin = n_thin,
      n_warmup = n_warmup,
      model = model_name,
      delta = delta,
      m_treed = m_treed,
      loo_fit = loo_fit,
      model_summary = NA
    )
  }

  return(model_object)
}

#' Fit Model Log
#'
#' Function that fits the logarithmic model to the data
#' @param model_data model_data
#' @param model refers to model selected
#' @param model_name name of the model selected
#' @param n_iters number of iterations. This value comes by default but it can be changed
#' @param n_thinThis value comes by default but it can be changed
#' @param delta This value comes by default but it can be changed
#' @param m_treed This value comes by default but it can be changed
#' @param decades The decades covered by the survey data
#' @return model_object
#' @examples
#' fit_model_log (model,
#'                model_data,
#'                model_name)
#' @export
fit_model_log <- function(model_data,
                          model_name,
                          n_iters = 3000,
                          n_thin = 2,
                          delta = 0.90,
                          m_treed = 10,
                          decades = 0) {
  model <- save_or_load_model(model_name)
  exposure_years <- get_exposure_years(model_data)
  exposure_years <- exposure_years[-length(exposure_years)]
  real_exposure_years <- (min(model_data$birth_year):model_data$tsur[1])[-1]
  exposure_matrix <- get_exposure_matrix(model_data, exposure_years)
  Nobs <- nrow(model_data)

  stan_data <- list(
    Nobs = nrow(model_data),
    Npos = model_data$counts,
    Ntotal = model_data$total,
    Age = model_data$age_mean_f,
    Ymax = max(exposure_years),
    AgeExpoMatrix = exposure_matrix,
    NDecades = decades
  )

  n_warmup <- floor(n_iters / 2)

  f_init <- function() {
    list(log_foi = rep(-3, length(exposure_years)))
  }

  fit <- rstan::sampling(
    model,
    data = stan_data,
    iter = n_iters,
    chains = 4,
    init = f_init,
    warmup = n_warmup,
    verbose = FALSE,
    refresh = 0,
    control = list(adapt_delta = delta,
                   max_treedepth = m_treed),
    seed = "12345",
    thin = n_thin
  )

  if (class(fit@sim$samples) != "NULL") {
    loo_fit <- loo::loo(fit, save_psis = TRUE, "logLikelihood")
    foi <- rstan::extract(fit, "foi", inc_warmup = FALSE)[[1]]

    foi_cent_est <- data.frame(
      year = real_exposure_years,
      lower = apply(foi, 2, function(x) quantile(x, 0.1)),

      upper = apply(foi, 2, function(x) quantile(x, 0.9)),

      medianv = apply(foi, 2, function(x) quantile(x, 0.5))
    )

    foi_post_s <- dplyr::sample_n(as.data.frame(foi), size = 1000)
    colnames(foi_post_s) <- real_exposure_years


    model_object <- list(
      fit = fit,
      model_data = model_data,
      stan_data = stan_data,
      real_exposure_years = real_exposure_years,
      exposure_years = exposure_years,
      n_iters = n_iters,
      n_thin = n_thin,
      n_warmup = n_warmup,
      model = model_name,
      delta = delta,
      m_treed = m_treed,
      loo_fit = loo_fit,
      foi_cent_est = foi_cent_est,
      foi_post_s = foi_post_s
    )

    model_object$model_summary <-
      extract_model_summary(model_object)
  } else {
    loo_fit <- c(-1e10, 0)
    model_object <- list(
      fit = "no model",
      stan_data = stan_data,
      real_exposure_years = real_exposure_years,
      exposure_years = exposure_years,
      n_iters = n_iters,
      n_thin = n_thin,
      n_warmup = n_warmup,
      model = model_name,
      delta = delta,
      m_treed = m_treed,
      loo_fit = loo_fit
    )
    model_object$model_summary <-
      extract_model_summary(model_object)
  }

  return(model_object)
}


#' Get exposure years
#'
#' Function that generates an atomic vector with the exposition years in model_data. The exposition years to the disease for each individual corresponds to the time from birth to the moment of the survey.
#' @param model_data A data frame containing the data from a seroprevalence survey. This data frame must contain the year of birth for each individual (birth_year) and the time of the survey (tsur). birth_year can be constructed by means of the <prepare_data> function.
#' @return exposure_years. An atomic vector with the numeration of the exposition years in model_data
#' @examples
#' model_data <- prepare_data(mydata)
#' exposure_years <- get_exposure_years (model_data)
#' @export
get_exposure_years <- function(model_data) {
  exposure_years <- (seq_along(min(model_data$birth_year):model_data$tsur[1]))
}


#' Get Exposure Matrix
#'
#' Function that generates the exposure matrix for a seroprevalence survey.
#' @param model_data A data frame containing the data from a seroprevalence survey. This data frame must contain the year of birth for each individual (birth_year) and the time of the survey (tsur). birth_year can be constructed by means of the <prepare_data> function.
#' @return exposure_output
#' @examples
#' model_data <- prepare_data(mydata)
#' exposure_matrix <- get_exposure_matrix(model_data)
#' @export
get_exposure_matrix <- function(model_data,
                                exposure_years) {
  age_class <- model_data$age_mean_f
  ly <- length(exposure_years)
  exposure <- matrix(0, nrow = length(age_class), ncol = ly)
  for (k in 1:length(age_class))
    exposure[k, (ly - age_class[k] + 1):ly] <- 1
  exposure_output <- exposure
  return(exposure_output)
}


#' Extract Model Summary
#'
#' Function that summarizes the models
#' @param model_object what the run model function returns
#' @param model_data refers to data of the model
#' @return summary of the models
#' @examples
#' extract_model_summary (model_object)
#' @export
extract_model_summary <- function(model_object) {
  model_name <- model_object$model
  #------- Loo estimates

  loo_fit <- model_object$loo_fit
  if (sum(is.na(loo_fit)) < 1) {
    lll <- as.numeric((round(loo_fit$estimates[1, ], 2)))
  } else {
    lll <- c(-1e10, 0)
  }

  model_data <- model_object$model_data
  summary_model <- data.frame(
    model = model_object$model,
    dataset = model_data$survey[1],
    country = model_data$country[1],
    year = model_data$tsur[1],
    test = model_data$test[1],
    antibody = model_data$antibody[1],
    n_sample = sum(model_data$total),
    n_agec = length(model_data$age_mean_f),
    n_iter = model_object$n_iters,
    elpd = lll[1],
    se = lll[2],
    converged = NA
  )

  rhats <- get_table_rhats(model_object)
  if (any(rhats$rhat > 1.1) == FALSE) {
    summary_model$converged <- "Yes"
  }

  return(summary_model)
}


#' Get Prevalence Expanded
#'
#' Function that generates the expanded prevalence
#' @param model_data refers to the model data that has been selected
#' @param foi force of infection
#' @return prev_final
#' @examples
#' get_prev_expanded <- function(foi, model_data)
#' @export


get_prev_expanded <- function(foi,
                              model_data) {
  dim_foi <- dim(foi)[2]
  if (dim_foi < 80) {
    oldest_year <- 80 - dim_foi + 1
    foin <- matrix(NA, nrow = dim(foi)[1], 80)
    foin[, oldest_year:80] <- foi
    foin[, 1:(oldest_year - 1)] <- rowMeans(foi[, 1:5])
  } else {
    foin <- foi
  }

  foi_expanded <- foin

  age_class <- 1:NCOL(foi_expanded)
  ly <- NCOL(foi_expanded)
  exposure <- matrix(0, nrow = length(age_class), ncol = ly)
  for (k in 1:length(age_class))
    exposure[k, (ly - age_class[k] + 1):ly] <- 1
  exposure_expanded <- exposure

  iterf <- NROW(foi_expanded)
  age_max <- NROW(exposure_expanded)
  prev_pn <- matrix(NA, nrow = iterf, ncol = age_max)
  for (i in 1:iterf) {
    prev_pn[i, ] <- 1 - exp(-exposure_expanded %*% foi_expanded[i, ])
  }

  lower <- apply(prev_pn, 2, function(x) quantile(x, 0.1))

  upper <- apply(prev_pn, 2, function(x) quantile(x, 0.9))

  medianv <- apply(prev_pn, 2, function(x) quantile(x, 0.5))

  predicted_prev <- data.frame(
    age = 1:80,
    predicted_prev = medianv,
    predicted_prev_lower = lower,
    predicted_prev_upper = upper
  )

  observed_prev <- model_data %>%
    dplyr::select(age_mean_f,
                  prev_obs,
                  prev_obs_lower,
                  prev_obs_upper,
                  total,
                  counts) %>%
    dplyr::rename(age = age_mean_f,
                  sample_by_age = total,
                  positives = counts)

  prev_expanded <-
    base::merge(predicted_prev,
                observed_prev,
                by = "age",
                all.x = TRUE) %>% dplyr::mutate(survey = model_data$survey[1])

  # I added this here for those cases when binned is prefered for plotting
  if (model_data$age_max[1] - model_data$age_min[1] < 3) {
  xx <- prepare_bin_data(model_data)
    prev_final <-
      base::merge(prev_expanded, xx, by = "age", all.x = TRUE)
  } else {
    prev_final <- prev_expanded %>% dplyr::mutate(
      cut_ages = "original",
      bin_size = .data$sample_by_age,
      bin_pos = .data$positives,
      p_obs_bin = .data$prev_obs,
      p_obs_bin_l = .data$prev_obs_lower,
      p_obs_bin_u = .data$prev_obs_upper
    )
  }

  return(prev_final)
}


#' Get Posterior Summary
#'
#' Function that gets a posterior summary
#' @param model_objects model_objects_chain
#' @return model_object
#' @examples
#' get_posterior_summary (model_objects_chain)
#' @export
get_posterior_summary <- function(model_objects_chain) {
  model_object <- sapply(model_objects_chain,
                         function(i) c(quantile(i, c(0.5, 0.025, 0.975))))
  row.names(model_object) <- c("Median", "Lower", "Upper")
  return(model_object)
}
