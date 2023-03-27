#' Run the specified stan model for the Force-of-Infection and estimates de seroprevalence based on the result of the fit.
#'
#' @param serodata A data frame containing the data from a seroprevalence survey.
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
#' The last six colums can be added to \code{serodata} by means of the function \code{\link{prepare_serodata}}.
#' @param seromodel_name Name of the selected model. Current version provides three options:
#' \describe{
#' \item{\code{"constant_foi_bi"}}{Runs a constant model}
#' \item{\code{"continuous_foi_normal_bi"}}{Runs a normal model}
#' \item{\code{"continuous_foi_normal_log"}}{Runs a normal logarithmic model}
#' }
#' @param n_iters Number of interations for eah chain including the warmup. \code{iter} in \link[rstan]{sampling}.
#' @param n_thin Positive integer specifying the period for saving samples. \code{thin} in \link[rstan]{sampling}.
#' @param delta Real number between 0 and 1 that represents the target average acceptance probability.
#' Increasing the value of \code{delta} will result in a smaller step size and fewer divergences.
#' For further details refer to the \code{control} parameter in \link[rstan]{sampling} or \href{https://mc-stan.org/rstanarm/reference/adapt_delta.html}{here}.
#' @param m_treed Maximum tree depth for the binary tree used in the NUTS stan sampler. For further details refer to the \code{control} parameter in \link[rstan]{sampling}.
#' @param decades Number of decades covered by the survey data.
#' @return \code{seromodel_object}. An object containing relevant information about the implementation of the model. For further details refer to \link{fit_seromodel}.
#' @examples
#' \dontrun{
#' serodata <- prepare_serodata(serodata)
#' run_seromodel (serodata,
#'            seromodel_name = "constant_foi_bi")
#' }
#' @export
run_seromodel <- function(serodata,
                      seromodel_name = "constant_foi_bi",
                      n_iters = 1000,
                      n_thin = 2,
                      delta = 0.90,
                      m_treed = 10,
                      decades = 0,
                      print_summary = TRUE) {
  survey <- unique(serodata$survey)
  if (length(survey) > 1) warning("You have more than 1 surveys or survey codes")
  seromodel_object <- fit_seromodel(serodata = serodata,
                                    seromodel_name = seromodel_name,
                                    n_iters = n_iters,
                                    n_thin = n_thin,
                                    delta = delta,
                                    m_treed = m_treed,
                                    decades = decades); print(paste0("serofoi model ",
                                                                      seromodel_name,
                                                                      " finished running ------"))
  if (print_summary){
    print(t(seromodel_object$model_summary))
  }
  return(seromodel_object)
}

#' Save or load model
#'
#' This function determines whether the corresponding .RDS file of the selected model exists or not.
#' In case the .RDS file exists, it is read and returned; otherwise, the object model is created through the
#' \link[rstan]{stan_model} function, saved as an .RDS file and returned as the output of the function.
#' @param seromodel_name Name of the selected model. Current version provides three options:
#' \describe{
#' \item{\code{"constant_foi_bi"}}{Runs a constant model}
#' \item{\code{"continuous_foi_normal_bi"}}{Runs a normal model}
#' \item{\code{"continuous_foi_normal_log"}}{Runs a normal logarithmic model}
#' }
#' @return \code{model}. The rstan model object corresponding to the selected model.
#' @examples
#' save_or_load_model(seromodel_name="constant_foi_bi")
#' @export

save_or_load_model <- function(seromodel_name = "constant_foi_bi") {
  base_path <- config::get("stan_models_base_path",
    file = system.file("config.yml", package = "serofoi", mustWork = TRUE))
  rds_path <- system.file(base_path, paste(seromodel_name, ".rds", sep = ""), package = getPackageName())
  if (!file.exists(rds_path)) {
    message(sprintf("\nNo rds file found for model %s. Compiling stan model...", seromodel_name))
  }
  stan_path <- system.file(base_path, paste(seromodel_name, ".stan", sep = ""), package = getPackageName())

  model <- rstan::stan_model(stan_path, auto_write = TRUE)

  return(model)
}


#' Fit Model
#'
#' Function that fits the selected model to the data
#' @param serodata A data frame containing the data from a seroprevalence survey. For further details refer to \link{run_seromodel}.
#' @param seromodel_name Name of the selected model. Current version provides three options:
#' \describe{
#' \item{\code{"constant_foi_bi"}}{Runs a constant model}
#' \item{\code{"continuous_foi_normal_bi"}}{Runs a normal model}
#' \item{\code{"continuous_foi_normal_log"}}{Runs a normal logarithmic model}
#' }
#' @param n_iters Number of interations for eah chain including the warmup. \code{iter} in \link[rstan]{sampling}.
#' @param n_thin Positive integer specifying the period for saving samples. \code{thin} in \link[rstan]{sampling}.
#' @param delta Real number between 0 and 1 that represents the target average acceptance probability.
#' Increasing the value of \code{delta} will result in a smaller step size and fewer divergences.
#' For further details refer to the \code{control} parameter in \link[rstan]{sampling} or \href{https://mc-stan.org/rstanarm/reference/adapt_delta.html}{here}.
#' @param m_treed Maximum tree depth for the binary tree used in the NUTS stan sampler. For further details refer to the \code{control} parameter in \link[rstan]{sampling}.
#' @param decades Number of decades covered by the survey data.
#' @return \code{seromodel_object}. An object containing relevant information about the implementation of the model. It contains the following:
#' \tabular{ll}{
#' \code{fit} \tab \code{stanfit} object returned by the function \link[rstan]{sampling} \cr \tab \cr
#' \code{serodata} \tab A data frame containing the data from a seroprevalence survey. For further details refer to \link{run_seromodel}.\cr \tab \cr
#' \code{stan_data} \tab List containing \code{Nobs}, \code{Npos}, \code{Ntotal}, \code{Age}, \code{Ymax}, \code{AgeExpoMatrix} and \code{NDecades}.
#' This object is used as an input for the \link[rstan]{sampling} function \cr \tab \cr
#' \code{exposure_years} \tab Integer atomic vector containing the actual exposure years (1946, ..., 2007 e.g.) \cr \tab \cr
#' \code{exposure_ages} \tab Integer atomic vector containing the numeration of the exposure ages. \cr \tab \cr
#' \code{n_iters} \tab Number of interations for eah chain including the warmup. \cr \tab \cr
#' \code{n_thin} \tab Positive integer specifying the period for saving samples. \cr \tab \cr
#' \code{n_warmup} \tab Number of warm up iterations. Set by default as n_iters/2. \cr \tab \cr
#' \code{seromodel_name} \tab The name of the model\cr \tab \cr
#' \code{delta} \tab Real number between 0 and 1 that represents the target average acceptance probability. \cr \tab \cr
#' \code{m_treed} \tab Maximum tree depth for the binary tree used in the NUTS stan sampler. \cr \tab \cr
#' \code{loo_fit} \tab Efficient approximate leave-one-out cross-validation. Refer to \link[loo]{loo} for further details. \cr \tab \cr
#' \code{foi_cent_est} \tab A data fram e containing \code{year} (corresponding to \code{exposure_years}), \code{lower}, \code{upper}, and \code{medianv} \cr \tab \cr
#' \code{foi_post_s} \tab Sample n rows from a table. Refer to \link[dplyr]{sample_n} for further details. \cr \tab \cr
#' \code{model_summary} \tab A data fram containing the summary of the model. Refer to \link{extract_seromodel_summary} for further details. \cr \tab \cr
#' }

#' @examples
#' \dontrun{
#' serodata <- prepare_serodata(serodata)
#' fit_seromodel (serodata,
#'            seromodel_name = "constant_foi_bi")
#' }
#'
#' @export
fit_seromodel <- function(serodata,
                      seromodel_name,
                      n_iters = 1000,
                      n_thin = 2,
                      delta = 0.90,
                      m_treed = 10,
                      decades = 0) {
  # add a warning because there are exceptions where a minimal amount of iterations need to be run
  model <- save_or_load_model(seromodel_name)
  exposure_ages <- get_exposure_ages(serodata)
  exposure_years <- (min(serodata$birth_year):serodata$tsur[1])[-1]
  exposure_matrix <- get_exposure_matrix(serodata)
  Nobs <- nrow(serodata)

  stan_data <- list(
    Nobs = Nobs,
    Npos = serodata$counts,
    Ntotal = serodata$total,
    Age = serodata$age_mean_f,
    Ymax = max(exposure_ages),
    AgeExpoMatrix = exposure_matrix,
    NDecades = decades
  )

  n_warmup <- floor(n_iters / 2)

  if (seromodel_name == "continuous_foi_normal_log") {
    f_init <- function() {
      list(log_foi = rep(-3, length(exposure_ages)))
    }
    lower_quantile = 0.1
    upper_quantile = 0.9
    medianv_quantile = 0.5
  }

  else {
  f_init <- function() {
    list(foi = rep(0.01, length(exposure_ages)))
  }
    lower_quantile = 0.05
    upper_quantile = 0.95
    medianv_quantile = 0.5

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
    thin = n_thin,
    chain_id = 0 # https://github.com/stan-dev/rstan/issues/761#issuecomment-647029649
  )

  if (class(fit@sim$samples) != "NULL") {
    loo_fit <- loo::loo(fit, save_psis = TRUE, "logLikelihood")
    foi <- rstan::extract(fit, "foi", inc_warmup = FALSE)[[1]]
    # foi <- rstan::extract(fit, "foi", inc_warmup = TRUE, permuted=FALSE)[[1]]
    # generates central estimations
    foi_cent_est <- data.frame(
      year = exposure_years,
      lower = apply(foi, 2, function(x) quantile(x, lower_quantile)),

      upper = apply(foi, 2, function(x) quantile(x, upper_quantile)),

      medianv = apply(foi, 2, function(x) quantile(x, medianv_quantile))
    )


    # generates a sample of iterations
    if (n_iters >= 2000) {
      foi_post_s <- dplyr::sample_n(as.data.frame(foi), size = 1000)
      colnames(foi_post_s) <- exposure_years
    } else {
      foi_post_s <- as.data.frame(foi)
      colnames(foi_post_s) <- exposure_years
    }

    seromodel_object <- list(
      fit = fit,
      serodata = serodata,
      stan_data = stan_data,
      exposure_years = exposure_years,
      exposure_ages = exposure_ages,
      n_iters = n_iters,
      n_thin = n_thin,
      n_warmup = n_warmup,
      seromodel_name = seromodel_name,
      delta = delta,
      m_treed = m_treed,
      loo_fit = loo_fit,
      foi_cent_est = foi_cent_est,
      foi_post_s = foi_post_s
    )
    seromodel_object$model_summary <-
      extract_seromodel_summary(seromodel_object)
  } else {
    loo_fit <- c(-1e10, 0)
    seromodel_object <- list(
      fit = "no model",
      serodata = serodata,
      stan_data = stan_data,
      exposure_years = exposure_years,
      exposure_ages = exposure_ages,
      n_iters = n_iters,
      n_thin = n_thin,
      n_warmup = n_warmup,
      model = seromodel_name,
      delta = delta,
      m_treed = m_treed,
      loo_fit = loo_fit,
      model_summary = NA
    )
  }

  return(seromodel_object)
}


#' Get exposure years
#'
#' Function that generates an atomic vector with the exposition years in serodata. The exposition years to the disease for each individual corresponds to the time from birth to the moment of the survey.
#' @param serodata A data frame containing the data from a seroprevalence survey. This data frame must contain the year of birth for each individual (birth_year) and the time of the survey (tsur). birth_year can be constructed by means of the \link{prepare_serodata} function.
#' @return \code{exposure_ages}. An atomic vector with the numeration of the exposition years in serodata
#' @examples
#' \dontrun{
#' serodata <- prepare_serodata(serodata = serodata, alpha = 0.05)
#' exposure_ages <- get_exposure_ages(serodata)
#' }
#' @export
get_exposure_ages <- function(serodata) {
  # TODO Verify if this change is correct
  return(seq_along(min(serodata$birth_year):(serodata$tsur[1] - 1)))
}


#' Get Exposure Matrix
#'
#' Function that generates the exposure matrix for a seroprevalence survey.
#' @param serodata A data frame containing the data from a seroprevalence survey. This data frame must contain the year of birth for each individual (birth_year) and the time of the survey (tsur). birth_year can be constructed by means of the \link{prepare_serodata} function.
#' @return \code{exposure_output}. An atomic matrix containing the expositions for each entry of \code{serodata} by year.
#' @examples
#' \dontrun{
#' serodata <- prepare_serodata(serodata = serodata, alpha = 0.05)
#' exposure_matrix <- get_exposure_matrix(serodata = serodata)
#' }
#' @export
get_exposure_matrix <- function(serodata) {
  age_class <- serodata$age_mean_f
  exposure_ages <- get_exposure_ages(serodata)
  ly <- length(exposure_ages)
  exposure <- matrix(0, nrow = length(age_class), ncol = ly)
  for (k in 1:length(age_class))
    exposure[k, (ly - age_class[k] + 1):ly] <- 1
  exposure_output <- exposure
  return(exposure_output)
}


#' Extract Model Summary
#'
#' Function to generate a summary of a model.
#' @param seromodel_object \code{seromodel_object}. An object containing relevant information about the implementation of the model. Refer to \link{fit_seromodel} for further details.
#' @return \code{model_summary}. Object with a summary of \code{seromodel_object} containing the following:
#' \tabular{ll}{
#' \code{seromodel_name} \tab Name of the selected model. For further details refer to \link{save_or_load_model}. \cr \tab \cr
#' \code{data_set} \tab Seroprevalence survey label.\cr \tab \cr
#' \code{country} \tab Name of the country were the survey was conducted in. \cr \tab \cr
#' \code{year} \tab Year in which the survey was conducted. \cr \tab \cr
#' \code{test} \tab Type of test of the survey. \cr \tab \cr
#' \code{antibody} \tab Antibody \cr \tab \cr
#' \code{n_sample} \tab Total number of samples in the survey. \cr \tab \cr
#' \code{n_agec} \tab Number of age groups considered. \cr \tab \cr
#' \code{n_iter} \tab Number of interations for eah chain including the warmup. \cr \tab \cr
#' \code{elpd} \tab elpd \cr \tab \cr
#' \code{se} \tab se \cr \tab \cr
#' \code{converged} \tab convergence \cr \tab \cr
#' }
#' @examples
#' \dontrun{
#' serodata <- prepare_serodata(serodata)
#' seromodel_object <- run_seromodel(serodata = serodata,
#'                           seromodel_name = "constant_foi_bi")
#' extract_seromodel_summary (seromodel_object)
#' }
#' @export
extract_seromodel_summary <- function(seromodel_object) {
  seromodel_name <- seromodel_object$seromodel_name
  serodata <- seromodel_object$serodata
  #------- Loo estimates

  loo_fit <- seromodel_object$loo_fit
  if (sum(is.na(loo_fit)) < 1) {
    lll <- as.numeric((round(loo_fit$estimates[1, ], 2)))
  } else {
    lll <- c(-1e10, 0)
  }
  model_summary <- data.frame(
    seromodel_name = seromodel_name,
    dataset = serodata$survey[1],
    country = serodata$country[1],
    year = serodata$tsur[1],
    test = serodata$test[1],
    antibody = serodata$antibody[1],
    n_sample = sum(serodata$total),
    n_agec = length(serodata$age_mean_f),
    n_iter = seromodel_object$n_iters,
    elpd = lll[1],
    se = lll[2],
    converged = NA
  )

  rhats <- get_table_rhats(seromodel_object)
  if (any(rhats$rhat > 1.1) == FALSE) {
    model_summary$converged <- "Yes"
  }

  return(model_summary)
}


#' Get Prevalence Expanded
#'
#' Function that generates the expanded prevalence
#' @param serodata A data frame containing the data from a seroprevalence survey. For further details refer to \link{run_seromodel}.
#' @param foi Object containing the information of the force of infection. It is obtained from \code{rstan::extract(seromodel_object$fit, "foi", inc_warmup = FALSE)[[1]]}.
#' @return \code{prev_final}. The expanded prevalence data. This is used for plotting purposes in the \code{visualization} module.
#' @examples
#' \dontrun{
#' serodata <- prepare_serodata(serodata)
#' seromodel_object <- run_seromodel(serodata = serodata,
#'                           seromodel_name = "constant_foi_bi")
#' foi <- rstan::extract(seromodel_object$fit, "foi", inc_warmup = FALSE)[[1]]
#' get_prev_expanded <- function(foi, serodata)
#' }
#' @export
get_prev_expanded <- function(foi,
                              serodata) {
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

  observed_prev <- serodata %>%
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
                all.x = TRUE) %>% dplyr::mutate(survey = serodata$survey[1])

  # I added this here for those cases when binned is prefered for plotting
  if (serodata$age_max[1] - serodata$age_min[1] < 3) {
  xx <- prepare_bin_data(serodata)
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
