# TODO Complete @param documentation


#' Function that runs the specified stan model for the Force-of-Infection and
#' estimates the seroprevalence based on the result of the fit
#'
#' This function runs the specified model for the Force-of-Infection `foi_model`
#' using the data from a seroprevalence survey `serodata` as the input data. See
#' [fit_seromodel] for further details.
#'
#' @param serodata A data frame containing the data from a seroprevalence
#'   survey. This data frame must contain the following columns:
#' \describe{
#'   \item{`survey`}{survey Label of the current survey}
#'   \item{`total`}{Number of samples for each age group}
#'   \item{`counts`}{Number of positive samples for each age group}
#'   \item{`age_min`}{age_min}
#'   \item{`age_max`}{age_max}
#'   \item{`tsur`}{Year in which the survey took place}
#'   \item{`country`}{The country where the survey took place}
#'   \item{`test`}{The type of test taken}
#'   \item{`antibody`}{antibody}
#'   \item{`age_mean_f`}{Floor value of the average between age_min and age_max}
#'   \item{`sample_size`}{The size of the sample}
#'   \item{`birth_year`}{The year in which the individuals of each age group
#'     were born}
#'   \item{`prev_obs`}{Observed prevalence}
#'   \item{`prev_obs_lower`}{Lower limit of the confidence interval for the
#'     observed prevalence}
#'   \item{`prev_obs_upper`}{Upper limit of the confidence interval for the
#'     observed prevalence}
#' }
#'   The last six columns can be added to `serodata` by means of the function
#'   [prepare_serodata()].
#' @param foi_model Name of the selected model. Current version provides three
#'   options:
#' \describe{
#' \item{`"constant"`}{Runs a constant model}
#' \item{`"tv_normal"`}{Runs a normal model}
#' \item{`"tv_normal_log"`}{Runs a normal logarithmic model}
#' }
#' @param n_iters Number of interactions for each chain including the warmup.
#'   `iter` in [sampling][rstan::sampling].
#' @param n_thin Positive integer specifying the period for saving samples.
#'   `thin` in [sampling][rstan::sampling].
#' @param delta Real number between 0 and 1 that represents the target average
#'   acceptance probability. Increasing the value of `delta` will result in a
#'   smaller step size and fewer divergences. For further details refer to the
#'   `control` parameter in [sampling][rstan::sampling] or
#'   [here](https://mc-stan.org/rstanarm/reference/adapt_delta.html).
#' @param m_treed Maximum tree depth for the binary tree used in the NUTS stan
#'   sampler. For further details refer to the `control` parameter in
#'   [sampling][rstan::sampling].
#' @param decades Number of decades covered by the survey data.
#' @param print_summary TBD
#' @return `seromodel_object`. An object containing relevant information about
#'   the implementation of the model. For further details refer to
#'   [fit_seromodel].
#' @examples
#' data(chagas2012)
#' serodata <- prepare_serodata(chagas2012)
#' run_seromodel(serodata,
#'   foi_model = "constant"
#' )
#' @export
run_seromodel <- function(serodata,
                          foi_model = c(
                            "constant", "tv_normal_log",
                            "tv_normal"
                          ),
                          n_iters = 1000,
                          n_thin = 2,
                          delta = 0.90,
                          m_treed = 10,
                          decades = 0,
                          print_summary = TRUE) {
  foi_model <- match.arg(foi_model)
  survey <- unique(serodata$survey)
  if (length(survey) > 1) {
    warning("You have more than 1 surveys or survey codes")
  }
  seromodel_object <- fit_seromodel(
    serodata = serodata,
    foi_model = foi_model,
    n_iters = n_iters,
    n_thin = n_thin,
    delta = delta,
    m_treed = m_treed,
    decades = decades
  )
  message(
    "serofoi model ",
    foi_model,
    " finished running ------"
  )
  if (print_summary) {
    model_summary <- extract_seromodel_summary(
      seromodel_object = seromodel_object,
      serodata = serodata
    )
    print(t(model_summary))
  }
  return(seromodel_object)
}

#' Function that fits the selected model to the specified seroprevalence survey
#' data
#'
#' This function fits the specified model `foi_model` to the serological survey
#' data `serodata` by means of the [sampling][rstan::sampling] method. The
#' function determines whether the corresponding stan model object needs to be
#' compiled by rstan.
#' @inheritParams run_seromodel
#' @param foi_model Name of the selected model. Current version provides three
#'   options:
#' \describe{
#' \item{`"constant"`}{Runs a constant model}
#' \item{`"tv_normal"`}{Runs a normal model}
#' \item{`"tv_normal_log"`}{Runs a normal logarithmic model}
#' }
#' @param n_iters Number of interactions for each chain including the warmup.
#'   `iter` in [sampling][rstan::sampling].
#' @param n_thin Positive integer specifying the period for saving samples.
#'   `thin` in [sampling][rstan::sampling].
#' @param delta Real number between 0 and 1 that represents the target average
#'   acceptance probability. Increasing the value of `delta` will result in a
#'   smaller step size and fewer divergences. For further details refer to the
#'   `control` parameter in [sampling][rstan::sampling] or
#'   [here](https://mc-stan.org/rstanarm/reference/adapt_delta.html).
#' @param m_treed Maximum tree depth for the binary tree used in the NUTS stan
#'   sampler. For further details refer to the `control` parameter in
#'   [sampling][rstan::sampling].
#' @param decades Number of decades covered by the survey data.
#' @return `seromodel_object`. `stanfit` object returned by the function
#'   [sampling][rstan::sampling]
#' @examples
#' data(chagas2012)
#' serodata <- prepare_serodata(chagas2012)
#' seromodel_fit <- fit_seromodel(
#'   serodata = serodata,
#'   foi_model = "constant"
#' )
#'
#' @export
fit_seromodel <- function(serodata,
                          foi_model = c(
                            "constant", "tv_normal_log",
                            "tv_normal"
                          ),
                          n_iters = 1000,
                          n_thin = 2,
                          delta = 0.90,
                          m_treed = 10,
                          decades = 0) {
  # TODO Add a warning because there are exceptions where a minimal amount of
  # iterations is needed
  foi_model <- match.arg(foi_model)
  model <- stanmodels[[foi_model]]
  cohort_ages <- get_cohort_ages(serodata = serodata)
  exposure_matrix <- get_exposure_matrix(serodata)
  Nobs <- nrow(serodata)

  stan_data <- list(
    Nobs = Nobs,
    Npos = serodata$counts,
    Ntotal = serodata$total,
    Age = serodata$age_mean_f,
    Ymax = max(cohort_ages$age),
    AgeExpoMatrix = exposure_matrix,
    NDecades = decades
  )

  n_warmup <- floor(n_iters / 2)
  if (foi_model == "tv_normal_log") {
    f_init <- function() {
      list(log_foi = rep(-3, nrow(cohort_ages)))
    }
  } else {
    f_init <- function() {
      list(foi = rep(0.01, nrow(cohort_ages)))
    }
  }

  seromodel_fit <- rstan::sampling(
    model,
    data = stan_data,
    iter = n_iters,
    chains = 4,
    init = f_init,
    warmup = n_warmup,
    verbose = FALSE,
    refresh = 0,
    control = list(
      adapt_delta = delta,
      max_treedepth = m_treed
    ),
    seed = "12345",
    thin = n_thin,
    # https://github.com/stan-dev/rstan/issues/761#issuecomment-647029649
    chain_id = 0
  )

  if (seromodel_fit@mode == 0) {
    seromodel_object <- seromodel_fit
    return(seromodel_object)
  } else {
    # This may happen for invalid inputs in rstan::sampling() (e.g. thin > iter)
    seromodel_object <- "no model"
    return(seromodel_object)
  }
}


#' Function that generates a data.frame containing the age of each cohort
#' corresponding to each birth year excluding the year of the survey.
#'
#' This function generates a data.frame containing the age of each cohort
#' corresponding to each `birth_year` excluding the year of the survey, for
#' which the cohort age is still 0. specified serological survey data `serodata`
#' excluding the year of the survey.
#' @inheritParams run_seromodel
#' @return `cohort_ages`. A data.frame containing the age of each cohort
#'   corresponding to each birth year
#' @examples
#' data(chagas2012)
#' serodata <- prepare_serodata(serodata = chagas2012, alpha = 0.05)
#' cohort_ages <- get_cohort_ages(serodata = serodata)
#' @export
get_cohort_ages <- function(serodata) {
  birth_year <- (min(serodata$birth_year):serodata$tsur[1])
  age <- (seq_along(min(serodata$birth_year):(serodata$tsur[1] - 1)))

  cohort_ages <- data.frame(
    birth_year = birth_year[-length(birth_year)],
    age = rev(age)
  )
  return(cohort_ages)
}

# TODO Is necessary to explain better what we mean by the exposure matrix.

#' Function that generates the exposure matrix corresponding to a serological
#' survey
#'
#' @inheritParams run_seromodel
#' @return `exposure_output`. An atomic matrix containing the expositions for
#'   each entry of `serodata` by year.
#' @examples
#' data(chagas2012)
#' serodata <- prepare_serodata(serodata = chagas2012)
#' exposure_matrix <- get_exposure_matrix(serodata = serodata)
#' @keywords internal
#' @noRd
get_exposure_matrix <- function(serodata) {
  age_class <- serodata$age_mean_f
  cohort_ages <- get_cohort_ages(serodata = serodata)
  ly <- nrow(cohort_ages)
  exposure <- matrix(0, nrow = length(age_class), ncol = ly)
  for (k in seq_along(age_class)) {
    exposure[k, (ly - age_class[k] + 1):ly] <- 1
  }
  exposure_output <- exposure
  return(exposure_output)
}

#' Function that generates the central estimates for the fitted forced FoI
#'
#' @param seromodel_object Stanfit object containing the results of fitting a
#'   model by means of [run_seromodel].
#' @param cohort_ages  A data.frame containing the age of each cohort
#'   corresponding to each birth year.
#' @return `foi_central_estimates`. Central estimates for the fitted forced FoI
#' @examples
#' data(chagas2012)
#' serodata <- prepare_serodata(chagas2012)
#' seromodel_object <- fit_seromodel(
#'   serodata = serodata,
#'   foi_model = "constant"
#' )
#' cohort_ages <- get_cohort_ages(serodata = serodata)
#' foi_central_estimates <- get_foi_central_estimates(
#'   seromodel_object = seromodel_object,
#'   cohort_ages = cohort_ages
#' )
#' @export
get_foi_central_estimates <- function(seromodel_object,
                                      cohort_ages) {
  if (seromodel_object@model_name == "tv_normal_log") {
    lower_quantile <- 0.1
    upper_quantile <- 0.9
    medianv_quantile <- 0.5
  } else {
    lower_quantile <- 0.05
    upper_quantile <- 0.95
    medianv_quantile <- 0.5
  }
  # extracts foi from stan fit
  foi <- rstan::extract(seromodel_object, "foi", inc_warmup = FALSE)[[1]]

  # generates central estimations
  foi_central_estimates <- data.frame(
    year = cohort_ages$birth_year,
    lower = apply(foi, 2, quantile, lower_quantile),
    upper = apply(foi, 2, quantile, upper_quantile),
    medianv = apply(foi, 2, quantile, medianv_quantile)
  )
  return(foi_central_estimates)
}

#' Function to extract a summary of the specified serological model object
#'
#' This function extracts a summary corresponding to a serological model object
#' containing information about the original serological survey data used to
#' fit the model, such as the year when the survey took place, the type of test
#' taken and the corresponding antibody, as well as information about the
#' convergence of the model, like the expected log pointwise predictive density
#' `elpd` and its corresponding standard deviation.
#' @inheritParams get_foi_central_estimates
#' @inheritParams run_seromodel
#' @return `model_summary`. Object with a summary of `seromodel_object`
#'   containing the following:
#' \describe{
#'   \item{`foi_model`}{Name of the selected model.}
#'   \item{`data_set`}{Seroprevalence survey label}
#'   \item{`country`}{Name of the country were the survey was conducted in.}
#'   \item{`year`}{Year in which the survey was conducted.}
#'   \item{`test`}{Type of test of the survey.}
#'   \item{`antibody`}{Antibody}
#'   \item{`n_sample`}{Total number of samples in the survey.}
#'   \item{`n_agec`}{Number of age groups considered.}
#'   \item{`n_iter`}{Number of interactions for each chain including the warmup.}
#'   \item{`elpd`}{elpd}
#'   \item{`se`}{se}
#'   \item{`converged`}{convergence}
#' }
#' @examples
#' data(chagas2012)
#' serodata <- prepare_serodata(chagas2012)
#' seromodel_object <- run_seromodel(
#'   serodata = serodata,
#'   foi_model = "constant"
#' )
#' extract_seromodel_summary(seromodel_object,
#'   serodata = serodata
#' )
#' @export
extract_seromodel_summary <- function(seromodel_object,
                                      serodata) {
  #------- Loo estimates
  # The argument parameter_name refers to the name given to the Log-likelihood
  # in the stan models. See loo::extract_log_lik() documentation for further
  # details
  loo_fit <- loo::loo(
    seromodel_object,
    save_psis = FALSE,
    pars = c(parameter_name = "logLikelihood")
  )
  if (sum(is.na(loo_fit)) < 1) {
    lll <- as.numeric((round(loo_fit$estimates[1, ], 2)))
  } else {
    lll <- c(-1e10, 0)
  }
  #-------
  model_summary <- data.frame(
    foi_model = seromodel_object@model_name,
    dataset = unique(serodata$survey),
    country = unique(serodata$country),
    year = unique(serodata$tsur),
    test = unique(serodata$test),
    antibody = unique(serodata$antibody),
    n_sample = sum(serodata$total),
    n_agec = length(serodata$age_mean_f),
    n_iter = seromodel_object@sim$iter,
    elpd = lll[1],
    se = lll[2],
    converged = NA
  )
  cohort_ages <- get_cohort_ages(serodata = serodata)
  rhats <- get_table_rhats(
    seromodel_object = seromodel_object,
    cohort_ages = cohort_ages
  )
  if (!any(rhats$rhat > 1.1)) {
    model_summary$converged <- "Yes"
  }

  return(model_summary)
}

# TODO Complete @param documentation

#' Function that generates an object containing the confidence interval based on
#' a Force-of-Infection fitting
#'
#' This function computes the corresponding binomial confidence intervals for
#' the obtained prevalence based on a fitting of the Force-of-Infection `foi`
#' for plotting an analysis purposes.
#' @param foi Object containing the information of the force of infection. It is
#'   obtained from `rstan::extract(seromodel_object$seromodel, "foi", inc_warmup
#'   = FALSE)[[1]]`.
#' @inheritParams run_seromodel
#' @param bin_data TBD
#' @return `prev_final`. The expanded prevalence data. This is used for plotting
#'   purposes in the `visualization` module.
#' @examples
#' data(chagas2012)
#' serodata <- prepare_serodata(chagas2012)
#' seromodel_object <- run_seromodel(
#'   serodata = serodata,
#'   foi_model = "constant"
#' )
#' foi <- rstan::extract(seromodel_object, "foi")[[1]]
#' get_prev_expanded(foi, serodata)
#' @keywords internal
#' @noRd
get_prev_expanded <- function(foi,
                              serodata,
                              bin_data = FALSE) {
  dim_foi <- dim(foi)[2]
  # TODO: check whether this conditional is necessary
  if (dim_foi < 80) {
    oldest_year <- 80 - dim_foi + 1
    foin <- matrix(NA, nrow = dim(foi)[1], 80)
    foin[, oldest_year:80] <- foi
    foin[, 1:(oldest_year - 1)] <- rowMeans(foi[, 1:5])
  } else {
    foin <- foi
  }

  foi_expanded <- foin

  ly <- NCOL(foi_expanded)
  exposure_expanded <- matrix(0, nrow = ly, ncol = ly)
  exposure_expanded[lower.tri(exposure_expanded, diag = TRUE)] <- 1

  iterf <- NROW(foi_expanded)
  age_max <- NROW(exposure_expanded)
  prev_pn <- matrix(NA, nrow = iterf, ncol = age_max)
  for (i in 1:iterf) {
    prev_pn[i, ] <- 1 - exp(-exposure_expanded %*% foi_expanded[i, ])
  }

  lower <- apply(prev_pn, 2, quantile, 0.1)

  upper <- apply(prev_pn, 2, quantile, 0.9)

  medianv <- apply(prev_pn, 2, quantile, 0.5)

  predicted_prev <- data.frame(
    age = 1:age_max,
    predicted_prev = medianv,
    predicted_prev_lower = lower,
    predicted_prev_upper = upper
  )

  observed_prev <- serodata %>%
    dplyr::select(
      "age_mean_f",
      "prev_obs",
      "prev_obs_lower",
      "prev_obs_upper",
      "total",
      "counts"
    ) %>%
    dplyr::rename(
      age = "age_mean_f",
      sample_by_age = "total",
      positives = "counts"
    )

  prev_expanded <-
    base::merge(predicted_prev,
      observed_prev,
      by = "age",
      all.x = TRUE
    ) %>% dplyr::mutate(survey = serodata$survey[1])
  if (bin_data) {
    # I added this here for those cases when binned is prefered for plotting
    if (serodata$age_max[1] - serodata$age_min[1] < 3) {
      xx <- prepare_bin_data(serodata)
      prev_expanded <-
        base::merge(prev_expanded, xx, by = "age", all.x = TRUE)
    } else {
      prev_expanded <- dplyr::mutate(
        prev_expanded,
        cut_ages = "original",
        bin_size = .data$sample_by_age,
        bin_pos = .data$positives,
        p_obs_bin = .data$prev_obs,
        p_obs_bin_l = .data$prev_obs_lower,
        p_obs_bin_u = .data$prev_obs_upper
      )
    }
  }

  return(prev_expanded)
}
