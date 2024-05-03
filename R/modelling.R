stop_if_missing <- function(serodata, must_have_cols) {
  if (
    !all(
      must_have_cols
      %in% colnames(serodata)
    )
  ) {
    missing <- must_have_cols[which(!(must_have_cols %in% colnames(serodata)))]
    stop(
      "The following mandatory columns in `serodata` are missing.\n",
      toString(missing)
    )
  }
}

stop_if_wrong_type <- function(serodata, col_types) {
  error_messages <- list()
  for (col in names(col_types)) {
    # valid_col_types <- ifelse(is.list(col_types[[col]]),
    #   col_types[[col]], as.list(col_types[[col]])
    # )
    valid_col_types <- as.list(col_types[[col]])

    # Only validates column type if the column exists in the dataframe
    if (col %in% colnames(serodata) &&
      !any(vapply(valid_col_types, function(type) {
        do.call(sprintf("is.%s", type), list(serodata[[col]]))
      }, logical(1)))) {
      error_messages <- append(
        error_messages,
        sprintf(
          "`%s` must be of any of these types: `%s`",
          col, toString(col_types[[col]])
        )
      )
    }
  }
  if (length(error_messages) > 0) {
    stop(
      "The following columns in `serodata` have wrong types:\n",
      toString(error_messages)
    )
  }
}

warn_missing <- function(
    serodata,
    optional_cols
    ) {
  if (
    !all(
      optional_cols
      %in% colnames(serodata)
    )
  ) {
    missing <- optional_cols[which(!(optional_cols %in% colnames(serodata)))]
    warning(
      "The following optional columns in `serodata` are missing. ",
      "Consider including them to get more information from this analysis:\n",
      toString(missing)
    )
    for (col in missing) {
      serodata[[col]] <- "None" # TODO Shouln't we use `NA` instead?
    }
  }

  return(serodata)
}


validate_serodata <- function(serodata) {
  col_types <- list(
    survey = c("character", "factor"),
    total = "numeric",
    counts = "numeric",
    tsur = "numeric",
    age_min = "numeric",
    age_max = "numeric"
  )

  stop_if_missing(serodata,
    must_have_cols = names(col_types)
  )

  stop_if_wrong_type(serodata, col_types)

  optional_col_types <- list(
    country = c("character", "factor"),
    test = c("character", "factor"),
    antibody = c("character", "factor")
  )

  # Add missing columns
  serodata <- warn_missing(
    serodata,
    optional_cols = names(optional_col_types)
  )

  # If any optional column is present, validates that is has the correct type
  stop_if_wrong_type(serodata, optional_col_types)

  return(serodata)
}

validate_prepared_serodata <- function(serodata) {
  col_types <- list(
    total = "numeric",
    counts = "numeric",
    tsur = "numeric",
    age_mean_f = "numeric",
    birth_year = "numeric",
    prev_obs = "numeric",
    prev_obs_lower = "numeric",
    prev_obs_upper = "numeric"
  )
  serodata <- validate_serodata(serodata)

  stop_if_missing(serodata, must_have_cols = names(col_types))
  stop_if_wrong_type(serodata, col_types)

  return(serodata)
}

#' Run specified stan model for the force-of-infection and
#' estimate the seroprevalence based on the result of the fit
#'
#' Starting on v.0.1.0, this function will be DEPRECATED. Use `fit_seromodel`
#' instead.
#' This function runs the specified model for the force-of-infection `foi_model`
#' using the data from a seroprevalence survey `serodata` as the input data. See
#' [fit_seromodel] for further details.
#'
#' @inheritParams fit_seromodel
#' @param print_summary Boolean. If `TRUE`, a table summarizing modelling
#' results is printed.
#' @return `seromodel_object`. An object containing relevant information about
#'   the implementation of the model. For further details refer to
#'   [fit_seromodel].
#' @examples
#' \dontrun{
#' data(chagas2012)
#' serodata <- prepare_serodata(chagas2012)
#' run_seromodel(
#'   serodata,
#'   foi_model = "constant"
#' )
#' }
#' @export
run_seromodel <- function(
    serodata,
    foi_model = c("constant", "tv_normal_log", "tv_normal"),
    foi_location = 0,
    foi_scale = 1,
    chunk_size = 1,
    chunks = NULL,
    iter = 1000,
    adapt_delta = 0.90,
    max_treedepth = 10,
    seed = 12345,
    print_summary = TRUE,
    ...) {
  .Deprecated("fit_seromodel")
  foi_model <- match.arg(foi_model)
  survey <- unique(serodata$survey)
  if (length(survey) > 1) {
    warning("You have more than 1 surveys or survey codes")
  }
  seromodel_object <- fit_seromodel(
    serodata = serodata,
    foi_model = foi_model,
    foi_location = foi_location,
    foi_scale = foi_scale,
    chunk_size = chunk_size,
    chunks = chunks,
    iter = iter,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    seed = seed,
    ...
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

#' Fit selected model to the specified seroprevalence survey
#' data
#'
#' This function fits the specified model `foi_model` to the serological survey
#' data `serodata` by means of [sampling][rstan::sampling]. The
#' function determines whether the corresponding stan model object needs to be
#' compiled by rstan.
#' @param serodata A data frame containing the data from a seroprevalence
#'   survey. This data frame must contain at least the following columns:
#' \describe{
#'   \item{`total`}{Number of samples for each age group}
#'   \item{`counts`}{Number of positive samples for each age group}
#'   \item{`tsur`}{Year in which the survey took place}
#'   \item{`age_mean_f`}{Floor value of the average between age_min and age_max}
#'   \item{`sample_size`}{The size of the sample}
#'   \item{`birth_year`}{The year in which the individuals of each age group
#'     were born}
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
#' @param foi_location Location parameter of the force-of-infection distribution
#' of the selected model. Depending on `foi_model`, the meaning may vary.
#' @param foi_scale Scale parameter of the force-of-infection distribution
#' of the selected model. Depending on `foi_model`, the meaning may vary.
#' @param chunks Numeric list specifying the chunk structure of the time
#' interval from the birth year of the oldest age cohort
#' `min(serodata$age_mean_f)` to the time when the serosurvey was conducted
#' `t_sur`. If `NULL`, the time interval is divided in chunks of size
#' `chunk_size`.
#' @param chunk_size Size of the chunks to be used in case that the chunk
#' structure `chunks` is not specified in [fit_seromodel].
#' Default is 1, meaning that one force of infection value is to be estimated
#' for every year in the time interval spanned by the serosurvey.
#' If the length of the time interval is not exactly divisible by `chunk_size`,
#' the remainder years are included in the last chunk.
#' @param iter Number of interactions for each chain including the warmup.
#'   `iter` in [sampling][rstan::sampling].
#' @param adapt_delta Real number between 0 and 1 that represents the target
#' average acceptance probability. Increasing the value of `adapt_delta` will
#' result in a smaller step size and fewer divergences. For further details
#' refer to the `control` parameter in [sampling][rstan::sampling] or
#' [here](https://mc-stan.org/rstanarm/reference/adapt_delta.html).
#' @param max_treedepth Maximum tree depth for the binary tree used in the NUTS
#' stan sampler. For further details refer to the `control` parameter in
#' [sampling][rstan::sampling].
#' @param seed For further details refer to the `seed` parameter in
#'   [sampling][rstan::sampling].
#' @param init_value optional list, each element of which is a named list 
#' containing initial values of the parameters for that MCMC chain.
#' @param ... Additional parameters for [sampling][rstan::sampling].
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
fit_seromodel <- function(
    serodata,
    foi_model = c("constant", "tv_normal_log", "tv_normal"),
    foi_location = 0,
    foi_scale = 1,
    chunks = NULL,
    chunk_size = 1,
    iter = 1000,
    adapt_delta = 0.90,
    max_treedepth = 10,
    seed = 12345,
    init_value = NULL,
    ...) {
  serodata <- validate_prepared_serodata(serodata)
  stopifnot(
    "foi_model must be either `constant`, `tv_normal_log`, or `tv_normal`" =
      foi_model %in% c("constant", "tv_normal_log", "tv_normal"),
    "iter must be numeric" = is.numeric(iter),
    "seed must be numeric" = is.numeric(seed)
  )
  model <- stanmodels[[foi_model]]
  exposure_matrix <- get_exposure_matrix(serodata)
  n_obs <- nrow(serodata)

  if (is.null(chunks)) {
    chunks <- get_chunk_structure(
      serodata = serodata,
      chunk_size = chunk_size
    )
  }
  checkmate::assert_class(chunks, "numeric")
  stopifnot(
    "`chunks` length must be equal to `max(serodata$age_mean_f)`" =
      length(chunks) == max(serodata$age_mean_f)
  )

  stan_data <- list(
    n_obs = n_obs,
    n_pos = serodata$counts,
    n_total = serodata$total,
    age_max = max(serodata$age_mean_f),
    observation_exposure_matrix = exposure_matrix,
    chunks = chunks,
    foi_location = foi_location,
    foi_scale = foi_scale
  )

  if (foi_model == "tv_normal_log") {
    f_init <- function() {
      list(log_fois = rep(-3, max(chunks)))
    }
  } else {
    f_init <- function() {
      list(fois = rep(0.01, max(chunks)))
    }
  }

  seromodel_fit <- rstan::sampling(
    model,
    data = stan_data,
    iter = iter,
    init = f_init,
    control = list(
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth
    ),
    seed = seed,
    # https://github.com/stan-dev/rstan/issues/761#issuecomment-647029649
    chain_id = 0,
    include = FALSE,
    pars = "fois_vector",
    verbose = FALSE,
    refresh = 0,
   ...
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


#' Generate list containing the chunk structure to be used in the retrospective
#' estimation of the force of infection.
#'
#' This function generates a numeric list specifying the chunk structure of the
#' time interval spanning from the year of birth of the oldest age cohort up to
#' the time when the serosurvey was conducted.
#' @inheritParams fit_seromodel
#' @examples
#' data(chagas2012)
#' serodata <- prepare_serodata(serodata = chagas2012, alpha = 0.05)
#' cohort_ages <- get_cohort_ages(serodata = serodata)
#' @export
get_chunk_structure <- function(
  serodata,
  chunk_size
  ) {
    checkmate::assert_int(
      chunk_size,
      lower = 1,
      upper = max(serodata$age_mean_f)
      )

    chunks <- unlist(
      purrr::map(
        seq(
          1,
          max(serodata$age_mean_f) / chunk_size,
          1),
        rep,
        times = chunk_size
      )
    )

    chunks <- append(
      chunks,
      rep(
        max(chunks),
        max(serodata$age_mean_f) - length(chunks)
      )
    )

  return(chunks)
}

#' Generate data frame containing the age of each cohort
#' corresponding to each birth year excluding the year of the survey.
#'
#' This function generates a data frame containing the age of each cohort
#' corresponding to each `birth_year` excluding the year of the survey, for
#' which the cohort age is still 0. specified serological survey data `serodata`
#' excluding the year of the survey.
#' @inheritParams run_seromodel
#' @return `cohort_ages`. A data frame containing the age of each cohort
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

#' Generate exposure matrix corresponding to a serological
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

#' Extract central estimates for the fitted forced FoI
#'
#' @param seromodel_object Stanfit object containing the results of fitting a
#'   model by means of [run_seromodel].
#' @param cohort_ages  A data frame containing the age of each cohort
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
  # extracts force-of-infection from stan fit
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
#'   \item{`n_iter`}{Number of iterations by chain including warmup.}
#'   \item{`elpd`}{elpd}
#'   \item{`se`}{se}
#'   \item{`converged`}{convergence}
#' }
#' @examples
#' data(chagas2012)
#' serodata <- prepare_serodata(chagas2012)
#' seromodel_object <- fit_seromodel(
#'   serodata = serodata,
#'   foi_model = "constant"
#' )
#' extract_seromodel_summary(seromodel_object,
#'   serodata = serodata
#' )
#' @export
extract_seromodel_summary <- function(seromodel_object,
                                      serodata) {
  serodata <- validate_prepared_serodata(serodata)
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


#' Generate data frame containing the confidence interval based on
#' a force-of-infection fitting
#'
#' This function computes the corresponding binomial confidence intervals for
#' the obtained prevalence based on a fitting of the force-of-infection `foi`
#' for plotting an analysis purposes.
#' @param foi Object containing the information of the force-of-infection. It is
#'   obtained from `rstan::extract(seromodel_object$seromodel, "foi", inc_warmup
#'   = FALSE)[[1]]`.
#' @param alpha Probability threshold for statistical significance used for both
#' the binomial confidence interval, and the lower and upper quantiles of the
#' estimated prevalence.
#' @inheritParams run_seromodel
#' @param bin_data If `TRUE`, `serodata` is binned by means of
#'   `prepare_bin_data`. Otherwise, age groups are kept as originally input.
#' @param bin_step Integer specifying the age groups bin size to be used when
#' `bin_data` is set to `TRUE`.
#' @return `prev_final`. The expanded prevalence data. This is used for plotting
#'   purposes in the `visualization` module.
#' @examples
#' data(chagas2012)
#' serodata <- prepare_serodata(chagas2012)
#' seromodel_object <- fit_seromodel(
#'   serodata = serodata,
#'   foi_model = "constant"
#' )
#' foi <- rstan::extract(seromodel_object, "foi")[[1]]
#' get_prev_expanded(foi, serodata)
#' @export
get_prev_expanded <- function(foi,
                              serodata,
                              alpha = 0.05,
                              bin_data = FALSE,
                              bin_step = 5
                              ) {

  if (bin_data && any(serodata$age_max - serodata$age_min > 2)) {
    warning("Make sure `serodata` is already grouped by age")
    bin_data <- FALSE
  }

  foi_expanded <- foi

  ly <- NCOL(foi_expanded)
  exposure_expanded <- matrix(0, nrow = ly, ncol = ly)
  exposure_expanded[apply(
    lower.tri(exposure_expanded, diag = TRUE),
    1, rev
  )] <- 1

  prev_pn <- t(1 - exp(-exposure_expanded %*% t(foi_expanded)))

  predicted_prev <- t(
    apply(
      prev_pn,
      2,
      function(x) {
        quantile(
          x,
          c(
            0.5,
            alpha,
            1 - alpha
          )
        )
      }
    )
  )
  colnames(predicted_prev) <- c(
    "predicted_prev",
    "predicted_prev_lower",
    "predicted_prev_upper"
  )
  predicted_prev <- as.data.frame(predicted_prev)
  predicted_prev$age <- 1:ly

  if (bin_data) {
      observed_prev <- prepare_bin_data(
        serodata = serodata,
        bin_step = bin_step,
        alpha = alpha
      )
      } else {
        observed_prev <- serodata
        }

  observed_prev <- observed_prev %>%
    dplyr::select(
      "age_mean_f",
      "prev_obs",
      "prev_obs_lower",
      "prev_obs_upper",
      "total",
      "counts"
    ) %>%
    rename(
      age = "age_mean_f"
    )

  prev_expanded <-
    base::merge(
      predicted_prev,
      observed_prev,
      by = "age",
      all.x = TRUE
    ) %>%
    dplyr::mutate(survey = unique(serodata$survey))

  return(prev_expanded)
}


#' Fit model to seroprevalence survey using optimization
#' Returns a single best fit parameter vector, rather than
#' samples from the posterior (which is what `fit_seromodel` does)
#' @inheritParams fit_seromodel
#' @param n_iters_max highest number of iterations for optimizer to run
#' @return `seromodel_optimization`. List of best fit optimized parameter value
#' @export
fit_seromodel_optimization <- function(
    serodata,
    foi_model = c("constant", "tv_normal_log", "tv_normal"),
    n_iters_max = 1000,
    seed = 12345,
    ...) {
  
  # Data processing is the same as for fit_seromodel
  # Validate data
  validate_prepared_serodata(serodata)
  stopifnot(
    "foi_model must be either `constant`, `tv_normal_log`, or `tv_normal`" =
      foi_model %in% c("constant", "tv_normal_log", "tv_normal"),
    "seed must be numeric" = is.numeric(seed),
    "n_iters_max must be numeric" = is.numeric(n_iters_max)
  )
  model <- stanmodels[[foi_model]]
  cohort_ages <- get_cohort_ages(serodata = serodata)
  exposure_matrix <- get_exposure_matrix(serodata)
  n_obs <- nrow(serodata)
  
  stan_data <- list(
    n_obs = n_obs,
    n_pos = serodata$counts,
    n_total = serodata$total,
    age_max = max(cohort_ages$age),
    observation_exposure_matrix = exposure_matrix
  )
  
  if (foi_model == "tv_normal_log") {
    f_init <- function() {
      list(log_foi = rep(-3, nrow(cohort_ages)))
    }
  } else {
    f_init <- function() {
      list(foi = rep(0.01, nrow(cohort_ages)))
    }
  }
  
  seromodel_optimization <- rstan::optimizing(
    model,
    data = stan_data,
    iter = n_iters,
    init = f_init,
    verbose = FALSE,
    refresh = 0,
    seed = seed,
    as_vector = FALSE
  )
  
  return(seromodel_optimization)
}

