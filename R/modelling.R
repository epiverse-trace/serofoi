#' Get Exposure Matrix
#'
#' Function that generates the exposure matrix
#' @param model_data refers to the model data that has been selected
#' @param yexpo what the make yexpo function returns
#' @return exposure_output
#' @examples
#' get_exposure_matrix <- function(model_data, yexpo)
#' @export
get_exposure_matrix <- function(model_data,
                                yexpo) {
  age_class <- model_data$age_mean_f
  ly <- length(yexpo)
  exposure <- matrix(0, nrow = length(age_class), ncol = ly)
  for (k in 1:length(age_class))
    exposure[k, (ly - age_class[k] + 1):ly] <- 1
  exposure_output <- exposure
  return(exposure_output)
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

#' Make yexpo
#'
#' Function that generates Yexpo
#' @param model_data refers to the model data that has been selected
#' @return yexpo
#' @examples
#' make_yexpo (model_data)
#' @export
make_yexpo <- function(model_data) {
  yexpo <- (seq_along(min(model_data$birth_year):model_data$tsur[1]))
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

#' Fit Model
#'
#' Function that fits the model to the data
#' @param model_data refers to the model data that has been selected
#' @param model refers to model selected
#' @param model_name name of the model selected
#' @param n_iters number of iterations. Each model has a default number.
#' @param n_thin Each model has a default number.
#' @param delta This value comes by default but it can be changed
#' @param m_treed This value comes by default but it can be changed
#' @param decades The decades covered by the survey data
#' @return model_object
#' @examples
#' fit_model (model,
#'            model_data,
#'            model_name,
#'            n_iters = n_iters,
#'            n_thin = n_thin,
#'            delta = delta,
#'            m_treed = m_treed,
#'            decades = decades)
#' @export
fit_model <- function(model,
                      model_data,
                      model_name,
                      n_iters = n_iters,
                      n_thin = n_thin,
                      delta = delta,
                      m_treed = m_treed,
                      decades = decades) {
  # add a warning to the warming process because there are exceptions where a minimal amount of iterations need to be run
  n_warmup <- floor(n_iters / 2)

  yexpo <- make_yexpo(model_data)
  yexpo <- yexpo[-length(yexpo)]
  real_yexpo <- (min(model_data$birth_year):model_data$tsur[1])[-1]
  exposure_matrix <- get_exposure_matrix(model_data, yexpo)
  Nobs <- nrow(model_data)

  stan_data <- list(
    Nobs = Nobs,
    Npos = model_data$counts,
    Ntotal = model_data$total,
    Age = model_data$age_mean_f,
    Ymax = max(yexpo),
    AgeExpoMatrix = exposure_matrix,
    NDecades = decades
  )

  f_init <- function() {
    list(foi = rep(0.01, length(yexpo)))
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
      year = real_yexpo,
      lower = apply(foi, 2, function(x) quantile(x, 0.05)),

      upper = apply(foi, 2, function(x) quantile(x, 0.95)),

      medianv = apply(foi, 2, function(x) quantile(x, 0.5))
    )


    # generates a sample of iterations
    if (n_iters >= 2000) {
      foi_post_s <- dplyr::sample_n(as.data.frame(foi), size = 1000)
      colnames(foi_post_s) <- real_yexpo
    } else {
      foi_post_s <- as.data.frame(foi)
      colnames(foi_post_s) <- real_yexpo
    }

    model_object <- list(
      fit = fit,
      model_data = model_data,
      stan_data = stan_data,
      real_yexpo = real_yexpo,
      yexpo = yexpo,
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
      extract_summary_model(model_object)
  } else {
    loo_fit <- c(-1e10, 0)
    model_object <- list(
      fit = "no model",
      model_data = model_data,
      stan_data = stan_data,
      real_yexpo = real_yexpo,
      yexpo = yexpo,
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
fit_model_log <- function(model,
                          model_data,
                          model_name,
                          n_iters = 3000,
                          n_thin = 2,
                          delta = 0.90,
                          m_treed = 10,
                          decades = 0) {
  yexpo <- make_yexpo(model_data)
  yexpo <- yexpo[-length(yexpo)]
  real_yexpo <- (min(model_data$birth_year):model_data$tsur[1])[-1]
  exposure_matrix <- get_exposure_matrix(model_data, yexpo)
  Nobs <- nrow(model_data)

  stan_data <- list(
    Nobs = nrow(model_data),
    Npos = model_data$counts,
    Ntotal = model_data$total,
    Age = model_data$age_mean_f,
    Ymax = max(yexpo),
    AgeExpoMatrix = exposure_matrix,
    NDecades = decades
  )

  n_warmup <- floor(n_iters / 2)

  f_init <- function() {
    list(log_foi = rep(-3, length(yexpo)))
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
      year = real_yexpo,
      lower = apply(foi, 2, function(x) quantile(x, 0.1)),

      upper = apply(foi, 2, function(x) quantile(x, 0.9)),

      medianv = apply(foi, 2, function(x) quantile(x, 0.5))
    )

    foi_post_s <- dplyr::sample_n(as.data.frame(foi), size = 1000)
    colnames(foi_post_s) <- real_yexpo


    model_object <- list(
      fit = fit,
      model_data = model_data,
      stan_data = stan_data,
      real_yexpo = real_yexpo,
      yexpo = yexpo,
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
      extract_summary_model(model_object)
  } else {
    loo_fit <- c(-1e10, 0)
    model_object <- list(
      fit = "no model",
      stan_data = stan_data,
      real_yexpo = real_yexpo,
      yexpo = yexpo,
      n_iters = n_iters,
      n_thin = n_thin,
      n_warmup = n_warmup,
      model = model_name,
      delta = delta,
      m_treed = m_treed,
      loo_fit = loo_fit
    )
    model_object$model_summary <-
      extract_summary_model(model_object)
  }

  return(model_object)
}

#' Save or read model
#' Function that saves the .RDS file of the model
#' @param model_name name of the model selected
#' @return model
#' @examples
#' save_or_read_model (model_name = "constant_foi_bi")
#' @export
save_or_read_model <- function(model_name = "constant_foi_bi") {
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


#' Run model
#' Runs a pre-specified stan model for the force-of-infection
#' @param model_data a dataframe after using <prepare_data> function
#' @param model_name name of the model selected. Current version provides three options:
#' <constant_foi_bi>, <continuous_foi_normal_bi> and <continuous_foi_normal_log>
#' @param n_iters number of iterations. There is a default value by model type.
#' @return model_object which is a list of various objects including stan objects and others
#' @examples
#' run_model (model_data,
#'            model_name = "constant_foi_bi",
#' @export
run_model <- function(model_data,
                      model_name = "constant_foi_bi",
                      n_iters = 1000,
                      n_thin = 2,
                      delta = 0.90,
                      m_treed = 10,
                      decades = 0) {
  model_data <- model_data %>%
    dplyr::arrange(.data$age_mean_f) %>%
    dplyr::mutate(birth_year = .data$tsur - .data$age_mean_f)
  survey <- unique(model_data$survey)
  if (length(survey) > 1) warning("You have more than 1 surveys or survey codes")

  if (model_name == "constant_foi_bi") {
    model_0 <- save_or_read_model(model_name = model_name)
    model_object <- fit_model(model = model_0,
                              model_data = model_data,
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
    model_1 <- save_or_read_model(model_name = model_name)
    model_object <- fit_model(model = model_1,
                              model_data = model_data,
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
    model_2 <- save_or_read_model(model_name = model_name)
    model_object <- fit_model_log(model = model_2,
                                  model_data = model_data,
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

#' Extract Summary Model
#'
#' Function that summarizes the models
#' @param model_object what the run model function returns
#' @param model_data refers to data of the model
#' @return summary of the models
#' @examples
#' extract_summary_model (model_object)
#' @export
extract_summary_model <- function(model_object) {
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
