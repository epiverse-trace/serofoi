#' Get Exposure Matrix
#'
#' Function that gets the exposure matrix
#' @param model_data refers to the model data that has been selected
#' @param yexpo what the make yexpo function returns
#' @return exposure_output
#' @export
get_exposure_matrix <- function(model_data,
                                yexpo) {
  age_class <- model_data$age_mean_f
  ly  <- length(yexpo)
  exposure       <- matrix(0,nrow = length(age_class), ncol = ly)
  for (k in 1:length(age_class)) exposure[k,(ly - age_class[k] + 1):ly] <-  1
  exposure_output <- exposure
  return(exposure_output)

}

#' Get Prevalence Expanded
#'
#' Function that obtains the expanded prevalence
#' @param model_data refers to the model data that has been selected
#' @param foi force of infection
#' @return prev_final
#' @export

get_prev_expanded <- function(foi,
                              model_data) {

  ndata <- data.frame(age = 1:80)
  dim_foi <- dim(foi)[2]
  if (dim_foi < 80)  {
    oldest_year <- 80 - dim_foi + 1
    foin <- matrix(NA, nrow = dim(foi)[1], 80)
    foin[, oldest_year:80] <- foi
    foin[, 1:(oldest_year - 1) ] <- rowMeans(foi[,1:5])
  } else {
    foin <- foi}

  foi_expanded <- foin

  age_class <- 1:NCOL(foi_expanded)
  ly  <- NCOL(foi_expanded)
  exposure       <- matrix(0, nrow = length(age_class), ncol = ly)
  for (k in 1:length(age_class)) exposure[k,(ly - age_class[k] + 1):ly] <-  1
  exposure_expanded <- exposure

  iterf <- NROW(foi_expanded)
  age_max <- NROW(exposure_expanded)
  PrevPn <- matrix(NA, nrow = iterf, ncol = age_max)
  for (i in 1:iterf) {
    PrevPn[i,] <- 1 - exp( -exposure_expanded %*% foi_expanded[i,])
  }

  lower <- apply(PrevPn, 2, function(x) quantile(x, 0.1))
  upper <- apply(PrevPn, 2, function(x) quantile(x, 0.9))
  medianv  <- apply(PrevPn, 2, function(x) quantile(x, 0.5))

  predicted_prev <- data.frame(age = 1:80,
                               predicted_prev = medianv,
                               predicted_prev_lower = lower,
                               predicted_prev_upper = upper)

  observed_prev <- model_data %>%
    dplyr::select(age_mean_f, prev_obs, prev_obs_lower, prev_obs_upper, total, counts) %>%
    dplyr::rename(age = age_mean_f, sample_by_age = total, positives = counts)

  prev_expanded <- base::merge(predicted_prev, observed_prev, by = "age", all.x = TRUE) %>% dplyr::mutate(survey = model_data$survey[1])

  # I added this here for those cases when binned is prefered for plotting
  if (model_data$age_max[1] - model_data$age_min[1] < 3) {
    model_data$cut_ages <- cut(as.numeric(model_data$age_mean_f), seq(1,101, by = 5), include.lowest = TRUE)
    xx <- model_data %>% dplyr::group_by(.data$cut_ages) %>% dplyr::summarise(bin_size = sum(.data$total), bin_pos = sum(.data$counts))
    labs <- read.table(text = gsub("[^.0-9]", " ", levels(xx$cut_ages)), col.names = c("lower", "upper")) %>%
      dplyr::mutate(lev = levels(xx$cut_ages), midAge = round((lower + upper)/2)) %>% dplyr::select(.data$midAge, .data$lev)
    xx$midAge <- labs$midAge[labs$lev %in% xx$cut_ages]
    conf <- data.frame(Hmisc::binconf(xx$bin_pos, xx$bin_size,method = "exact"))
    xx  <- cbind(xx, conf) %>% dplyr::rename(age = .data$midAge, p_obs_bin = .data$PointEst,
                                      p_obs_bin_l = .data$Lower,p_obs_bin_u = .data$Upper)

    prev_final <- base::merge(prev_expanded, xx, by = "age", all.x = TRUE)

  } else {

    prev_final <- prev_expanded %>% dplyr::mutate(cut_ages = "original", bin_size = .data$sample_by_age,
                                            bin_pos  = .data$positives, p_obs_bin = .data$prev_obs,
                                            p_obs_bin_l = .data$prev_obs_lower, p_obs_bin_u = .data$prev_obs_upper)

  }

  return(prev_final)

}

#' Make yexpo
#'
#' Function thats make Yexpo
#' @param model_data refers to the model data that has been selected
#' @return yexpo
#' @export

make_yexpo <- function(model_data) {
  yexpo <- (seq_along(min(model_data$birth_year):model_data$tsur[1]))
}

#' Get Posterior Summary
#'
#' Function that gets a posterior summary
#' @param model_objects model_objects_chain
#' @return model_object
#' @export
get_posterior_summary <- function(model_objects_chain) {
  model_object <- sapply(model_objects_chain,
                function(i) c(quantile(i, c(0.5, 0.025, 0.975))))
  row.names(model_object) <- c("Median", "Lower", "Upper")
  return(model_object)
}

#' Obtain Prevalence Extended
#'
#' Function that obtains the extended prevalence
#' @param model_data refers to the model data that has been selected
#' @param exposure
#' @param ly
#' @param nbreaks
#' @param lambdaYexpo
#' @return new_PPP
#' @export
obtain_prevalence_extended <- function(model_data,
                                       exposure,
                                       ly,
                                       nbreaks,
                                       lambdaYexpo) {

  ly            <- length(yexpo)

  if (nbreaks == 0 & ly < 100)
  {
    ly = 100
    lambdaYexpo <- matrix(lambdaYexpo, nrow = length(lambdaYexpo), ncol = ly )
  }

  new_data       <- data.frame(age = 1:ly)
  exposure_new  <- matrix(0,nrow = length(new_data$age), ncol = ly)
  for (k in 1:length(new_data$age)) exposure_new[k,(ly - new_data$age[k] + 1):ly] <-  1

  olders <- data.frame(age = (ly + 1):99)
  exposure_olders <- matrix(1,nrow = length(olders$age), ncol = ly)

  exposure_total <- rbind(exposure_new, exposure_olders)
  new_age_clases <- c(new_dat$age,   olders$age)

  iterf <- nrow(lambdaYexpo)
  PrevPn <- matrix(NA, nrow = iterf, ncol = length(new_age_clases))
  for (i in 1:iterf) {
    PrevPn[i,] <- 1 - exp( -exposure_total %*% lambdaYexpo[i,])}
  new_PPP <- matrix(NA, ncol = 3, nrow = length(new_age_clases))
  for (j in seq_along(new_age_clases)) {
    new_PPP[j,] <- quantile(PrevPn[,j], c(.025, .5, .975))}
  new_PPP        <- as.data.frame(new_PPP)
  new_PPP$age    <- new_age_clases
  names(new_PPP) <- c("L", "M", "U", "age")

  return(new_PPP)

}

#' Make Thin Chain
#'
#' Function that does the thinning of the number of iterations
#' @param model_object_chain
#' @param thin by default the value 10 is taken but it can be changed
#' @return model_objects_chain
#' @export
make_thin_chain <- function(model_objects_chain, thin = 10)
{
  model_objects_chain <- model_objects_chain[seq(1, nrow(model_objects_chain), thin),]
  return(model_objects_chain)
}

#' Get Residuals
#'
#' Function that gets the residuals
#' @param model_data refers to the model data that has been selected
#' @param fit refers to fit of the model
#' @return merged_prev
#' @export
get_residuals <- function(fit, model_data)
{
  P_sim <- rstan::extract(fit, "P_sim")[[1]]
  colnames(P_sim) <- dat$age_mean_f
  P_sim <- as.data.frame(P_sim)
  P_sim$iteration <- seq_along(P_sim[,1])
  P_sim <- P_sim %>%
    reshape2::melt(id.vars = "iteration") %>%
    dplyr::rename(prev_predicted = .data$value, age = .data$variable)
  P_sim <- P_sim %>% dplyr::mutate(age = as.numeric(as.character(.data$age)))
  P_obs <- dat %>% dplyr::select(age = .data$age_mean_f, .data$prev_obs)

  merged_prev <-  merge(P_sim, P_obs, by = "age") %>%
    dplyr::mutate(residuals = .data$prev_predicted - .data$prev_obs)

  return(merged_prev)
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
#' @export
fit_model <- function(model,
                      model_data,
                      model_name,
                      n_iters = n_iters,
                      n_thin = n_thin,
                      delta = delta,
                      m_treed = m_treed,
                      decades = decades){

  # add a warning to the warming process because there are exceptions where a minimal amount of iterations need to be run
  n_warmup = floor(n_iters/2)

  yexpo <- make_yexpo(model_data)
  yexpo <- yexpo[-length(yexpo)]
  real_yexpo <- (min(model_data$birth_year):model_data$tsur[1])[-1]
  exposure_matrix <- get_exposure_matrix(model_data, yexpo)
  Nobs <- nrow(model_data)

  stan_data <- list(
    Nobs   = Nobs,
    Npos   = model_data$counts,
    Ntotal = model_data$total,
    Age    = model_data$age_mean_f,
    Ymax   = max(yexpo),
    AgeExpoMatrix = exposure_matrix,
    NDecades = decades)

  f_init <- function(){
    list(foi = rep(0.01, length(yexpo)))
  }

  fit <- rstan::sampling(model,
                  data = stan_data ,
                  iter = n_iters,
                  chains = 4,
                  init = f_init,
                  warmup =  n_warmup,
                  verbose = FALSE,
                  refresh = 0,
                  control = list(adapt_delta = delta,
                               max_treedepth = m_treed),
                  seed = "12345",
                  thin = n_thin
  )

  if (class(fit@sim$samples)  != "NULL") {
    loo_fit <- loo::loo(fit, save_psis = TRUE, "logLikelihood")
    foi <- rstan::extract(fit, "foi", inc_warmup = FALSE)[[1]]

    # generates central estimations
    foi_cent_est <- data.frame(year  = real_yexpo,
                               lower = apply(foi, 2, function(x) quantile(x, 0.05)),
                               upper = apply(foi, 2, function(x) quantile(x, 0.95)),
                               medianv = apply(foi, 2, function(x) quantile(x, 0.5)))


    # generates a sample of iterations
    if (n_iters >= 2000){
      foi_post_s <- dplyr::sample_n(as.data.frame(foi), size = 1000)
      colnames(foi_post_s) <- real_yexpo
    }
    else{
      foi_post_s <- as.data.frame(foi)
      colnames(foi_post_s) <- real_yexpo
    }

    model_object <- list(fit = fit,
                         model_data = model_data,
                         stan_data = stan_data,
                         real_yexpo = real_yexpo,
                         yexpo     = yexpo,
                         n_iters   = n_iters,
                         n_thin    = n_thin,
                         n_warmup  = n_warmup,
                         model     = model_name,
                         delta     = delta,
                         m_treed    = m_treed,
                         loo_fit   = loo_fit,
                         foi_cent_est   = foi_cent_est,
                         foi_post_s  = foi_post_s)
  model_object$model_summary <- extract_summary_model(model_object)
  } else {
    loo_fit <- c(-1e10, 0)
    model_object <- list(fit = "no model",
                        model_data = model_data,
                        stan_data = stan_data,
                        real_yexpo = real_yexpo,
                        yexpo     = yexpo,
                        n_iters   = n_iters,
                        n_thin    = n_thin,
                        n_warmup  = n_warmup,
                        model     = model_name,
                        delta     = delta,
                        m_treed    = m_treed,
                        loo_fit   = loo_fit,
                        model_summary = NA)
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
#' @export
fit_model_log <- function(model,
                          model_data,
                          model_name,
                          n_iters = 3000,
                          n_thin = 2,
                          delta = 0.90,
                          m_treed = 10,
                          decades = 0){

  yexpo <- make_yexpo(model_data)
  yexpo <- yexpo[-length(yexpo)]
  real_yexpo <- (min(model_data$birth_year):model_data$tsur[1])[-1]
  exposure_matrix <- get_exposure_matrix(model_data, yexpo)
  Nobs <- nrow(model_data)

  stan_data <- list(
    Nobs   = nrow(model_data),
    Npos   = model_data$counts,
    Ntotal = model_data$total,
    Age    = model_data$age_mean_f,
    Ymax   = max(yexpo),
    AgeExpoMatrix = exposure_matrix,
    NDecades = decades)

  n_warmup = floor(n_iters/2)

  f_init <- function(){
    list(log_foi = rep(-3, length(yexpo)))
  }

  fit <- rstan::sampling(model,
                  data = stan_data ,
                  iter = n_iters,
                  chains = 4,
                  init = f_init,
                  warmup =  n_warmup,
                  verbose = FALSE,
                  refresh = 0,
                  control = list(adapt_delta = delta,
                               max_treedepth = m_treed),
                  seed = "12345",
                  thin = n_thin
  )

  if (class(fit@sim$samples)  != "NULL") {

    loo_fit <- loo::loo(fit, save_psis = TRUE, "logLikelihood")
    foi <- rstan::extract(fit, "foi", inc_warmup = FALSE)[[1]]

    foi_cent_est <- data.frame(year  = real_yexpo,
                               lower = apply(foi, 2, function(x) quantile(x, 0.1)),
                               upper = apply(foi, 2, function(x) quantile(x, 0.9)),
                               medianv = apply(foi, 2, function(x) quantile(x, 0.5)))

    foi_post_s <- dplyr::sample_n(as.data.frame(foi), size = 1000)
    colnames(foi_post_s) <- real_yexpo


    model_object <- list(fit = fit,
                        model_data = model_data,
                        stan_data = stan_data,
                        real_yexpo = real_yexpo,
                        yexpo     = yexpo,
                        n_iters   = n_iters,
                        n_thin    = n_thin,
                        n_warmup  = n_warmup,
                        model     = model_name,
                        delta     = delta,
                        m_treed    = m_treed,
                        loo_fit   = loo_fit,
                        foi_cent_est   = foi_cent_est,
                        foi_post_s  = foi_post_s)

    model_object$model_summary <- extract_summary_model(model_object)
  } else {

    loo_fit <- c(-1e10, 0)
    model_object <- list(fit = "no model",
                          stan_data = stan_data,
                          real_yexpo = real_yexpo,
                          yexpo     = yexpo,
                          n_iters   = n_iters,
                          n_thin    = n_thin,
                          n_warmup  = n_warmup,
                          model     = model_name,
                          delta     = delta,
                          m_treed    = m_treed,
                          loo_fit   = loo_fit)
    model_object$model_summary <- extract_summary_model(model_object)
  }

  return(model_object)

}

#' Save or read model
#' Function that saves the .RDS file of the model
#' @param model_name name of the model selected
save_or_read_model <- function(model_name="constant_foi_bi") {

  rds_path <- config::get(model_name)$rds_path
  stan_path <- config::get(model_name)$stan_path

  if (!file.exists(rds_path)){
    model <- rstan::stan_model(stan_path)
    saveRDS(model, rds_path)
  }
  else{
    model <- readRDS(rds_path)
  }
  return(model)
}


#' Run model
#' Function that runs the specified model
#' @param model_data model_data
#' @param survey survey
#' @param model refers to model selected
#' @param model_name name of the model selected
#' @param n_iters number of iterations. Each model has a value by default.
#' @return model_object of model
#' @export
run_model <- function(model_data,
                      model_name="constant_foi_bi",
                      n_iters=1000,
                      n_thin = 2,
                      delta = 0.90,
                      m_treed = 10,
                      decades = 0) {

  my_dir <- paste0(config::get("test_files_path"), epitrix::clean_labels(paste0("tests_", Sys.time())))

  model_data <- model_data %>% dplyr::arrange(.data$age_mean_f) %>% dplyr::mutate(birth_year = .data$tsur - .data$age_mean_f)
  # survey <- model_data$survey[1] # Revisar la mejor opción para el warning de número de surveys
  survey <- unique(model_data$survey)
  if (length(survey)>0) warning("WARNING!! You have more than 1 surveys or survey codes")

  if (model_name == "constant_foi_bi"){
    model_0 <- save_or_read_model(model_name = model_name)
    model_object <- fit_model(model = model_0,
                              model_data = model_data,
                              model_name = model_name,
                              n_iters = n_iters,
                              n_thin = n_thin,
                              delta = delta,
                              m_treed = m_treed,
                              decades = decades); print(paste0(survey, "finished ------ model_0"))
  }
  if (model_name == "continuous_foi_normal_bi"){
    model_1 <- save_or_read_model(model_name = model_name)
    model_object <- fit_model(model = model_1,
                              model_data = model_data,
                              model_name = model_name,
                              n_iters = n_iters,
                              n_thin = n_thin,
                              delta = delta,
                              m_treed = m_treed,
                              decades = decades); print(paste0(survey, "finished ------ model_1"))
  }
  if (model_name == "continuous_foi_normal_log"){
    model_2   <- save_or_read_model(model_name = model_name)
    model_object <- fit_model_log(model = model_2,
                                  model_data = model_data,
                                  model_name = model_name,
                                  n_iters = n_iters,
                                  n_thin = n_thin,
                                  delta = delta,
                                  m_treed = m_treed,
                                  decades = decades); print(paste0(survey, "finished ------ model_2"))
  }
  print(model_object$model_summary)
  return(model_object)
}

#' Extract Summary Model
#'
#' Function that summarizes the models
#' @param model_object what the run model function returns
#' @param model_data refers to data of the model
#' @return summary of the models
#' @export
extract_summary_model <- function(model_object) {

  model_name <- model_object$model
  #------- Loo estimates

  loo_fit <- model_object$loo_fit
  if (sum(is.na(loo_fit)) < 1)
  {
    lll <- as.numeric((round(loo_fit$estimates[1,],2)))} else
    {
      lll <- c(-1e10, 0)
    }

  model_data <- model_object$model_data
  summary_model <- data.frame(model = model_object$model,
                              dataset = model_data$survey[1],
                              country = model_data$country[1],
                              year    = model_data$tsur[1],
                              test    = model_data$test[1],
                              antibody = model_data$antibody[1],
                              n_sample = sum(model_data$total),
                              n_agec  = length(model_data$age_mean_f),
                              n_iter  = model_object$n_iters,
                              performance = "_____",
                              elpd = lll[1],
                              se = lll[2],
                              converged = NA
  )

  rhats <- get_table_rhats(model_object)
  if (any(rhats$rhat > 1.1 ) == FALSE) {
    summary_model$converged = "Yes"  }

  return(summary_model)
}
