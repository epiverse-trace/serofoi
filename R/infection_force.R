#' Dir Results
#'
#' Función guarda los resultados en una dirección
#' Function saves the results in a place
#' @param name_dir
#' @return dir
#' @export
get_exposure_matrix <- function(dat, yexpo) {
  age_class <- dat$age_mean_f
  ly  <- length(yexpo)
  exposure       <- matrix(0,nrow = length(age_class), ncol = ly)
  for (k in 1:length(age_class)) exposure[k,(ly - age_class[k] + 1):ly] <-  1
  exposure_output <- exposure
  return(exposure_output)

}

#' Get Prevalence Expanded
#'
#' Función que obtiene la prevalencia expandida
#' Function that obtains the expanded prevalence
#' @param data dat
#' @param foi
#' @return prev_final
#' @export
get_prev_expanded <- function(foi, dat) {

  ndat <- data.frame(age = 1:80)
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


  observed_prev <- dat %>%
    dplyr::select(age_mean_f, prev_obs, prev_obs_lower, prev_obs_upper, total, counts) %>%
    dplyr::rename(age = age_mean_f, sample_by_age = total, positives = counts)

  prev_expanded <- merge(predicted_prev, observed_prev, by = 'age', all.x = TRUE) %>% dplyr::mutate(survey = dat$survey[1])

  # I added this here for those cases when binned is prefered for plotting
  if (dat$age_max[1] - dat$age_min[1] < 3) {
    dat$cut_ages <- cut(as.numeric(dat$age_mean_f), seq(1,101, by = 5), include.lowest = TRUE)
    xx <- dat %>% dplyr::group_by(.data$cut_ages) %>% dplyr::summarise(bin_size = sum(.data$total), bin_pos = sum(.data$counts))
    labs <- read.table(text = gsub("[^.0-9]", " ", levels(xx$cut_ages)), col.names = c("lower", "upper")) %>%
      dplyr::mutate(lev = levels(xx$cut_ages), midAge = round((lower + upper)/2)) %>% dplyr::select(.data$midAge, .data$lev)
    xx$midAge <- labs$midAge[labs$lev %in% xx$cut_ages]
    conf <- data.frame(Hmisc::binconf(xx$bin_pos, xx$bin_size,method = "exact"))
    xx  <- cbind(xx, conf) %>% dplyr::rename(age = .data$midAge, p_obs_bin = .data$PointEst,
                                      p_obs_bin_l = .data$Lower,p_obs_bin_u = .data$Upper)

    prev_final <- merge(prev_expanded, xx, by = 'age', all.x = TRUE)

  } else {

    prev_final <- prev_expanded %>% dplyr::mutate(cut_ages = "original", bin_size = .data$sample_by_age,
                                            bin_pos  = .data$positives, p_obs_bin = .data$prev_obs,
                                            p_obs_bin_l = .data$prev_obs_lower, p_obs_bin_u = .data$prev_obs_upper)

  }

  return(prev_final)

}

#' Make Yexpo
#'
#' Función que hace Yexpo
#' Function thats make Yexpo
#' @param data dat
#' @return Yexpo
#' @export
make_yexpo <- function(dat) {
  yexpo <- (seq_along(min(dat$birth_year):dat$tsur[1]))
}

#' Get Posterior Summary
#'
#' Función que consigue un resumen posterior
#' Function that gets a post summary
#' @param results results_chain
#' @return results
#' @export
get_posterior_summary <- function(results_chain) {
  res <- sapply(results_chain,
                function(i) c(quantile(i, c(0.5, 0.025, 0.975))))
  row.names(res) <- c('Median', 'Lower', 'Upper')
  return(res)
}

#' Obtain Prevalence Extended
#'
#' Función que obtiene la prevalencia extendida
#' Function that obtains the extended prevalence
#' @param data dat0
#' @param exposure
#' @param ly
#' @param nbreaks
#' @param lambdaYexpo
#' @return new_PPP
#' @export
obtain_prevalence_extended <- function(dato, exposure, ly, nbreaks, lambdaYexpo) {

  ly            <- length(yexpo)

  if (nbreaks == 0 & ly < 100)
  {
    ly = 100
    lambdaYexpo <- matrix(lambdaYexpo, nrow = length(lambdaYexpo), ncol = ly )
  }

  new_dat       <- data.frame(age = 1:ly)
  exposure_new  <- matrix(0,nrow = length(new_dat$age), ncol = ly)
  for (k in 1:length(new_dat$age)) exposure_new[k,(ly - new_dat$age[k] + 1):ly] <-  1

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
#' Función que hace el adelgazamiento del número de iteraciones
#' Function that does the thinning of the number of iterations
#' @param result_chain
#' @param thin
#' @return results_chain
#' @export
make_thin_chain <- function(results_chain, thin = 10)
{
  results_chain <- results_chain[seq(1, nrow(results_chain), thin),]
  return(results_chain)
}

#' Get Residuals
#'
#' Función que obtiene los residuos
#' Function that gets the residuals
#' @param data dat
#' @param fit
#' @return merged_prev
#' @export
get_residuals <- function(fit, dat)
{
  P_sim <- rstan::extract(fit, 'P_sim')[[1]]
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
#' Función que ajusta el modelo a los datos
#' Function that fits the model to the data
#' @param data dat
#' @param model
#' @param m_name
#' @param n_iters
#' @param n_thin
#' @param delta
#' @param mtreed
#' @param mDecades
#' @return res
#' @export
fit_model <- function(model, dat, m_name,
                      n_iters = 3000,
                      n_thin = 2,
                      delta = 0.90,
                      mtreed = 10,
                      Decades = 0){


  n_warmup = floor(n_iters/2)


  yexpo <- make_yexpo(dat)
  yexpo <- yexpo[-length(yexpo)]
  RealYexpo <- (min(dat$birth_year):dat$tsur[1])[-1]
  ExposureMatrix <- get_exposure_matrix(dat, yexpo)
  Nobs <- nrow(dat)


  stan_data <- list(
    Nobs   = nrow(dat),
    Npos   = dat$counts,
    Ntotal = dat$total,
    Age    = dat$age_mean_f,
    Ymax   = max(yexpo),
    AgeExpoMatrix = ExposureMatrix,
    NDecades = Decades)

  fInit <- function(){
    list(foi = rep(0.01, length(yexpo)))
  }

  fit <- rstan::sampling(model,
                  data = stan_data ,
                  iter = n_iters,
                  chains = 4,
                  init = fInit,
                  warmup =  n_warmup,
                  verbose = FALSE,
                  refresh = 0,
                  control = list(adapt_delta = delta,
                               max_treedepth = mtreed),
                  seed = "12345",
                  thin = n_thin
  )


  if (class(fit@sim$samples)  != "NULL") {
    loo_fit <- loo::loo(fit, save_psis = TRUE, 'logLikelihood')
    foi <- rstan::extract(fit, 'foi', inc_warmup = FALSE)[[1]]

    foi_cent_est <- data.frame(year  = RealYexpo,
                               lower = apply(foi, 2, function(x) quantile(x, 0.1)),
                               upper = apply(foi, 2, function(x) quantile(x, 0.9)),
                               medianv = apply(foi, 2, function(x) quantile(x, 0.5)))

    foi_post_1000s <- dplyr::sample_n(as.data.frame(foi), size = 1000)
    colnames(foi_post_1000s) <- RealYexpo

    res <- list(fit = fit,
                stan_data = stan_data,
                RealYexpo = RealYexpo,
                yexpo     = yexpo,
                n_iters   = n_iters,
                n_thin    = n_thin,
                n_warmup  = n_warmup,
                model     = m_name,
                delta     = delta,
                mtreed    = mtreed,
                loo_fit   = loo_fit,
                foi_cent_est   = foi_cent_est,
                foi_post_1000s  = foi_post_1000s)


  } else {

    loo_fit <- c(-1e10, 0)
    res <- list(fit = "no model",
                stan_data = stan_data,
                RealYexpo = RealYexpo,
                yexpo     = yexpo,
                n_iters   = n_iters,
                n_thin    = n_thin,
                n_warmup  = n_warmup,
                model     = m_name,
                delta     = delta,
                mtreed    = mtreed,
                loo_fit   = loo_fit)

  }

  return(res)

}

#' Fit Model Log
#'
#' Función que ajusta el modelo logarítmico a los datos
#' Function that fits the logarithmic model to the data
#' @param data dat
#' @param model
#' @param m_name
#' @param n_iters
#' @param n_thin
#' @param delta
#' @param mtreed
#' @param mDecades
#' @return res
#' @export
fit_model_log <- function(model, dat, m_name, n_iters = 3000,
                          n_thin = 2, delta = 0.90, mtreed = 10, Decades = 0){

  yexpo <- make_yexpo(dat)
  yexpo <- yexpo[-length(yexpo)]
  RealYexpo <- (min(dat$birth_year):dat$tsur[1])[-1]
  ExposureMatrix <- get_exposure_matrix(dat, yexpo)
  Nobs <- nrow(dat)

  stan_data <- list(
    Nobs   = nrow(dat),
    Npos   = dat$counts,
    Ntotal = dat$total,
    Age    = dat$age_mean_f,
    Ymax   = max(yexpo),
    AgeExpoMatrix = ExposureMatrix,
    NDecades = Decades)

  fInit <- function(){
    list(log_foi = rep(-3, length(yexpo)))
  }

  n_warmup = floor(n_iters/2)


  fit <- rstan::sampling(model,
                  data = stan_data ,
                  iter = n_iters,
                  chains = 4,
                  init = fInit,
                  warmup =  n_warmup,
                  verbose = FALSE,
                  refresh = 0,
                  control = list(adapt_delta = delta,
                               max_treedepth = mtreed),
                  seed = "12345",
                  thin = n_thin
  )

  if (class(fit@sim$samples)  != "NULL") {

    loo_fit <- loo::loo(fit, save_psis = TRUE, 'logLikelihood')
    foi <- rstan::extract(fit, 'foi', inc_warmup = FALSE)[[1]]

    foi_cent_est <- data.frame(year  = RealYexpo,
                               lower = apply(foi, 2, function(x) quantile(x, 0.1)),
                               upper = apply(foi, 2, function(x) quantile(x, 0.9)),
                               medianv = apply(foi, 2, function(x) quantile(x, 0.5)))

    foi_post_1000s <- dplyr::sample_n(as.data.frame(foi), size = 1000)
    colnames(foi_post_1000s) <- RealYexpo


    res <- list(fit = fit,
                stan_data = stan_data,
                RealYexpo = RealYexpo,
                yexpo     = yexpo,
                n_iters   = n_iters,
                n_thin    = n_thin,

                n_warmup  = n_warmup,
                model     = m_name,
                delta     = delta,
                mtreed    = mtreed,
                loo_fit   = loo_fit,
                foi_cent_est   = foi_cent_est,
                foi_post_1000s  = foi_post_1000s)


  } else {

    loo_fit <- c(-1e10, 0)
    res <- list(fit = "no model",
                stan_data = stan_data,
                RealYexpo = RealYexpo,
                yexpo     = yexpo,
                n_iters   = n_iters,
                n_thin    = n_thin,

                n_warmup  = n_warmup,
                model     = m_name,
                delta     = delta,
                mtreed    = mtreed,
                loo_fit   = loo_fit)

  }

  return(res)

}