
#------------------------- Fit Function

fFitModel <- function(model, dat, m_name, 
                      n_iters = 3000, 
                      # n_warmup = 1000,
                      n_thin = 2, 
                      delta = 0.90, 
                      mtreed = 10, 
                      Decades = 0){
  
  # Temporary solution for non convergence models
  # delta = 0.92
  # mtreed = 12
  
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
    list(foi=rep(0.01, length(yexpo)))
  } 
  
  
  # n_warmup = floor(n_iters/2)
  
  
  
  fit <- sampling(model, 
                  data=stan_data , 
                  iter=n_iters,
                  chains=4,
                  init=fInit,
                  warmup =  n_warmup,
                  verbose=FALSE,
                  refresh=0,
                  control=list(adapt_delta= delta,
                               max_treedepth = mtreed),
                  seed = "12345",
                  thin = n_thin
  )
  
  
  if(class(fit@sim$samples)  != "NULL") {
    
    loo_fit <- loo(fit, save_psis = TRUE, 'logLikelihood')
    foi <- rstan::extract(fit, 'foi', inc_warmup = FALSE)[[1]]
    
    foi_cent_est <- data.frame(year  = RealYexpo,
                               lower = apply(foi, 2, function(x) quantile(x, 0.1)),
                               upper = apply(foi, 2, function(x) quantile(x, 0.9)),
                               medianv = apply(foi, 2, function(x) quantile(x, 0.5)))
    
    foi_post_1000s <- dplyr::sample_n(as.data.frame(foi), size = 1000)
    colnames(foi_post_1000s) <- RealYexpo
    
    
    res <- list(fit=fit,
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
    res <- list(fit='no model',
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






# ---- log
# 
fFitModel_log <- function(model, dat, m_name, n_iters = 3000, 
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
    list(log_foi=rep(-3, length(yexpo)))
  } 
  
  
  n_warmup = floor(n_iters/2)
  # n_warmup = 1000
  
  fit <- sampling(model, 
                  data=stan_data , 
                  iter=n_iters,
                  chains=4,
                  init=fInit,
                  warmup =  n_warmup,
                  verbose=FALSE,
                  refresh = 0,
                  control=list(adapt_delta= delta,
                               max_treedepth = mtreed),
                  seed = "12345",
                  thin = n_thin
  )
  
  
  if(class(fit@sim$samples)  != "NULL") {
    
    loo_fit <- loo(fit, save_psis = TRUE, 'logLikelihood')
    foi <- rstan::extract(fit, 'foi', inc_warmup = FALSE)[[1]]
    
    foi_cent_est <- data.frame(year  = RealYexpo,
                               lower = apply(foi, 2, function(x) quantile(x, 0.1)),
                               upper = apply(foi, 2, function(x) quantile(x, 0.9)),
                               medianv = apply(foi, 2, function(x) quantile(x, 0.5)))
    
    foi_post_1000s <- dplyr::sample_n(as.data.frame(foi), size = 1000)
    colnames(foi_post_1000s) <- RealYexpo
    
    
    res <- list(fit=fit,
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
    res <- list(fit='no model',
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


