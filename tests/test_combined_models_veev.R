library(serofoi)
library(tidyverse)
library(cowplot)
data("veev2012")

veev2012p <- prepare_serodata(veev2012)

#### RESULTS OBTAINED DIRECTLY USING SEROFOI ####
m1_veev <- run_seromodel(serodata = veev2012p,
                        foi_model = "constant",
                        n_iters = 500,
                        n_thin = 2)

m2_veev <- run_seromodel(serodata = veev2012p,
                        foi_model = "tv_normal",
                        n_iters = 500,
                        n_thin = 2)

m3_veev <- run_seromodel(serodata = veev2012p,
                        foi_model = "tv_normal_log",
                        n_iters = 500,
                        n_thin = 2)

# Visualisation of the results
p1_veev <- plot_seromodel(m1_veev, serodata = veev2012p, size_text = 6)
p2_veev <- plot_seromodel(m2_veev, serodata = veev2012p, size_text = 6)
p3_veev <- plot_seromodel(m3_veev, serodata = veev2012p, size_text = 6)

####################################################################
####################################################################
#### RESULTS USING COMBINED MODELS STAN ####

## Test of time_dependent.stan using tv log normal for veev
foi_model <- "tv_normal_log" #this must be unified with the prior_choice (for tv lognormal use 1)

### 1: forward random walk <==> tv log normal model
### 2: backward random walk
### 3: forward random walk with Student-t
### 4: backward random walk with Student t
### 5: uniform
### 6: weakly informative
### 7: Laplace (sparsity-inducing)
### 8: Normal <==> tv normal model

# Load and prepare data with serofoi functions
serodata_p <- prepare_serodata(veev2012)

cohort_ages <- get_cohort_ages(serodata_p)
exposure_years <- cohort_ages$birth_year
exposure_matrix <- get_exposure_matrix(serodata_p)

# Indexing data
n_fois_exposed_per_obs <- rowSums(exposure_matrix)
foi_index_start_per_obs <- c(1, 1 + cumsum(n_fois_exposed_per_obs))
foi_index_start_per_obs <- foi_index_start_per_obs[-length(foi_index_start_per_obs)]
foi_indices <- map(seq(1, nrow(exposure_matrix), 1), ~which(exposure_matrix[., ] == 1)) %>%
  unlist()

# Stan data dictionary

#chunks <- map(seq(1, 8, 1), ~rep(., 10)) %>%
#  unlist()
chunks <- map(seq(1, 55, 1), ~rep(., 1)) %>%
  unlist()

data_stan <- list(
  n_obs = length(serodata_p$counts),
  n_pos = serodata_p$counts,
  n_total = serodata_p$total,
  age_max = max(serodata_p$age_mean_f),
  observation_exposure_matrix = exposure_matrix,
  n_fois_exposed_per_obs = n_fois_exposed_per_obs,
  foi_index_start_per_obs = foi_index_start_per_obs,
  include_seroreversion = 0,
  n_fois_exposed = sum(n_fois_exposed_per_obs),
  foi_indices = foi_indices,
  chunks = chunks,
  prior_choice = 1,
  prior_a = -6,
  prior_b = 4
)

# Run model with sampling

## Model
model <- rstan::stan_model("inst/stan/time_dependent.stan")

## Feats
n_iters <- 1000
n_thin <- 2
delta <- 0.90
m_treed <- 10
decades <- 0

n_warmup <- floor(n_iters / 2)

## Init
if (foi_model == "tv_normal_log") {
  f_init <- function() {
    list(log_foi = rep(-3, nrow(cohort_ages)))
  }
} else {
  f_init <- function() {
    list(foi = rep(0.01, nrow(cohort_ages)))
  }
}

## Sampling
m3_veev_combined <- rstan::sampling(
  model,
  data = data_stan,
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
  chain_id = 0
)

plt1 <- plot_foi(m3_veev_combined, cohort_ages)
plt2 <- plot_foi(m3_veev, cohort_ages)

plot_grid(plt1, plt2, ncol = 2)