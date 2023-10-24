### This script tests the combined models posted on issue 69
library(tidyverse)
library(rstan)
library(serofoi)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

set.seed(1234)
age_max <- 80

#### SIMULATED DATA ####
# generate synthetic data assuming we have data for all individuals aged between 1-age_max
year <- 2013
ages <- seq(1, age_max, 1)
#using a huge sample size in each age group to check
#that model recovers true values in large sample limit
sample_size <- 100000

# Staircase-like FOI
foi <- c(rep(0.03, age_max / 4), rep(0.05, age_max / 4), rep(0.01, age_max / 2))
mu <- 0.005

# Generate exposure matrix: years exposed depending on age 1=exposed, 0=non-exposed
exposure_matrix_sim <- matrix(nrow = length(ages),
                              ncol = length(ages))
for (i in seq_along(ages)) {
  n_zeros <- length(ages) - i
  n_ones <- i
  exposure_matrix_sim[i, ] <- c(rep(0, n_zeros), rep(1, n_ones))
}

# Map exposure matrix indexes with 1 (exposed years) into vector
foi_indices <- unlist(map(seq(1, nrow(exposure_matrix_sim), 1), ~ which(exposure_matrix_sim[., ] == 1)))
fois_long <- foi[foi_indices] # Map FOI

# Analytical solution to ODE
simulate_foi_exact <- function(a, fois, mu) {
  I <- 0
  # solves ODE exactly within pieces
  for (i in 1:a) {
    lambda <- fois[i]
    I <- (1 / (lambda + mu)) * exp(- (lambda + mu)) * (lambda * (exp(lambda + mu)  - 1) + I * (lambda + mu))
  }
  I
}

# Declare prob_infected:
# analytical prob of infected calculated with exact solution (simulate_foi_exact)
prob_infected <- vector(length = length(ages))

# Indexes for FOIs depending on number of years exposed
n_fois_exposed_per_obs <- rowSums(exposure_matrix_sim)
foi_index_start_per_obs <- c(1, 1 + cumsum(n_fois_exposed_per_obs))
foi_index_start_per_obs <- foi_index_start_per_obs[-length(foi_index_start_per_obs)]

# Calculate Prob. infected using indexes and analytical solution
for (a in seq_along(ages)) {
  start <- foi_index_start_per_obs[a]
  len <- n_fois_exposed_per_obs[a]
  end <- start + len - 1
  fois <- fois_long[start:end] #Reconstruction of foi
  prob_infected[a] <- simulate_foi_exact(a, fois, mu)
}

# Check reconstructed foi
all(fois == foi)

# Plot simulated data
plot(prob_infected)
plot(fois)

# Binomial sample for each age using prob_infected
n_positive_age <- vector(length = length(ages))
for (i in seq_along(n_positive_age)) {
  n_positive_age[i] <- rbinom(1, size = sample_size, prob = prob_infected[i])
}

####################################################################
#### DATA STRUCTURE SEROFOI ####

# Dummy dataframe structured for serofoi
serodata_sim <- data.frame(
  survey = rep("AAAA", length(ages)),
  country = rep("AAAA", length(ages)),
  test = rep("XXXX", length(ages)),
  antibody = rep("IgG", length(ages)),
  total = rep(sample_size, length(ages)),
  counts = n_positive_age,
  age_min = ages,
  age_max = ages,
  tsur = rep(year, length(ages))
)

# Test serofoi EDA functions
plot_seroprev(serodata_sim)
serodata_sim_p <- prepare_serodata(serodata_sim)
exposure_matrix <- get_exposure_matrix(serodata_sim_p)
all(exposure_matrix == exposure_matrix_sim) #Same exposure matrix

cohort_ages <- get_cohort_ages(serodata_sim_p)

#### TODO: Create function to generate this variables

prob_infected <- vector(length = length(serodata_sim_p$age_mean_f))
n_fois_exposed_per_obs <- rowSums(exposure_matrix)
foi_index_start_per_obs <- c(1, 1 + cumsum(n_fois_exposed_per_obs))
foi_index_start_per_obs <- foi_index_start_per_obs[-length(foi_index_start_per_obs)]
foi_indices <- unlist(map(seq(1, nrow(exposure_matrix), 1), ~which(exposure_matrix[., ] == 1)))
fois_long <- foi[foi_indices]
####################################################################

#### TODO: Create function to generate data dictionary
chunks <- unlist(map(seq(1, 8, 1), ~rep(., 10)))
data_stan <- list(
  n_obs = length(serodata_sim_p$counts),
  n_pos = serodata_sim_p$counts,
  n_total = serodata_sim_p$total,
  age_max = max(serodata_sim_p$age_mean_f),
  observation_exposure_matrix = exposure_matrix,
  n_fois_exposed_per_obs = n_fois_exposed_per_obs,
  foi_index_start_per_obs = foi_index_start_per_obs,
  include_seroreversion = 1,
  n_fois_exposed = sum(n_fois_exposed_per_obs),
  foi_indices = foi_indices,
  chunks = chunks,
  prior_choice = 5,
  prior_a = 0,
  prior_b = 5
)

model <- stan_model("inst/stan/time_dependent.stan")
init_fn <- function() {
  list(log_foi = log(rep(0.01, max(data_stan$chunks))))
}
fit <- optimizing(model, data = data_stan, init = init_fn, as_vector = FALSE)

# check fois are similar between true and estimated values
foi_est <- fit$par$foi[chunks]
tibble(estimated=foi_est,
       true=foi) %>% 
  mutate(age=80 - ages) %>% 
  pivot_longer(-age) %>% 
  ggplot(aes(x=age, y=value, colour=name)) +
  geom_line() +
  xlab("years ago")

# check that the estimated probability of infection across ages looks same as true
probs <- fit$par$prob_infected

tibble(true=prob_infected, estimated=probs, age=ages) %>%
  pivot_longer(-age) %>% 
  ggplot(aes(x=age, y=value, colour=name)) +
  geom_line()