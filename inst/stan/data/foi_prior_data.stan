  // prior index
  int foi_prior_index;
  // foi indexes (chunks)
  int foi_index[age_max];
  // uniform
  real<lower=0> foi_min;
  real<lower=foi_min> foi_max;
  // normal
  real<lower=0> foi_mean;
  real<lower=0> foi_sd;
  // cauchy
  real<lower=0> foi_sigma_rw_loc;
  real<lower=0> foi_sigma_rw_sc;
