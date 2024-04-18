functions {
  #include functions/prob_infected.stan
}

data {
  int<lower=0> n_obs;
  int n_pos[n_obs];
  int n_total[n_obs];
  int <lower=1>age_max;
  matrix[n_obs, age_max] observation_exposure_matrix;

  // prior choices
  int chunks[age_max];
  real<lower=0> foi_location;
  real<lower=0> foi_scale;
}

transformed data {
  int n_chunks = max(chunks);
}

parameters {
   row_vector<lower=0>[n_chunks] fois;
   real<lower=0> sigma;
}

transformed parameters {
  vector[n_chunks] fois_vector;
  vector[n_obs] prob_infected;

  fois_vector = to_vector(fois);

  prob_infected = prob_infected_calculate(
    fois_vector,
    observation_exposure_matrix,
    n_obs,
    age_max,
    chunks
  );
}

model {
  n_pos ~ binomial(n_total, prob_infected);
  sigma ~ cauchy(0, 1);

  fois[1] ~ normal(foi_location, foi_scale);
  for(i in 2:n_chunks)
    fois[i] ~ normal(fois[i - 1], sigma);
}

generated quantities{
  vector[n_obs] n_pos_sim;
  vector[n_obs] P_sim;
  vector[n_obs] logLikelihood;
  vector[age_max] foi;

  for(i in 1:age_max) {
    foi[i] = fois_vector[chunks[i]];
  }

  for(i in 1:n_obs){
    n_pos_sim[i] = binomial_rng(n_total[i], prob_infected[i]);
    P_sim[i] = n_pos_sim[i] / n_total[i];
    logLikelihood[i] = binomial_lpmf(n_pos[i] | n_total[i], prob_infected[i]);
  }
}
