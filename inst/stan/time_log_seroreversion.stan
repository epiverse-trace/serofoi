functions {
  #include functions/prob_infected_time.stan
}

data {
  #include data/basic_data.stan
  #include data/foi_prior_data.stan
  #include data/seroreversion_prior_data.stan
}

transformed data {
  int n_foi = max(foi_index);
}

parameters {
  vector[n_foi] log_foi_vector;
  real<lower=0> seroreversion_rate;
  real<lower=0> sigma;
}

transformed parameters {
  vector<lower=0>[n_foi] foi_vector;
  vector[n_observations] prob_infected;

  foi_vector = exp(log_foi_vector);

  prob_infected = prob_infected_time_model(
    age_groups,
  	n_observations,
  	age_max,
  	foi_vector,
  	foi_index,
  	seroreversion_rate
  );
}

model {
  n_seropositive ~ binomial(n_sample, prob_infected);
  sigma ~ normal(foi_sigma_rw_loc, foi_sigma_rw_sc);

  // force of infection prior
  if (foi_prior_index == 0)
    foi_vector[1] ~ uniform(foi_min, foi_max);
  if (foi_prior_index == 1)
    foi_vector[1] ~ normal(foi_mean, foi_sd);
  target += log(exp(log_foi_vector[1])); // Jacobian term

  for(i in 2:n_foi)
    log_foi_vector[i] ~ normal(log_foi_vector[i - 1], sigma);

  // seroreversion prior
  if (seroreversion_prior_index == 0)
    seroreversion_rate ~ uniform(seroreversion_min, seroreversion_max);
  if (seroreversion_prior_index == 1)
    seroreversion_rate ~ normal(seroreversion_mean, seroreversion_sd);
}

generated quantities{
  #include generated_quantities/log_likelihood.stan

  vector[age_max] prob_infected_expanded;
  vector[age_max] foi_expanded;

  for(i in 1:age_max) {
    int time_index = age_max - i + 1;
    foi_expanded[time_index] = foi_vector[foi_index[i]];
  }

	prob_infected_expanded = prob_infected_time_model(
    ages,
  	age_max, // number of ages
  	age_max, // number of years
  	foi_vector,
  	foi_index,
  	seroreversion_rate
  );
}
