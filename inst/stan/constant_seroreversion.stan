functions {
  #include functions/prob_infected_constant.stan
}

data {
  #include data/basic_data.stan
  #include data/foi_prior_data.stan
  #include data/seroreversion_prior_data.stan
}

parameters {
   real<lower=0> foi;
   real<lower=0> seroreversion_rate;
}

transformed parameters {
  vector[n_observations] prob_infected;

  prob_infected = prob_infected_constant_model(
    age_groups,
		n_observations,
    foi,
    seroreversion_rate
  );
}

model {
  n_seropositive ~ binomial(sample_size, prob_infected);

  // force of infection prior
  if (foi_prior_index == 0)
    foi ~ uniform(foi_min, foi_max);
  if (foi_prior_index == 1)
    foi ~ normal(foi_mean, foi_sd);

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
    foi_expanded[i] = foi;
  }

	prob_infected_expanded = prob_infected_constant_model(
		ages,
		age_max,
		foi,
		seroreversion_rate
	);
}
