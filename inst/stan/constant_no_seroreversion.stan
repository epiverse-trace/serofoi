functions {
  #include functions/prob_infected_constant.stan
}
data {
  #include data/basic_data.stan
  #include data/foi_prior_data.stan
}

parameters {
   real<lower=0> foi;
}

transformed parameters {
  vector[n_observations] prob_infected;

  prob_infected = prob_infected_constant_model(
    age_groups,
		n_observations,
    foi,
    0.0
  );
}

model {
  n_seropositive ~ binomial(n_sample, prob_infected);

  // force of infection prior
  if (foi_prior_index == 0)
    foi ~ uniform(foi_min, foi_max);
  if (foi_prior_index == 1)
    foi ~ normal(foi_mean, foi_sd);
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
		0.0
	);
}
