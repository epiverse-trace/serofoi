functions {
  #include functions/prob_infected_av.stan
}

data {
	int<lower=0> n_obs;
	int<lower=1> age_max;
	int n_pos[n_obs];
	int n_total[n_obs];
	int ages[n_obs];

	// prior choices
	int chunks[age_max];
	real<lower=0> foi_location;
	real<lower=0> foi_scale;
}

transformed data {
  int n_chunks = max(chunks);
}

parameters {
  vector<lower=0>[n_chunks] foi_vector;
  real<lower=0> sigma;
}

transformed parameters {
  vector[n_obs] prob_infected;

  for (i in 1:n_obs){
    int age = ages[i];
    prob_infected[i] = prob_infected_age_varying(
      age,
      foi_vector,
      chunks,
      0.0
    );
  }
}

model {
  n_pos ~ binomial(n_total, prob_infected);
  sigma ~ cauchy(0, 1);

	foi_vector[1] ~ normal(foi_location, foi_scale);
  for(i in 2:n_chunks)
    foi_vector[i] ~ normal(foi_vector[i - 1], sigma);
}

generated quantities{
  vector[n_obs] n_pos_sim;
  vector[n_obs] P_sim;
  vector[n_obs] logLikelihood;
  vector[age_max] foi;
  vector[age_max] prob_infected_expanded;

  for(age in 1:age_max) {
    foi[age] = foi_vector[chunks[age]];
    prob_infected_expanded[age] = prob_infected_age_varying(
      age,
      foi_vector,
      chunks,
      0.0
    );
  }

  for(i in 1:n_obs){
    n_pos_sim[i] = binomial_rng(n_total[i], prob_infected[i]);
    P_sim[i] = n_pos_sim[i] / n_total[i];
    logLikelihood[i] = binomial_lpmf(n_pos[i] | n_total[i], prob_infected[i]);
  }
}
