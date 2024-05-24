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

  int<lower=0, upper=1> serorev_prior;
  real<lower=0> serorev_a;
  real<lower=0> serorev_b;
}

transformed data {
  int n_chunks = max(chunks);
}

parameters {
  vector<lower=0>[n_chunks] foi_vector;
  real<lower=0> sigma;
  real<lower=0> seroreversion_rate;
}

transformed parameters {
  vector[n_obs] prob_infected;

  for (i in 1:n_obs){
    int age = ages[i];
    prob_infected[i] = prob_infected_age_varying(
      age,
      foi_vector,
      chunks,
      seroreversion_rate
    );
  }
}

model {
  n_pos ~ binomial(n_total, prob_infected);
  sigma ~ cauchy(0, 1);

	foi_vector[1] ~ normal(foi_location, foi_scale);
  for(i in 2:n_chunks)
    foi_vector[i] ~ normal(foi_vector[i - 1], sigma);

  if(serorev_prior == 0) {
    seroreversion_rate ~  uniform(serorev_a, serorev_b);
  } else if( serorev_prior == 1) {
    seroreversion_rate ~  normal(serorev_a, serorev_b);
  }
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
      seroreversion_rate
    );
  }

  for(i in 1:n_obs){
    n_pos_sim[i] = binomial_rng(n_total[i], prob_infected[i]);
    P_sim[i] = n_pos_sim[i] / n_total[i];
    logLikelihood[i] = binomial_lpmf(n_pos[i] | n_total[i], prob_infected[i]);
  }
}
