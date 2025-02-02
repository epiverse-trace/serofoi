data {
     int<lower=0> n_obs;
     array[n_obs] int n_pos;
     array[n_obs] int n_total;
     int<lower=1> age_max;
     matrix[n_obs, age_max] observation_exposure_matrix;
}

parameters {
   vector[age_max] log_foi;
   real<lower=0> sigma;
}

transformed parameters {
  vector[n_obs] prob_infected;
  vector[n_obs] scalar_dot_product;
  vector<lower=0>[age_max] foi;
  foi = exp(log_foi);

 for (i in 1:n_obs){
   scalar_dot_product[i] = dot_product(observation_exposure_matrix[i,], foi);
   prob_infected[i] = 1 - exp(-scalar_dot_product[i]);
 }

}


model {
    n_pos ~ binomial(n_total, prob_infected);
    sigma ~ cauchy(0, 1);

  for(i in 2:age_max)
    log_foi[i] ~ normal(log_foi[i - 1], sigma);
    log_foi[1] ~ normal(-6, 4);
 }


generated quantities{
  vector[n_obs] n_pos_sim;
  vector[n_obs] P_sim;
  vector[n_obs] logLikelihood;
  for (i in 1:n_obs) {
    n_pos_sim[i] = binomial_rng(n_total[i], prob_infected[i]);
    P_sim[i] = n_pos_sim[i] / n_total[i];
    logLikelihood[i] = binomial_lpmf(n_pos[i] | n_total[i], prob_infected[i]);
  }
}
