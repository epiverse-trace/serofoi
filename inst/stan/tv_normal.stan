data {
     int<lower=0> n_obs;
     int n_pos[n_obs];
     int n_total[n_obs];
     int <lower=1>age_max;
     matrix[n_obs, age_max] observation_exposure_matrix;
}

parameters {
   row_vector<lower=0>[age_max] foi;
   real<lower=0> sigma;
}

transformed parameters {
  real prob_infected[n_obs];
  real scalar_dot_product[n_obs];
  
 for (i in 1:n_obs){
   scalar_dot_product[i] = dot_product(observation_exposure_matrix[i,], foi);
   prob_infected[i] = 1 - exp(-scalar_dot_product[i]);
 }
}

model {
  for (i in 1:n_obs)
    n_pos[i] ~ binomial(n_total[i], prob_infected[i]) ;
    sigma ~ cauchy(0, 1);
  
  for(i in 2:age_max)
    foi[i] ~ normal(foi[i - 1], sigma);
  foi[1] ~ normal(0, 1);
}


generated quantities{
  vector[n_obs] n_pos_sim;
  vector[n_obs] P_sim;
  vector[n_obs] logLikelihood;
  for(i in 1:n_obs){
    n_pos_sim[i] = binomial_rng(n_total[i], prob_infected[i]);
    P_sim[i] = n_pos_sim[i] / n_total[i];
    logLikelihood[i] = binomial_lpmf(n_pos[i] | n_total[i], prob_infected[i]);
  }
}
