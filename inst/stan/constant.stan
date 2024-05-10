data {
     int<lower=0> n_obs;
     int<lower=1> age_max;
     int n_pos[n_obs];
     int n_total[n_obs];
     matrix[n_obs, age_max] observation_exposure_matrix;

     // prior choices
     real<lower=0> foi_a;
     real<lower=0, upper=2> foi_b;
}

parameters {
   real<lower=0,upper=2> lambda0;
}

transformed parameters {
  real prob_infected[n_obs];
  real scalar_dot_product[n_obs];
  row_vector[age_max] foi;
  
 for (j in 1:age_max) {
   foi[j] = lambda0;
   }
  
 for (i in 1:n_obs){
   scalar_dot_product[i] = dot_product(observation_exposure_matrix[i,], foi);
   prob_infected[i] = 1 - exp(-scalar_dot_product[i]);
 }
}

model {
  for (i in 1:n_obs)
   n_pos[i] ~ binomial(n_total[i], prob_infected[i]);
   lambda0 ~ uniform (foi_a, foi_b);
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
