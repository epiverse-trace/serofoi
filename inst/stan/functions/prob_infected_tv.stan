vector  prob_infected_noseroreversion(
  vector foi_vector,
  int[] chunks,
  matrix observation_exposure_matrix,
  int n_obs,
  int age_max
  ) {
    real scalar_dot_product;
    vector[n_obs] prob_infected;
    vector[age_max] foi_every_age = foi_vector[chunks];

    for(i in 1:n_obs){
      scalar_dot_product = dot_product(
        observation_exposure_matrix[i,],
        foi_every_age
        );
      prob_infected[i] = 1 - exp(-scalar_dot_product);
    }

    return prob_infected;
}

  real prob_infected_age_varying(
    int age,
    vector foi_vector,
    int[] chunks,
    real seroreversion_rate
  ) {
    real prob = 0.0;
    for(j in 1:age){
      real foi = foi_vector[chunks[j]];
      real lambda_over_both = foi / (foi + seroreversion_rate);
      real e_lower = exp(-(foi + seroreversion_rate));

      prob = lambda_over_both + e_lower * (prob - lambda_over_both);
    }
    return prob;
  }
