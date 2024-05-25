  real prob_infected_age_varying(
    int age,
    real age_max,
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
