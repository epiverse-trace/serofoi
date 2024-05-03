  vector  prob_infected_noseroreversion(
    vector fois_vector,
    matrix observation_exposure_matrix,
    int n_obs,
    int age_max,
    int[] chunks
    ) {
      real scalar_dot_product;
      vector[n_obs] prob_infected;
      vector[age_max] foi_every_age = fois_vector[chunks];

      for(i in 1:n_obs){
        scalar_dot_product = dot_product(
          observation_exposure_matrix[i,],
          foi_every_age
          );
        prob_infected[i] = 1 - exp(-scalar_dot_product);
      }

      return prob_infected;
  }

  vector prob_infected_calculate(
    vector fois_vector,
    matrix observation_exposure_matrix,
    int n_obs,
    int age_max,
    int[] chunks
    ) {
      vector[n_obs] prob_infected;

      prob_infected = prob_infected_noseroreversion(
        fois_vector,
        observation_exposure_matrix,
        n_obs,
        age_max,
        chunks
        );

      return prob_infected;
    }
