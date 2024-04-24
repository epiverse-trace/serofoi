real prob_infected_age_varying(
  vector fois_vector,
  int[] chunks,
  int age
) {
  real prob = 0.0;
  for(j in 1:age){
      real foi = fois_vector[chunks[j]];
      prob = 1/foi * exp(-foi) * (foi * (exp(foi)  - 1) + prob * foi);
  }
  return prob;
}

vector prob_infected_calculate(
  vector fois_vector,
  int[] chunks,
  int[] ages,
  int n_obs
) {
  vector[n_obs] prob_infected;

  for(i in 1:n_obs){
    int age = ages[i];
    prob_infected[i] = prob_infected_age_varying(
      fois_vector,
      chunks,
      age
      );
  }

  return prob_infected;
}
