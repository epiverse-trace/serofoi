real prob_infected_time_model_single_age(
  int age,
  int age_max,
  vector foi_vector,
  int[] foi_index,
  real seroreversion_rate
) {
  real prob = 0.0;
  int birth_year = age_max - age;
  for(j in 1:age) {
    real foi = foi_vector[foi_index[birth_year + j]];
    real lambda_over_both = foi / (foi + seroreversion_rate);
    real e_lower = exp(-(foi + seroreversion_rate));

    prob = lambda_over_both + e_lower * (prob - lambda_over_both);
  }
  return prob;
}

vector prob_infected_time_model(
  int[] ages,
	int n_ages,
	int age_max,
	vector foi_vector,
	int[] foi_index,
	real seroreversion_rate
) {
	vector[n_ages] prob_infected;

  for (i in 1:n_ages) {
    int age = ages[i];
    prob_infected[i] = prob_infected_time_model_single_age(
      age,
      age_max,
      foi_vector,
      foi_index,
      seroreversion_rate
    );
  }
  return prob_infected;
}
