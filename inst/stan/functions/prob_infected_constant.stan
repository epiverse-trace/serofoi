real prob_infected_constant_single_age(
    int age,
    real foi,
    real seroreversion_rate
) {
    real foi_over_both = foi / (foi + seroreversion_rate);
    real e_lower = exp(-(foi + seroreversion_rate));

    real prob = 0.0;
    for(i in 1:age) { 
    	prob = foi_over_both + e_lower * (prob - foi_over_both);
    }
    return prob;
}

vector prob_infected_constant(
	int[] ages,
	int n_ages,
	real foi,
	real seroreversion_rate
) {
	vector[n_ages] prob_infected;

	for (i in 1:n_ages) {
		prob_infected[i] = prob_infected_constant_single_age(
			ages[i],
			foi,
			seroreversion_rate
		);
	}
	return prob_infected;
}
