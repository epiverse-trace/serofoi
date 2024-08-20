vector[n_observations] log_likelihood;

for(i in 1:n_observations){
  log_likelihood[i] = binomial_lpmf(n_seropositive[i] | sample_size[i], prob_infected[i]);
}
