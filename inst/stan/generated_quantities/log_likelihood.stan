vector[n_observations] n_seropositive_sim;
vector[n_observations] prob_infected_sim;
vector[n_observations] log_likelihood;

for(i in 1:n_observations){
  n_seropositive_sim[i] = binomial_rng(sample_size[i], prob_infected[i]);
  prob_infected_sim[i] = n_seropositive_sim[i] / sample_size[i];
  log_likelihood[i] = binomial_lpmf(n_seropositive[i] | sample_size[i], prob_infected[i]);
}
