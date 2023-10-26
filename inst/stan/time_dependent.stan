functions {

  vector  prob_infected_noseroreversion(
    vector foi, 
    matrix observation_exposure_matrix, 
    int n_obs, 
    int age_max, 
    int[] chunks) {
      real scalar_dot_product;
      vector[n_obs] prob_infected;
      vector[age_max] foi_every_age = foi[chunks];
    
      for(i in 1:n_obs){
        scalar_dot_product = dot_product(observation_exposure_matrix[i,], foi_every_age);
        prob_infected[i] = 1 - exp(-scalar_dot_product);
      }
    
      return prob_infected;
  }
  
  real solve_ode_one_year(
    real prob_previous, 
    real foi, 
    real seroreversion_rate) {
      real e_upper = exp(foi + seroreversion_rate);
      real e_lower = exp(-(foi + seroreversion_rate));
      
      return e_lower * 
        ((prob_previous - 1 + e_upper) * foi + 
          prob_previous * seroreversion_rate) / 
        (foi + seroreversion_rate);
  }
  
  real prob_infected_seroreversion_single_age(
    vector fois_per_year_for_one_age_group, 
    real seroreversion_rate) {
      int n_iter = num_elements(fois_per_year_for_one_age_group);
      real prob_infected = 0;
    
      for(i in 1:n_iter) {
        prob_infected = solve_ode_one_year(prob_infected, 
          fois_per_year_for_one_age_group[i], 
          seroreversion_rate);
      }
    
    return prob_infected;
  }
  
  vector get_fois(
    int i, 
    vector fois_stacked, 
    int n_fois, 
    int [] foi_index_start_per_obs) {
      int start = foi_index_start_per_obs[i];
      vector[n_fois] fois = segment(fois_stacked, start, n_fois);
      
      return fois;
  }
  
  vector prob_infected_seroreversion(
    vector foi, 
    real seroreversion_rate,
    int n_obs,
    int age_max,
    int[] chunks,
    int [] n_fois_exposed_per_obs,
    int [] foi_index_start_per_obs,
    int [] foi_indices,
    int n_fois_exposed) {
      vector[n_obs] prob_infected;
      vector[age_max] fois_every_age = foi[chunks];
      vector[n_fois_exposed] fois_stacked=fois_every_age[foi_indices];
      
      for(i in 1:n_obs) {
        int n_fois = n_fois_exposed_per_obs[i];
        vector[n_fois] fois_per_year_for_one_age_group = get_fois(
          i, fois_stacked, n_fois, foi_index_start_per_obs);
        prob_infected[i] = prob_infected_seroreversion_single_age(
          fois_per_year_for_one_age_group, seroreversion_rate);
      }
      
      return prob_infected;
  }
  
  vector prob_infected_calculate(vector foi, 
    real seroreversion_rate,
    int include_seroreversion, 
    matrix observation_exposure_matrix,
    int n_obs, 
    int age_max, 
    int[] chunks, 
    int [] n_fois_exposed_per_obs,
    int [] foi_index_start_per_obs, 
    int [] foi_indices, 
    int n_fois_exposed) {
      vector[n_obs] prob_infected;
      
      if(include_seroreversion == 1)
        prob_infected = prob_infected_seroreversion(
            foi, seroreversion_rate, n_obs, age_max, chunks,
            n_fois_exposed_per_obs, foi_index_start_per_obs,
            foi_indices, n_fois_exposed);
      else
        prob_infected = prob_infected_noseroreversion(
            foi, observation_exposure_matrix, n_obs, age_max, chunks);
        
      return prob_infected;
    }
}

data {
     int<lower=0> n_obs;
     int n_pos[n_obs];
     int n_total[n_obs];
     int<lower=1> age_max;
     matrix[n_obs, age_max] observation_exposure_matrix; // only used for non-seroreverting models
     int<lower=1> n_fois_exposed_per_obs[n_obs]; // the number of fois each age group is exposed to
     int<lower=1> foi_index_start_per_obs[n_obs]; // this gives the starting point of the foi corresponding to that age group when stacking the foi vectors for all ages on top of one another
     int n_fois_exposed; // total number of foi-years across all observations
     int<lower=1> foi_indices[n_fois_exposed]; // stacked vector of indices giving location of foi within list, corresponding to fois for each age group
     
     // model type
     int<lower=0, upper=1> include_seroreversion;

     // prior choices
     int chunks[age_max];
     int<lower=1, upper=7> prior_choice;
     real prior_a; #<lower=0>
     real prior_b; #<lower=0>
}

transformed data {
  int is_random_walk = prior_choice <= 4 ? 1 : 0;
  int n_chunks = max(chunks);
  matrix[n_obs, age_max] foi_index_matrix;
}

parameters {
   row_vector[n_chunks] log_foi; // log foi
   real<lower=0> sigma[is_random_walk ? 1 : 0]; // normal distribution scale parameter
   real<lower=0> nu[is_random_walk ? 1 : 0]; // student-t degrees of freedom
   real<lower=0> seroreversion_rate[include_seroreversion ? 1 : 0]; // rate of seroreversion
}

transformed parameters {
  vector[n_chunks] foi = to_vector(exp(log_foi));
  real mu;
  vector[n_obs] prob_infected;
  
  if(include_seroreversion)
    mu = seroreversion_rate[1];
  else
    mu = 0;

  prob_infected = prob_infected_calculate(
    foi, mu, include_seroreversion,
    observation_exposure_matrix, n_obs, age_max,
    chunks, n_fois_exposed_per_obs, foi_index_start_per_obs,
    foi_indices, n_fois_exposed);
}

model {
  // likelihood
  n_pos ~ binomial(n_total, prob_infected) ;

  // priors
  if(include_seroreversion)
    seroreversion_rate ~ cauchy(0, 1);
  
  if(prior_choice == 1) { // forward random walk == tv normal log
    sigma ~ cauchy(0, 1);
    log_foi[1] ~ normal(prior_a, prior_b);

    for(i in 2:n_chunks)
      log_foi[i] ~ normal(log_foi[i - 1], sigma);

  } else if (prior_choice == 2) { // backward random walk
    sigma ~ cauchy(0, 1);
    log_foi[n_chunks] ~ normal(prior_a, prior_b);

    for(i in 1:(n_chunks - 1))
      log_foi[n_chunks - i] ~ normal(log_foi[n_chunks - i + 1], sigma);

  } else if(prior_choice == 3){ // forward random walk with Student-t
    sigma ~ cauchy(0, 1);
    nu ~ cauchy(0, 1);
    log_foi[1] ~ normal(prior_a, prior_b);

    for(i in 2:n_chunks)
      log_foi[i] ~ student_t(nu, log_foi[i - 1], sigma);

  } else if (prior_choice == 4) { // backward random walk with Student t
    sigma ~ cauchy(0, 1);
    nu ~ cauchy(0, 1);
    log_foi[n_chunks] ~ normal(prior_a, prior_b);

    for(i in 1:(n_chunks - 1))
      log_foi[n_chunks - i] ~ student_t(nu, log_foi[n_chunks - i + 1], sigma);

  } else if(prior_choice == 5){ // uniform
    foi ~ uniform(prior_a, prior_b);
    target += sum(foi);

  } else if(prior_choice == 6) { // weakly informative
    foi ~ cauchy(prior_a, prior_b);
    target += sum(foi);

 } else if(prior_choice == 7) { // Laplace (sparsity-inducing)
    foi ~ double_exponential(prior_a, prior_b);
    target += sum(foi);
 
 } else if(prior_choice == 8) { // forward random walk - foi
    sigma ~ cauchy(0, 1);
    foi[1] ~ normal(prior_a, prior_b);
    for(i in 2:n_chunks)
      foi[i] ~ normal(foi[i - 1], sigma);
 }
}

generated quantities {
  vector[age_max] fois_by_year;
  for(i in 1:age_max) {
    fois_by_year[i] = foi[chunks[i]];
  }
}
