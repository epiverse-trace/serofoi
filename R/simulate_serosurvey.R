
#' Computes the probability of being seropositive when
#' Forces-of-Infection (FoIs) vary by age
#'
#' @param ages Integer indicating the ages of the exposed cohorts
#' @param fois Numeric atomic vector corresponding to the age-varying
#' Force-of-Infection to simulate from
#' @param seroreversion_rate Non-negative seroreversion rate. Default is 0.
#' @return vector of probabilities of being seropositive for age-varying FoI
#' including seroreversion (ordered from youngest to oldest individuals)
#' @export
probability_exact_age_varying <- function(
  ages,
  fois,
  seroreversion_rate = 0
) {

  probabilities <- vector(length = length(ages))
  # solves ODE exactly within pieces
  for (i in seq_along(ages)) {
    foi_tmp <- fois[i]
    if (i == 1)
      probability_previous <- 0
    else
      probability_previous <- probabilities[i - 1]
    if (foi_tmp == 0)
      lambda_over_both <- 0
    else
      lambda_over_both <- (foi_tmp / (foi_tmp + seroreversion_rate))
    probabilities[i] <- lambda_over_both +
      (probability_previous - lambda_over_both) *
      exp(-(foi_tmp + seroreversion_rate))
  }
  return(probabilities)
}

#' Computes the probability of being seropositive when
#' Forces-of-Infection (FoIs) vary by time
#'
#' @param years Integer indicating the years covering the birth ages of the
#' sample
#' @param fois Numeric atomic vector corresponding to the age-varying
#' FoI to simulate from
#' @param seroreversion_rate Non-negative seroreversion rate. Default is 0.
#' @return vector of probabilities of being seropositive for age-varying FoI
#' including seroreversion (ordered from youngest to oldest individuals)
#' @export
probability_exact_time_varying <- function(
  years,
  fois,
  seroreversion_rate = 0
) {

  n_years <- length(years)

  probabilities <- vector(length = length(years))
  # solves ODE exactly within pieces
  for (i in seq_along(years)) { # birth cohorts
    probability_previous <- 0
    for (j in 1:(n_years - i + 1)) { # exposure during lifetime
      foi_tmp <- fois[i + j - 1]
      if (foi_tmp == 0)
        lambda_over_both <- 0
      else
        lambda_over_both <- (foi_tmp / (foi_tmp + seroreversion_rate))
      probability_previous <- lambda_over_both +
        (probability_previous - lambda_over_both) *
        exp(-(foi_tmp + seroreversion_rate))
    }
    probabilities[i] <- probability_previous
  }
  probabilities_oldest_age_last <- rev(probabilities)
  return(probabilities_oldest_age_last)
}

#' Generate probabilities of seropositivity by age based on a time-varying
#' Force-of-Infection (FoI) model.
#'
#' This function calculates the probabilities of seropositivity by age based on
#' a time-varying FoI model.
#' It takes into account the FoI and the rate of seroreversion.
#'
#' @param foi A dataframe containing the FoI values
#' for different years. It should have two columns: 'year' and 'foi'.
#' @param seroreversion_rate A non-negative numeric value representing the
#' rate of seroreversion.
#' @return A dataframe with columns 'age' and 'seropositivity'.
#' @export
prob_seroprev_time_by_age <- function(
  foi,
  seroreversion_rate
) {

  years <- foi$year

  probabilities <- probability_exact_time_varying(
    years = years,
    fois = foi$foi,
    seroreversion_rate = seroreversion_rate
  )

  seroprev_df <- data.frame(
    age = seq_along(years),
    seropositivity = probabilities
  )

  seroprev_df
}


#' Generate probabilities of seropositivity by age based on an age-varying
#' Force-of-Infection (FoI) model.
#'
#' This function calculates the probabilities of seropositivity by age based on
#' an age-varying FoI model.
#' It takes into account the FoI and the rate of seroreversion.
#'
#' @param foi A dataframe containing the FoI values for
#' different ages. It should have two columns: 'age' and 'foi'.
#' @param seroreversion_rate A non-negative numeric value representing the rate
#' of seroreversion.
#' @return A dataframe with columns 'age' and 'seropositivity'.
#' @export
prob_seroprev_age_by_age <- function(
  foi,
  seroreversion_rate
) {

  ages <- seq_along(foi$age)

  probabilities <- probability_exact_age_varying(
    ages = ages,
    fois = foi$foi,
    seroreversion_rate = seroreversion_rate
  )

  seroprev_df <- data.frame(
    age = ages,
    seropositivity = probabilities
  )

  seroprev_df
}

#' Generate probabilities of seropositivity by age based on an age-and-time
#' varying Force-of-Infection (FoI) model.
#'
#' This function calculates the probabilities of seropositivity by age based on
#' an age-and-time-varying FoI model.
#' It takes into account the FoI and the rate of seroreversion.
#'
#' @param foi A dataframe containing the FoI values
#' for different ages. It should have three columns: 'year', 'age' and 'foi'.
#' @param seroreversion_rate A non-negative numeric value representing
#' the rate of seroreversion.
#' @return A dataframe with columns 'age' and 'seropositivity'.
#' @export
prob_seroprev_age_time_by_age <- function(
  foi,
  seroreversion_rate
) {

  foi_matrix <- tidyr::pivot_wider(
    foi,
    values_from = foi,
    names_from = dplyr::all_of("year")
    ) |>
    tibble::column_to_rownames("age") |>
    as.matrix()

  years <- unique(foi$year)
  n_years <- length(years)
  ages <- unique(foi$age)
  n_ages <- length(ages)

  probabilities <- vector(length = length(n_ages))
  # solves ODE exactly within pieces
  for (i in seq_along(years)) { # birth cohorts
    probability_previous <- 0
    foi_matrix_subset <- as.matrix( # only to handle single element matrix case
      foi_matrix[1:(n_ages - i + 1), i:n_ages]
    )
    foi_diag <- diag(foi_matrix_subset)
    for (j in 1:(n_years - i + 1)) { # exposure during lifetime
      foi_tmp <- foi_diag[j]
      if (foi_tmp == 0)
        lambda_over_both <- 0
      else
        lambda_over_both <- (foi_tmp / (foi_tmp + seroreversion_rate))
      probability_previous <- lambda_over_both +
        (probability_previous - lambda_over_both) *
        exp(-(foi_tmp + seroreversion_rate))
    }
    probabilities[i] <- probability_previous
  }
  probabilities_oldest_age_last <- rev(probabilities)

  seroprev_df <- data.frame(
    age = ages,
    seropositivity = probabilities_oldest_age_last
  )

  seroprev_df
}

#' Generate probabilities of seropositivity by age based on model choice.
#'
#' This function generates seropositivity probabilities based on either a
#' time-varying Force-of-Infection (FoI) model, an age-varying FoI model,
#' or an age-and-time-varying FoI model.
#' In all cases, it is possible to optionally include seroreversion.
#'
#' @inheritParams simulate_serosurvey
#' @param seroreversion_rate A non-negative value determining the rate of
#'                           seroreversion (per year). Default is 0.
#'
#' @return A dataframe with columns 'age' and 'seropositivity'.
#' @examples
#' prob_seroprev_by_age(
#'   model = "age",
#'   foi = data.frame(
#'     age = 1:80,
#'     foi = rep(0.01, 80)
#'   )
#' )
#' @export
prob_seroprev_by_age <- function(
  model,
  foi,
  seroreversion_rate = 0
) {

  if (model == "time") {
    probability_function <- prob_seroprev_time_by_age
  } else if (model == "age") {
    probability_function <- prob_seroprev_age_by_age
  } else if (model == "age-time" || model == "time-age") {
    probability_function <- prob_seroprev_age_time_by_age
  }

  seroprev_df <- probability_function(
    foi = foi,
    seroreversion_rate = seroreversion_rate
  )

  return(seroprev_df)
}

#' Generate probabilities of seropositivity by age based on a general
#' Force-of-Infection (FoI) model.
#'
#' This function calculates the probabilities of seropositivity by age based on
#' an abstract model of the serocatalytic system.
#'
#' @param construct_A_fun A function that constructs a matrix that defines the
#' multiplier term in the linear ODE system.
#' @param calculate_seroprev_fun A function which takes the state
#' vector and returns the seropositive fraction.
#' @param initial_conditions The initial state vector proportions for each
#' birth cohort.
#' @param max_age The maximum age to simulate seropositivity for.
#' @param ... Additional parameters for `construct_A_fun`
#' @return A dataframe with columns 'age' and 'seropositivity'.
#' @examples
#' # define age- and time-specific multipliers
#' foi_df_time <- data.frame(
#'   year = seq(1946, 2025, 1),
#'   foi = c(rep(0, 40), rep(1, 40))
#' )
#'
#' foi_df_age <- data.frame(
#'   age = 1:80,
#'   foi = 2 * dlnorm(1:80, meanlog = 3.5, sdlog = 0.5)
#' )
#'
#' u <- foi_df_age$foi
#' v <- foi_df_time$foi
#'
#' # function to construct A matrix for one piece
#' construct_A <- function(t, tau, u, v) {
#'   u_bar <- u[t - tau]
#'   v_bar <- v[t]
#'
#'   A <- diag(-1, ncol = 12, nrow = 12)
#'   A[row(A) == (col(A) + 1)] <- 1
#'   A[1, 1] <- -u_bar * v_bar
#'   A[2, 1] <- u_bar * v_bar
#'   A[12, 12] <- 0
#'
#'   A
#' }
#'
#' # determines the sum of seropositive compartments of those still alive
#' calculate_seropositivity_fn <- function(Y) {
#'   sum(Y[2:11]) / (1 - Y[12])
#' }
#'
#' # initial conditions in 12D state vector
#' initial_conditions <- rep(0, 12)
#' initial_conditions[1] <- 1
#'
#' # calculate probability
#' seropositive_hiv <- prob_seroprev_gen_by_age(
#'   construct_A,
#'   calculate_seropositivity_fn,
#'   initial_conditions,
#'   max_age = 80,
#'   u,
#'   v
#' )
#' @export
prob_seroprev_gen_by_age <- function(
  construct_A_fun,
  calculate_seroprev_fun,
  initial_conditions,
  max_age,
  ...
) {
  # Auxiliar function to compute the sum of A
  sum_of_A <- function(t, tau, construct_A_fun, ...) {
    k <- 1
    for (t_primed in (tau + 1):t) {
      if (k == 1) {
        A <- construct_A_fun(t_primed, tau, ...)
      } else {
        tmp <- construct_A_fun(t_primed, tau, ...)
        A <- A + tmp
      }
      k <- k + 1
    }
    A
  }

  probabilities <- vector(length = max_age)
  for (i in seq_along(probabilities)) {
    A_sum <- sum_of_A(max_age, max_age - i, construct_A_fun, ...)
    Y <- expm::expm(A_sum) %*% initial_conditions
    probabilities[i] <- calculate_seroprev_fun(Y)
  }

  seroprev_df <- data.frame(
    age = 1:max_age,
    seropositivity = probabilities
  )

  seroprev_df
}


#' Add bins based on age intervals.
#'
#' It generates a new column 'group' in the survey_features dataframe,
#' representing the group interval for each row based on the age_min and
#' age_max columns.
#'
#' @param survey_features A dataframe containing age_min and age_max columns
#' representing the minimum and maximum age boundaries for each group.
#'
#' @return A dataframe with an additional 'group' column representing the group
#' interval for each row based on the age_min and age_max columns.
#' @noRd
add_age_bins <- function(survey_features) {
  intervals <- vector(length = nrow(survey_features))
  for (i in seq_along(intervals)) {
    age_min <- survey_features$age_min[i]
    age_max <- survey_features$age_max[i]
    intervals[i] <- paste0("[", age_min, ",", age_max, "]")
  }
  survey_features <- dplyr::mutate(survey_features, group = intervals)

  survey_features
}

#' Create a survey dataframe with per individual age information.
#'
#'
#' @param survey_features A dataframe containing information about age groups
#' and sample sizes.
#' @param age_df A dataframe containing 'age' and 'group'.
#'
#' @return A dataframe with overall sample sizes calculated by joining
#' survey_features and age_df.
#' This dataframe has columns including 'age' and 'overall_sample_size'.
#' @noRd
survey_by_individual_age <- function(survey_features, age_df) {
  overall_sample_size_df <- dplyr::left_join(
      age_df, survey_features,
      by = "group"
    ) |>
    dplyr::rename(overall_sample_size = "n_sample")

  overall_sample_size_df
}

#' Generate random sample sizes using multinomial sampling.
#'
#' This function generates random sample sizes for each age group using
#' multinomial sampling. It takes the total sample size and the number of age
#' groups as input and returns a vector of sample sizes for each age group.
#'
#' @param n_sample The total sample size to be distributed among age groups.
#' @param n_ages The number of age groups.
#'
#' @return A vector containing random sample sizes for each age group.
#' @noRd
multinomial_sampling_group <- function(n_sample, n_ages) {
  prob_value <- 1 / n_ages
  probs <- rep(prob_value, n_ages)
  sample_size_by_age <- as.vector(
    stats::rmultinom(1, n_sample, prob = probs)
  )

  sample_size_by_age
}

#' Generate random sample sizes for each age group.
#'
#' This function generates random sample sizes for each age group based on the
#' overall sample size and the distribution of individuals across age groups.
#' It uses multinomial sampling to allocate the total sample size to each age
#' group proportionally.
#'
#' @param survey_df_long A dataframe with columns 'age', 'group' and
#' 'overall_sample_size'.
#'
#' @return A dataframe with random sample sizes generated for each age based on
#' the overall sample size.
#' @noRd
generate_random_sample_sizes <- function(survey_df_long) {

  df_new <- NULL
  intervals <- unique(survey_df_long$group)
  for (interval_aux in stats::na.omit(intervals)) {
    df_tmp <- dplyr::filter(
      survey_df_long,
      .data$group == interval_aux
    )
    n_sample <- df_tmp$overall_sample_size[1]
    sample_size_by_age <- multinomial_sampling_group(n_sample, nrow(df_tmp))
    df_tmp <- dplyr::mutate(
      df_tmp,
      n_sample = sample_size_by_age
    )

    if (is.null(df_new)) {
      df_new <- df_tmp
    } else {
      df_new <- dplyr::bind_rows(df_new, df_tmp)
    }
  }

  df_new
}

#' Generate random sample sizes for each individual age based on survey
#' features.
#'
#' This function generates random sample sizes for each individual age based on
#' the provided survey features. It first creates age bins, assigns each
#' individual in the survey features to an age bin, calculates the overall
#' sample size by group, and then generates random sample sizes for each age
#' group. Finally, it returns a dataframe with the random sample sizes for each
#' individual age.
#'
#' @param survey_features A dataframe containing information about individuals'
#' age ranges and sample sizes.
#'
#' @return A dataframe with random sample sizes generated for each individual
#' age based on the provided survey features.
#' @noRd
sample_size_by_individual_age <- function(survey_features) {

  ages <- seq(1, max(survey_features$age_max), 1)

  age_bins <- rep(NA, length(ages))
  for (i in seq_along(ages)) {
    # Find index of the row in survey_features where age falls within the range
    age_min <- survey_features$age_min
    age_max <- survey_features$age_max
    idx <- which(ages[i] >= age_min & ages[i] <= age_max)
    if (length(idx) > 0) {
      # Assign the group based on the index found
      age_bins[i] <- paste0("[", age_min[idx], ",", age_max[idx], "]")
    }
  }

  age_df <- data.frame(
    age = ages,
    group = age_bins
    )

  survey_features <- add_age_bins(survey_features)

  survey_features_all_ages <- survey_by_individual_age(
    survey_features,
    age_df)

  df_new <- generate_random_sample_sizes(survey_features_all_ages)

  df_new
}

#' Generate seropositivity counts by bin given the probability and sample size
#' per age group bin
#'
#' @param prob_seroprev_by_age Probability of seropositivity by age
#' @param sample_size_by_age_random Random sample size by age
#' @inheritParams simulate_serosurvey
#' @keywords internal
get_seroprev_counts_by_bin <- function(
    prob_seroprev_by_age,
    sample_size_by_age_random,
    survey_features
) {

  combined_df <- dplyr::left_join(
      prob_seroprev_by_age, sample_size_by_age_random,
      by = "age"
    ) |>
    dplyr::mutate(
      n_seropositive = stats::rbinom(
        nrow(prob_seroprev_by_age),
        .data$n_sample,
        .data$seropositivity)
    )

  grouped_df <- dplyr::group_by(
    combined_df,
    .data$age_min, .data$age_max
  ) |>
  dplyr::summarise(
    n_sample = sum(.data$n_sample),
    n_seropositive = sum(.data$n_seropositive),
    .groups = "drop"
  ) |>
  dplyr::left_join(
    survey_features,
    by = c("age_min", "age_max", "n_sample")
  )

  grouped_df
}

#' Simulate serosurvey data based on a time-varying
#' Force-of-Infection (FoI) model.
#'
#' This function generates binned serosurvey data based on a time-varying FoI
#' model, optionally including seroreversion. This function allows construction
#' of serosurveys with binned age groups, and it generates uncertainty in the
#' distribution of a sample size within an age bin through multinomial sampling.
#'
#' @inheritParams simulate_serosurvey
#' @return A dataframe with simulated serosurvey data, including age group
#' information, overall sample sizes, the number of seropositive individuals,
#' and other survey features.
#' @examples
#' # specify FOIs for each year
#' foi_df <- data.frame(
#'   year = seq(1990, 2009, 1),
#'   foi = rnorm(20, 0.1, 0.01)
#' )
#' survey_features <- data.frame(
#'   age_min = c(1, 3, 15),
#'   age_max = c(2, 14, 20),
#'   n_sample = c(1000, 2000, 1500))
#' serosurvey <- simulate_serosurvey_time(
#' foi_df, survey_features)
#' @export
simulate_serosurvey_time <- function(
  foi,
  survey_features,
  seroreversion_rate = 0
) {

  # Input validation
  validate_foi_df(foi, "year")
  validate_survey_features(survey_features)
  validate_seroreversion_rate(seroreversion_rate)
  validate_simulation_age(
    survey_features,
    foi
    )

  probability_serop_by_age <- prob_seroprev_time_by_age(
    foi = foi,
    seroreversion_rate = seroreversion_rate
  )

  sample_size_by_age_random <- sample_size_by_individual_age(
    survey_features = survey_features
  )

  grouped_df <- get_seroprev_counts_by_bin(
    probability_serop_by_age,
    sample_size_by_age_random,
    survey_features
  )

  return(grouped_df)
}


#' Simulate serosurvey data based on an age-varying
#' Force-of-Infection (FoI) model.
#'
#' This function generates binned serosurvey data based on an age-varying FoI
#' model, optionally including seroreversion. This function allows construction
#' of serosurveys with binned age groups, and it generates uncertainty in the
#' distribution of a sample size within an age bin through multinomial sampling.
#'
#' @inheritParams simulate_serosurvey
#' @return A dataframe with simulated serosurvey data, including age group
#' information, overall sample sizes, the number of seropositive individuals,
#' and other survey features.
#' @examples
#' # specify FOIs for each year
#' foi_df <- data.frame(
#'   age = seq(1, 20, 1),
#'   foi = rnorm(20, 0.1, 0.01)
#' )
#' survey_features <- data.frame(
#'   age_min = c(1, 3, 15),
#'   age_max = c(2, 14, 20),
#'   n_sample = c(1000, 2000, 1500))
#' serosurvey <- simulate_serosurvey_age(
#' foi_df, survey_features)
#' @export
simulate_serosurvey_age <- function(
  foi,
  survey_features,
  seroreversion_rate = 0
) {

  # Input validation
  validate_foi_df(foi, "age")
  validate_survey_features(survey_features)
  validate_seroreversion_rate(seroreversion_rate)
  validate_simulation_age(
    survey_features,
    foi
  )

  probability_serop_by_age <- prob_seroprev_age_by_age(
    foi = foi,
    seroreversion_rate = seroreversion_rate
  )

  sample_size_by_age_random <- sample_size_by_individual_age(
    survey_features = survey_features
  )

  grouped_df <- get_seroprev_counts_by_bin(
    probability_serop_by_age,
    sample_size_by_age_random,
    survey_features
  )

  return(grouped_df)
}

#' Simulate serosurvey data based on an age-and-time-varying
#' Force-of-Infection (FoI) model.
#'
#' This function generates binned serosurvey data based on an
#' age-and-time-varying FoI model, optionally including seroreversion.
#' This function allows construction of serosurveys with binned age groups, and
#' it generates uncertainty in the distribution of a sample size within an age
#' bin through multinomial sampling.
#'
#' @inheritParams simulate_serosurvey
#' @return A dataframe with simulated serosurvey data, including age group
#' information, overall sample sizes, the number of seropositive individuals,
#' and other survey features.
#' @examples
#' # specify FOIs for each year
#' foi_df <- expand.grid(
#'   year = seq(1990, 2009, 1),
#'   age = seq(1, 20, 1)
#' )
#' foi_df$foi <- rnorm(20 * 20, 0.1, 0.01)
#' survey_features <- data.frame(
#'   age_min = c(1, 3, 15),
#'   age_max = c(2, 14, 20),
#'   n_sample = c(1000, 2000, 1500))
#' serosurvey <- simulate_serosurvey_age_time(
#' foi_df, survey_features)
#' @export
simulate_serosurvey_age_time <- function(
  foi,
  survey_features,
  seroreversion_rate = 0
) {

  # Input validation
  validate_foi_df(foi, c("age", "year"))
  validate_survey_features(survey_features)
  validate_seroreversion_rate(seroreversion_rate)
  validate_simulation_age_time(
    survey_features,
    foi
  )

  probability_serop_by_age <-
    prob_seroprev_age_time_by_age(
      foi = foi,
      seroreversion_rate = seroreversion_rate
      )

  sample_size_by_age_random <- sample_size_by_individual_age(
    survey_features = survey_features
  )

  grouped_df <- get_seroprev_counts_by_bin(
    probability_serop_by_age,
    sample_size_by_age_random,
    survey_features
  )

  return(grouped_df)
}


#' Simulate serosurvey data based on various Force-of-Infection (FoI) models.
#'
#' This function generates binned serosurvey data based on either a time-varying
#' FoI model, an age-varying FoI model, or an age-and-time-varying FoI model.
#' In all cases, it is possible to optionally include seroreversion. This
#' function allows construction of serosurveys with binned age groups, and it
#' generates uncertainty in the distribution of a sample size within an age bin
#' through multinomial sampling.
#'
#' @param model A string specifying the model type which can be either
#' '"age"', '"time"', '"age-time"'.
#' @param foi A dataframe containing the FoI values.
#' For time-varying models the columns should be:
#' \describe{
#'  \item{year}{Calendar years starting at the birth year of the oldest person
#'              and up to the time of the serosurvey}
#'  \item{foi}{Corresponding values of the FoI by year}
#' }
#' For age-varying models the columns should be:.
#' \describe{
#'  \item{age}{Ages starting at 1 and up to the age of the oldest person in the
#'             serosurvey}
#'  \item{foi}{Corresponding values of the FoI by age}
#' }
#' For age-and-time-varying models the columns should be:
#' \describe{
#'  \item{age}{Ages starting at 1 and up to the age of the oldest person in the
#'             serosurvey}
#'  \item{time}{Calendar years starting at the birth year of the oldest person
#'              and up to the time of the serosurvey}
#'  \item{foi}{Corresponding values of FoI by age and year}
#' }
#' @param survey_features A dataframe containing information about the binned
#' age groups and sample sizes for each.
#' It should contain columns:
#' \describe{
#'  \item{age_min}{Left limits of the age groups to be considered in the
#'                 serosurvey}
#'  \item{age_max}{Right limits of the age groups to be considered in the
#'                 serosurvey}
#'  \item{n_sample}{Number of samples by age group}
#' }
#' The resulting age intervals are closed to the left `[` and
#' open to the right `)`.
#' @param seroreversion_rate A non-negative value determining the rate of
#' seroreversion (per year). Default is 0.
#'
#' @return A dataframe with simulated serosurvey data, including age group
#' information, overall sample sizes, the number of seropositive individuals,
#' and other survey features.
#' @examples
#' # time-varying model
#' foi_df <- data.frame(
#'   year = seq(1990, 2009, 1),
#'   foi = rnorm(20, 0.1, 0.01)
#' )
#' survey_features <- data.frame(
#'   age_min = c(1, 3, 15),
#'   age_max = c(2, 14, 20),
#'   n_sample = c(1000, 2000, 1500))
#' serosurvey <- simulate_serosurvey(
#' model = "time",
#' foi = foi_df,
#' survey_features = survey_features)
#'
#' # age-varying model
#' foi_df <- data.frame(
#'   age = seq(1, 20, 1),
#'   foi = rnorm(20, 0.1, 0.01)
#' )
#' survey_features <- data.frame(
#'   age_min = c(1, 3, 15),
#'   age_max = c(2, 14, 20),
#'   n_sample = c(1000, 2000, 1500))
#' serosurvey <- simulate_serosurvey(
#' model = "age",
#' foi = foi_df,
#' survey_features = survey_features)
#'
#' # age-and-time varying model
#' foi_df <- expand.grid(
#'   year = seq(1990, 2009, 1),
#'   age = seq(1, 20, 1)
#' )
#' foi_df$foi <- rnorm(20 * 20, 0.1, 0.01)
#' survey_features <- data.frame(
#'   age_min = c(1, 3, 15),
#'   age_max = c(2, 14, 20),
#'   n_sample = c(1000, 2000, 1500))
#' serosurvey <- simulate_serosurvey(
#' model = "age-time",
#' foi = foi_df,
#' survey_features = survey_features)
#' @export
simulate_serosurvey <- function(
  model,
  foi,
  survey_features,
  seroreversion_rate = 0
) {

  # don't advertise time-age in case people think this is something else
  if (!model %in% c("age", "time", "age-time", "time-age"))
    stop("model must be one of 'age', 'time', or 'age-time'.", call. = FALSE)

  if (model == "time") {
    serosurvey <- simulate_serosurvey_time(
      foi,
      survey_features,
      seroreversion_rate
    )
  } else if (model == "age") {
    serosurvey <- simulate_serosurvey_age(
      foi,
      survey_features,
      seroreversion_rate
    )
  } else if (model == "age-time" || model == "time-age") {
    serosurvey <- simulate_serosurvey_age_time(
      foi,
      survey_features,
      seroreversion_rate
    )
  }

  return(serosurvey)
}

#' Simulate serosurvey data based on general serocatalytic model.
#'
#' This simulation method assumes only that the model system can be written as
#' a piecewise-linear ordinary differential equation system.
#'
#' @inheritParams prob_seroprev_gen_by_age
#' @inheritParams simulate_serosurvey
#'
#' @return A dataframe with simulated serosurvey data, including age group
#' information, overall sample sizes, the number of seropositive individuals,
#' and other survey features.
#'
#' @examples
#' foi_df_time <- data.frame(
#'   year = seq(1946, 2025, 1),
#'   foi = c(rep(0, 40), rep(1, 40))
#' )
#'
#' foi_df_age <- data.frame(
#'   age = 1:80,
#'   foi = 2 * dlnorm(1:80, meanlog = 3.5, sdlog = 0.5)
#' )
#'
#' # generate age and time dependent FoI from multipliers
#' foi_age_time <- expand.grid(
#'   year = foi_df_time$year,
#'   age = foi_df_age$age
#' ) |>
#'   dplyr::left_join(foi_df_age, by = "age") |>
#'   dplyr::rename(foi_age = foi) |>
#'   dplyr::left_join(foi_df_time, by = "year") |>
#'   dplyr::rename(foi_time = foi) |>
#'   dplyr::mutate(foi = foi_age * foi_time) |>
#'   dplyr::select(-c("foi_age", "foi_time"))
#'
#' # create survey features for simulating
#' max_age <- 80
#' n_sample <- 50
#' survey_features <- data.frame(
#'   age_min = seq(1, max_age, 5),
#'   age_max = seq(5, max_age, 5)) |>
#'   dplyr::mutate(n_sample = rep(n_sample, length(age_min))
#'   )
#'
#' # simulate survey from age and time FoI
#' serosurvey <- simulate_serosurvey(
#'   model = "age-time",
#'   foi = foi_age_time,
#'   survey_features = survey_features
#' )
#' @export
simulate_serosurvey_general <- function( #nolint
  construct_A_fun,
  calculate_seroprev_fun, #nolint
  initial_conditions,
  survey_features,
  ...
) {

  # Input validation
  validate_survey_features(survey_features)

  probability_serop_by_age <- prob_seroprev_gen_by_age(
    construct_A_fun,
    calculate_seroprev_fun,
    initial_conditions,
    max(survey_features$age_max),
    ...
  )

  sample_size_by_age_random <- sample_size_by_individual_age(
    survey_features = survey_features
  )

  grouped_df <- get_seroprev_counts_by_bin(
    probability_serop_by_age,
    sample_size_by_age_random,
    survey_features
  )

  grouped_df
}
