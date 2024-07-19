
#' Computes the probability of being seropositive when FOIs vary by age
#'
#' @param ages Integer indicating the ages of the exposed cohorts
#' @param fois Numeric atomic vector corresponding to the age-varying
#' force-of-infection to simulate from
#' @param seroreversion_rate Non-negative seroreversion rate. Default is 0.
#' @return vector of probabilities of being seropositive for age-varying FoI
#' including seroreversion (ordered from youngest to oldest individuals)
probability_exact_age_varying <- function(
    ages,
    fois,
    seroreversion_rate = 0
) {

  probabilities <- vector(length = length(ages))
  # solves ODE exactly within pieces
  for (i in seq_along(ages)) {
    foi_tmp <- fois[i]
    if(i == 1)
      probability_previous <- 0
    else
      probability_previous <- probabilities[i - 1]
    if(foi_tmp == 0)
      lambda_over_both <- 0
    else
      lambda_over_both <- (foi_tmp / (foi_tmp + seroreversion_rate))
    probabilities[i] <- lambda_over_both +
      (probability_previous - lambda_over_both) *
      exp(-(foi_tmp + seroreversion_rate))
  }
  return(probabilities)
}

#' Computes the probability of being seropositive when FOIs vary by time
#'
#' @param years Integer indicating the years covering the birth ages of the sample
#' @param fois Numeric atomic vector corresponding to the age-varying
#' force-of-infection to simulate from
#' @param seroreversion_rate Non-negative seroreversion rate. Default is 0.
#' @return vector of probabilities of being seropositive for age-varying FoI
#' including seroreversion (ordered from youngest to oldest individuals)
probability_exact_time_varying <- function(
    years,
    fois,
    seroreversion_rate = 0
) {

  n_years <- length(years)
  ages <- seq(1, n_years, 1)

  probabilities <- vector(length = length(years))
  # solves ODE exactly within pieces
  for (i in seq_along(years)) { # birth cohorts
    probability_previous <- 0
    for(j in 1:(n_years - i + 1)) { # exposure during lifetime
      foi_tmp <- fois[i + j - 1]
      if(foi_tmp == 0)
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

#' Generate probabilities of seropositivity by age based on a time-varying FOI model.
#'
#' This function calculates the probabilities of seropositivity by age based on a time-varying FOI model.
#' It takes into account the FOI and the rate of seroreversion.
#'
#' @param foi A dataframe containing the force of infection (FOI) values for different years.
#'            It should have two columns: 'year' and 'foi'.
#' @param seroreversion_rate A non-negative numeric value representing the rate of seroreversion.
#'
#' @return A dataframe with columns 'age' and 'seropositivity'.
probability_seropositive_time_model_by_age <- function(
    foi,
    seroreversion_rate) {

  years <- foi$year

  probabilities <- probability_exact_time_varying(
    years = years,
    fois = foi$foi,
    seroreversion_rate = seroreversion_rate
  )

  df <- data.frame(
    age = seq_along(years),
    seropositivity = probabilities
  )

  return(df)
}


#' Generate probabilities of seropositivity by age based on an age-varying FOI model.
#'
#' This function calculates the probabilities of seropositivity by age based on an age-varying FOI model.
#' It takes into account the FOI and the rate of seroreversion.
#'
#' @param foi A dataframe containing the force of infection (FOI) values for different ages.
#'            It should have two columns: 'age' and 'foi'.
#' @param seroreversion_rate A non-negative numeric value representing the rate of seroreversion.
#'
#' @return A dataframe with columns 'age' and 'seropositivity'.
probability_seropositive_age_model_by_age <- function(
    foi,
    seroreversion_rate) {

  ages <- seq_along(foi$age)

  probabilities <- probability_exact_age_varying(
    ages = ages,
    fois = foi$foi,
    seroreversion_rate = seroreversion_rate
  )

  df <- data.frame(
    age = ages,
    seropositivity = probabilities
  )

  return(df)
}

#' Generate probabilities of seropositivity by age based on an age-and-time varying FOI model.
#'
#' This function calculates the probabilities of seropositivity by age based on an age-and-time-varying FOI model.
#' It takes into account the FOI and the rate of seroreversion.
#'
#' @param foi A dataframe containing the force of infection (FOI) values for different ages.
#'            It should have three columns: 'year', 'age' and 'foi'.
#' @param seroreversion_rate A non-negative numeric value representing the rate of seroreversion.
#'
#' @return A dataframe with columns 'age' and 'seropositivity'.
probability_seropositive_age_and_time_model_by_age <- function(
    foi,
    seroreversion_rate
    ) {

  foi_matrix <- foi %>%
    tidyr::pivot_wider(
      values_from = foi,
      names_from = c(year)) %>%
    tibble::column_to_rownames("age") %>%
    as.matrix()

  years <- unique(foi$year)
  n_years <- length(years)
  ages <- unique(foi$age)
  n_ages <- length(ages)

  probabilities <- vector(length = length(n_ages))
  # solves ODE exactly within pieces
  for (i in seq_along(years)) { # birth cohorts
    probability_previous <- 0
    foi_matrix_subset <- foi_matrix[1:(n_ages - i + 1), i:n_ages] %>%
      as.matrix() # only to handle single element matrix case
    foi_diag <- diag(foi_matrix_subset)
    for(j in 1:(n_years - i + 1)) { # exposure during lifetime
      foi_tmp <- foi_diag[j]
      if(foi_tmp == 0)
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

  df <- data.frame(
    age = ages,
    seropositivity = probabilities_oldest_age_last
  )

  return(df)
}

#' Generate probabilities of seropositivity by age based user's choice of model.
#'
#' This function generates seropositivity probabilities based on either a
#' time-varying FOI model, an age-varying FOI model, or an age-and-time-varying
#' FOI model. In all cases, it is possible to optionally include seroreversion.
#'
#' @param model A string specifying the model type which can be one of ['age', 'time', 'age-time'].
#' @param foi A dataframe containing the force of infection (FOI) values.
#'            For time-varying models the columns should be ['year', 'foi'].
#'            For age-varying models the columns should be ['age', 'foi'].
#'            For age-and-time-varying models the columns should be ['age', 'time', 'foi'].
#' @param seroreversion_rate A non-negative value determining the rate of seroreversion (per year).
#'                           Default is 0.
#'
#' @return A dataframe with columns 'age' and 'seropositivity'.
#' @export
probability_seropositive_by_age <- function(
    model,
    foi,
    seroreversion_rate = 0) {

  if(model == "time") {
    probability_function <- probability_seropositive_time_model_by_age
  } else if(model == "age") {
    probability_function <- probability_seropositive_age_model_by_age
  } else if(model == "age-time" || model == "time-age") {
    probability_function <- probability_seropositive_age_and_time_model_by_age
  }

  df <- probability_function(
    foi = foi,
    seroreversion_rate = seroreversion_rate
  )

  return(df)
}


#' Add bins based on age intervals.
#'
#' It generates a new column 'group' in the survey_features dataframe, representing
#' the group interval for each row based on the age_min and age_max columns.
#'
#' @param survey_features A dataframe containing age_min and age_max columns representing
#'                        the minimum and maximum age boundaries for each group.
#'
#' @return A dataframe with an additional 'group' column representing the group interval
#'         for each row based on the age_min and age_max columns.
add_age_bins <- function(survey_features) {
  intervals <- vector(length = nrow(survey_features))
  for(i in seq_along(intervals)) {
    age_min <- survey_features$age_min[i]
    age_max <- survey_features$age_max[i]
    intervals[i] <- paste0("[", age_min, ",", age_max, "]")
  }
  survey_features <- survey_features %>%
    mutate(group = intervals)
  return(survey_features)
}

#' Create a survey dataframe with per individual age information.
#'
#'
#' @param survey_features A dataframe containing information about age groups and sample sizes.
#' @param age_df A dataframe containing 'age' and 'group'.
#'
#' @return A dataframe with overall sample sizes calculated by joining survey_features and age_df.
#'         This dataframe has columns including 'age' and 'overall_sample_size'.
survey_by_individual_age <- function(survey_features, age_df) {
  age_df %>%
    left_join(survey_features, by = "group") %>%
    rename(overall_sample_size = sample_size)
}

#' Generate random sample sizes using multinomial sampling.
#'
#' This function generates random sample sizes for each age group using multinomial sampling.
#' It takes the total sample size and the number of age groups as input and returns a vector
#' of sample sizes for each age group.
#'
#' @param sample_size The total sample size to be distributed among age groups.
#' @param n_ages The number of age groups.
#'
#' @return A vector containing random sample sizes for each age group.
multinomial_sampling_group <- function(sample_size, n_ages) {
  prob_value <- 1 / n_ages
  probs <- rep(prob_value, n_ages)
  sample_size_by_age <- as.vector(
    rmultinom(1, sample_size, prob=probs)
  )
  return(sample_size_by_age)
}

#' Generate random sample sizes for each age group.
#'
#' This function generates random sample sizes for each age group based on the overall sample size
#' and the distribution of individuals across age groups. It uses multinomial sampling to allocate
#' the total sample size to each age group proportionally.
#'
#' @param survey_df_long A dataframe with columns 'age', 'group' and 'overall_sample_size'.
#'
#' @return A dataframe with random sample sizes generated for each age based on the overall
#'         sample size.
generate_random_sample_sizes <- function(survey_df_long) {

  df_new <- NULL
  intervals <- unique(survey_df_long$group)
  for (interval_aux in na.omit(intervals)) {
    df_tmp <- survey_df_long %>%
      filter(group == interval_aux)
    sample_size <- df_tmp$overall_sample_size[1]
    sample_size_by_age <- multinomial_sampling_group(sample_size, nrow(df_tmp))
    df_tmp <- df_tmp %>%
      mutate(sample_size = sample_size_by_age)

    if (is.null(df_new)) {
      df_new <- df_tmp
    } else {
      df_new <- bind_rows(df_new, df_tmp)
    }
  }
  return(df_new)
}

#' Generate random sample sizes for each individual age based on survey features.
#'
#' This function generates random sample sizes for each individual age based on the provided
#' survey features. It first creates age bins, assigns each individual in the survey features
#' to an age bin, calculates the overall sample size by group, and then generates random sample
#' sizes for each age group. Finally, it returns a dataframe with the random sample sizes for
#' each individual age.
#'
#' @param survey_features A dataframe containing information about individuals' age ranges and
#'                        sample sizes.
#'
#' @return A dataframe with random sample sizes generated for each individual age based on the
#'         provided survey features.
sample_size_by_individual_age_random <- function(survey_features) {

  ages <- seq(1, max(survey_features$age_max), 1)

  age_bins <- rep(NA, length(ages))
  for (i in seq_along(ages)) {
    # Find the index of the row in survey_features where age falls within the range
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

  survey_features_by_individual_age <- survey_by_individual_age(
    survey_features,
    age_df)

  df_new <- generate_random_sample_sizes(survey_features_by_individual_age)

  return(df_new)
}

check_age_constraints <- function(df) {
  for (i in 1:nrow(df)) {
    for (j in 1:nrow(df)) {
      if (i != j && df$age_max[i] == df$age_min[j]) {
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

validate_survey <- function(survey_features) {

  if (!is.data.frame(survey_features) || !all(c("age_min", "age_max", "sample_size") %in% names(survey_features))) {
    stop("survey_features must be a dataframe with columns 'age_min', 'age_max', and 'sample_size'.")
  }

  # check that the age_max of a bin does not coincide with
  # the age min of a different bin
  is_age_ok <- check_age_constraints(survey_features)
  if(!is_age_ok)
    stop("Age bins in a survey are inclusive of both bounds, so the age_max of
         one bin cannot equal the age_min of another.")
}

validate_foi_df <- function(foi_df, cnames_additional) {
  if (!is.data.frame(foi_df) || !all(cnames_additional %in% names(foi_df)) || ncol(foi_df) != (1 + length(cnames_additional))) {
    if(length(cnames_additional) == 1)
      message_end <- paste0(" and ", cnames_additional, ".")
    else
      message_end <- paste0(", ", paste(cnames_additional, collapse=" and "), ".")
    message_beginning <- "foi must be a dataframe with columns foi"
    stop(paste0(message_beginning, message_end))
  }
}

validate_seroreversion_rate <- function(seroreversion_rate) {
  if (!is.numeric(seroreversion_rate) || seroreversion_rate < 0) {
    stop("seroreversion_rate must be a non-negative numeric value.")
  }
}

validate_survey_and_foi_consistency <- function(
    survey_features,
    foi_df
) {
  max_age_foi_df <- nrow(foi_df)
  if(max_age_foi_df > max(survey_features$age_max))
    stop("maximum age implicit in foi_df should not exceed max age in survey_features.")
}

validate_survey_and_foi_consistency_age_time <- function(
    survey_features,
    foi_df
) {
  max_age_foi_df <- max(foi_df$year) - min(foi_df$year) + 1
  if(max_age_foi_df > max(survey_features$age_max))
    stop("maximum age implicit in foi_df should not exceed max age in survey_features.")
}

generate_seropositive_counts_by_age_bin <- function(
    probability_seropositive_by_age,
    sample_size_by_age_random,
    survey_features
    ) {

  combined_df <- probability_seropositive_by_age %>%
    dplyr::left_join(sample_size_by_age_random, by="age") %>%
    dplyr::mutate(
      n_seropositive=rbinom(nrow(probability_seropositive_by_age),
                            sample_size,
                            seropositivity))

  grouped_df <- combined_df %>%
    dplyr::group_by(age_min, age_max) %>%
    dplyr::summarise(
      sample_size=sum(sample_size),
      n_seropositive=sum(n_seropositive),
      .groups = "drop"
    ) %>%
    left_join(
      survey_features,
      by = c("age_min", "age_max", "sample_size")
    )

  return(grouped_df)
}

#' Simulate serosurvey data based on a time-varying FOI model.
#'
#' This function generates binned serosurvey data based on a time-varying FOI model,
#' optionally including seroreversion. This function allows construction of serosurveys
#' with binned age groups, and it generates uncertainty in the distribution of a sample size
#' within an age bin through multinomial sampling.
#'
#' @param foi A dataframe containing the force of infection (FOI) values for different years.
#'            It should have two columns: 'year' and 'foi'.
#' @param survey_features A dataframe containing information about the binned age groups and sample
#'                        sizes for each. It should contain columns: ['age_min', 'age_max', 'sample_size'].
#' @param seroreversion_rate A non-negative value determining the rate of seroreversion (per year).
#'                           Default is 0.
#'
#' @return A dataframe with simulated serosurvey data, including age group information, overall
#'         sample sizes, the number of seropositive individuals, and other survey features.
#' @examples
#' # specify FOIs for each year
#' foi_df <- data.frame(
#'   year = seq(1990, 2009, 1)
#' ) %>%
#' mutate(foi = rnorm(length(year), 0.1, 0.01))
#' survey_features <- data.frame(
#'   age_min = c(1, 3, 15),
#'   age_max = c(2, 14, 20),
#'   sample_size = c(1000, 2000, 1500))
#' serosurvey <- simulate_serosurvey_time_model(
#' foi_df, survey_features)
#' @export
simulate_serosurvey_time_model <- function(
    foi,
    survey_features,
    seroreversion_rate=0
) {

  # Input validation
  validate_foi_df(foi, c("year"))
  validate_survey(survey_features)
  validate_seroreversion_rate(seroreversion_rate)
  validate_survey_and_foi_consistency(
    survey_features,
    foi
    )

  probability_serop_by_age <- probability_seropositive_time_model_by_age(
    foi = foi,
    seroreversion_rate = seroreversion_rate
  )

  sample_size_by_age_random <- sample_size_by_individual_age_random(
    survey_features = survey_features
  )

  grouped_df <- generate_seropositive_counts_by_age_bin(
    probability_serop_by_age,
    sample_size_by_age_random,
    survey_features
  )

  return(grouped_df)
}


#' Simulate serosurvey data based on an age-varying FOI model.
#'
#' This function generates binned serosurvey data based on an age-varying FOI model,
#' optionally including seroreversion. This function allows construction of serosurveys
#' with binned age groups, and it generates uncertainty in the distribution of a sample size
#' within an age bin through multinomial sampling.
#'
#' @param foi A dataframe containing the force of infection (FOI) values for different ages.
#'            It should have two columns: 'age' and 'foi'.
#' @param survey_features A dataframe containing information about the binned age groups and sample
#'                        sizes for each. It should contain columns: ['age_min', 'age_max', 'sample_size'].
#' @param seroreversion_rate A non-negative value determining the rate of seroreversion (per year).
#'                           Default is 0.
#'
#' @return A dataframe with simulated serosurvey data, including age group information, overall
#'         sample sizes, the number of seropositive individuals, and other survey features.
#' @examples
#' # specify FOIs for each year
#' foi_df <- data.frame(
#'   age = seq(1, 20, 1)
#' ) %>%
#' mutate(foi = rnorm(length(year), 0.1, 0.01))
#' survey_features <- data.frame(
#'   age_min = c(1, 3, 15),
#'   age_max = c(2, 14, 20),
#'   sample_size = c(1000, 2000, 1500))
#' serosurvey <- simulate_serosurvey_age_model(
#' foi_df, survey_features)
#' @export
simulate_serosurvey_age_model <- function(
    foi,
    survey_features,
    seroreversion_rate=0
) {

  # Input validation
  validate_foi_df(foi, c("age"))
  validate_survey(survey_features)
  validate_seroreversion_rate(seroreversion_rate)
  validate_survey_and_foi_consistency(
    survey_features,
    foi
  )

  probability_serop_by_age <- probability_seropositive_age_model_by_age(
    foi = foi,
    seroreversion_rate = seroreversion_rate
  )

  sample_size_by_age_random <- sample_size_by_individual_age_random(
    survey_features = survey_features
  )

  grouped_df <- generate_seropositive_counts_by_age_bin(
    probability_serop_by_age,
    sample_size_by_age_random,
    survey_features
  )

  return(grouped_df)
}

#' Simulate serosurvey data based on an age-and-time-varying FOI model.
#'
#' This function generates binned serosurvey data based on an age-and-time-varying FOI model,
#' optionally including seroreversion. This function allows construction of serosurveys
#' with binned age groups, and it generates uncertainty in the distribution of a sample size
#' within an age bin through multinomial sampling.
#'
#' @param foi A dataframe containing the force of infection (FOI) values for different ages.
#'            It should have two columns: 'year', 'age' and 'foi'.
#' @param survey_features A dataframe containing information about the binned age groups and sample
#'                        sizes for each. It should contain columns: ['age_min', 'age_max', 'sample_size'].
#' @param seroreversion_rate A non-negative value determining the rate of seroreversion (per year).
#'                           Default is 0.
#'
#' @return A dataframe with simulated serosurvey data, including age group information, overall
#'         sample sizes, the number of seropositive individuals, and other survey features.
#' @examples
#' # specify FOIs for each year
#' foi_df <- data.frame(
#'   year = seq(1990, 2009, 1),
#'   age = seq(1, 20, 1)
#' ) %>%
#' mutate(foi = rnorm(20 * 20, 0.1, 0.01))
#' survey_features <- data.frame(
#'   age_min = c(1, 3, 15),
#'   age_max = c(2, 14, 20),
#'   sample_size = c(1000, 2000, 1500))
#' serosurvey <- simulate_serosurvey_age_and_time_model(
#' foi_df, survey_features)
#' @export
simulate_serosurvey_age_and_time_model <- function(
    foi,
    survey_features,
    seroreversion_rate=0
) {

  # Input validation
  validate_foi_df(foi, c("age", "year"))
  validate_survey(survey_features)
  validate_seroreversion_rate(seroreversion_rate)
  validate_survey_and_foi_consistency_age_time(
    survey_features,
    foi
  )

  probability_serop_by_age <- probability_seropositive_age_and_time_model_by_age(
    foi = foi,
    seroreversion_rate = seroreversion_rate
  )

  sample_size_by_age_random <- sample_size_by_individual_age_random(
    survey_features = survey_features
  )

  grouped_df <- generate_seropositive_counts_by_age_bin(
    probability_serop_by_age,
    sample_size_by_age_random,
    survey_features
  )

  return(grouped_df)
}


#' Simulate serosurvey data based on various FOI models.
#'
#' This function generates binned serosurvey data based on either a time-varying FOI model,
#' an age-varying FOI model, or an age-and-time-varying FOI model. In all cases, it is possible
#' to optionally include seroreversion. This function allows construction of serosurveys
#' with binned age groups, and it generates uncertainty in the distribution of a sample size
#' within an age bin through multinomial sampling.
#'
#' @param model A string specifying the model type which can be one of ['age', 'time', 'age-time'].
#' @param foi A dataframe containing the force of infection (FOI) values.
#'            For time-varying models the columns should be ['year', 'foi'].
#'            For age-varying models the columns should be ['age', 'foi'].
#'            For age-and-time-varying models the columns should be ['age', 'time', 'foi'].
#' @param survey_features A dataframe containing information about the binned age groups and sample
#'                        sizes for each. It should contain columns: ['age_min', 'age_max', 'sample_size'].
#' @param seroreversion_rate A non-negative value determining the rate of seroreversion (per year).
#'                           Default is 0.
#'
#' @return A dataframe with simulated serosurvey data, including age group information, overall
#'         sample sizes, the number of seropositive individuals, and other survey features.
#' @examples
#' # time-varying model
#' foi_df <- data.frame(
#'   year = seq(1990, 2009, 1)
#' ) %>%
#' mutate(foi = rnorm(length(year), 0.1, 0.01))
#' survey_features <- data.frame(
#'   age_min = c(1, 3, 15),
#'   age_max = c(2, 14, 20),
#'   sample_size = c(1000, 2000, 1500))
#' serosurvey <- simulate_serosurvey(
#' model = "time",
#' foi = foi_df,
#' survey_features = survey_features)
#'
#' # age-varying model
#' foi_df <- data.frame(
#'   age = seq(1, 20, 1)
#' ) %>%
#' mutate(foi = rnorm(length(year), 0.1, 0.01))
#' survey_features <- data.frame(
#'   age_min = c(1, 3, 15),
#'   age_max = c(2, 14, 20),
#'   sample_size = c(1000, 2000, 1500))
#' serosurvey <- simulate_serosurvey(
#' model = "age",
#' foi = foi_df,
#' survey_features = survey_features)
#'
#' # age-and-time varying model
#' foi_df <- expand.grid(
#'   year = seq(1990, 2009, 1),
#'   age = seq(1, 20, 1)
#' ) %>%
#' mutate(foi = rnorm(20 * 20, 0.1, 0.01))
#' survey_features <- data.frame(
#'   age_min = c(1, 3, 15),
#'   age_max = c(2, 14, 20),
#'   sample_size = c(1000, 2000, 1500))
#' serosurvey <- simulate_serosurvey(
#' model = "age",
#' foi = foi_df,
#' survey_features = survey_features)
#' @export
simulate_serosurvey <- function(
    model,
    foi,
    survey_features,
    seroreversion_rate=0
) {

  # don't advertise time-age in case people think this is something else
  if(!model %in% c("age", "time", "age-time", "time-age"))
    stop("model must be one of 'age', 'time', or 'age-time'.")

  if(model == "time") {
    serosurvey <- simulate_serosurvey_time_model(
      foi,
      survey_features,
      seroreversion_rate
    )
  } else if(model == "age") {
    serosurvey <- simulate_serosurvey_age_model(
      foi,
      survey_features,
      seroreversion_rate
    )
  } else if(model == "age-time" || model == "time-age") {
    serosurvey <- simulate_serosurvey_age_and_time_model(
      foi,
      survey_features,
      seroreversion_rate
    )
  }

  return(serosurvey)
}
