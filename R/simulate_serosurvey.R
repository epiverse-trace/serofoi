
#' Create an exposure matrix based on age groups.
#'
#' @param ages A vector containing the ages for which to generate the exposure matrix.
#'             The age groups should be numeric and ordered.
#' @return An exposure matrix indicating exposure status between different age groups.
#'         Rows and columns correspond to ages, and elements indicate exposure status.
create_exposure_matrix <- function(ages) {

  n_ages <- length(ages)
  exposure_matrix <- matrix(0, n_ages, n_ages)
  exposure_matrix <- lower.tri(exposure_matrix, diag = TRUE)
  exposure_matrix[exposure_matrix==TRUE] <- 1

  return(exposure_matrix)
}

#' Computes the probability of being seropositive for age-varying
#' force-of-infection including seroreversion
#'
#' @param ages Integer indicating the ages of the exposed cohorts
#' @param fois Numeric atomic vector corresponding to the age-varying
#' force-of-infection to simulate from
#' @param seroreversion_rate Non-negative seroreversion rate. Default is 0.
#' @return vector of probabilities of being seropositive for age-varying FoI
#' including seroreversion (ordered from youngest to oldest individuals)
probability_exact_time_varying <- function(
    ages,
    fois,
    seroreversion_rate = 0
) {

  exposure_matrix <- create_exposure_matrix(ages)
  probabilities <-
    (fois / (fois + seroreversion_rate)) * (1 - exp(-exposure_matrix %*% (fois + seroreversion_rate)))
  return(probabilities)
}

#' Returns the probability of being seropositive for age-varying
#' force-of-infection including seroreversion
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
    probabilities[i] <-
      (1 / (foi_tmp + seroreversion_rate)) *
      (foi_tmp + exp(-(foi_tmp + seroreversion_rate)) * (probability_previous * (foi_tmp + seroreversion_rate) - foi_tmp))
  }
  return(probabilities)
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

  ages <- seq_along(foi$year)

  probabilities <- probability_exact_time_varying(
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

#' Create a group interval string based on age boundaries.
#'
#' This function generates a group interval string based on the specified age boundaries.
#' It constructs the interval string in the format '[age_min, age_max]' or '(age_min, age_max]',
#' depending on whether it's the first row of a dataframe or not.
#'
#' @param age_min The minimum age of the interval.
#' @param age_max The maximum age of the interval.
#' @param is_first_row Logical indicating whether it's the first row. Default is FALSE.
#'
#' @return A string representing the group interval.
create_group_interval <- function(age_min, age_max, is_first_row=FALSE) {

  first_element <- dplyr::if_else(is_first_row, "[", "(")
  last_element <- "]"
  age_min <- dplyr::if_else(is_first_row, age_min, age_min - 1)
  interval <- paste0(first_element, age_min, ",", age_max, last_element)

  return(interval)
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
    is_first_row <- i == 1
    intervals[i] <- create_group_interval(age_min, age_max, is_first_row)
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
  for (interval_aux in intervals) {
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
  age_bins <- cut(
    ages,
    breaks = c(1, survey_features$age_max),
    include.lowest = TRUE
    )
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
  if (!is.data.frame(foi) || !all(c("year", "foi") %in% names(foi))) {
    stop("foi must be a dataframe with columns 'year' and 'foi'.")
  }
  if (!is.data.frame(survey_features) || !all(c("age_min", "age_max", "sample_size") %in% names(survey_features))) {
    stop("survey_features must be a dataframe with columns 'age_min', 'age_max', and 'sample_size'.")
  }
  if (!is.numeric(seroreversion_rate) || seroreversion_rate < 0) {
    stop("seroreversion_rate must be a non-negative numeric value.")
  }

  probability_serop_by_age <- probability_seropositive_time_model_by_age(
    foi = foi,
    seroreversion_rate = seroreversion_rate
  )

  sample_size_by_age_random <- sample_size_by_individual_age_random(
    survey_features = survey_features
  )

  combined_df <- probability_serop_by_age %>%
    dplyr::left_join(sample_size_by_age_random, by="age") %>%
    dplyr::mutate(
      n_seropositive=rbinom(nrow(probability_serop_by_age),
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
  if (!is.data.frame(foi) || !all(c("age", "foi") %in% names(foi))) {
    stop("foi must be a dataframe with columns 'age' and 'foi'.")
  }
  if (!is.data.frame(survey_features) || !all(c("age_min", "age_max", "sample_size") %in% names(survey_features))) {
    stop("survey_features must be a dataframe with columns 'age_min', 'age_max', and 'sample_size'.")
  }
  if (!is.numeric(seroreversion_rate) || seroreversion_rate < 0) {
    stop("seroreversion_rate must be a non-negative numeric value.")
  }

  probability_serop_by_age <- probability_seropositive_age_model_by_age(
    foi = foi,
    seroreversion_rate = seroreversion_rate
  )

  sample_size_by_age_random <- sample_size_by_individual_age_random(
    survey_features = survey_features
  )

  combined_df <- probability_serop_by_age %>%
    dplyr::left_join(sample_size_by_age_random, by="age") %>%
    dplyr::mutate(
      n_seropositive=rbinom(nrow(probability_serop_by_age),
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
#' # age-and-time-varying model TODO
#' @export
simulate_serosurvey <- function(
    model,
    foi,
    survey_features,
    seroreversion_rate=0
) {

  # don't advertise time-age in case people think this is something else
  if(model %in% c("age", "time", "age-time", "time-age"))
    stop("model must be one of ['age', 'time', 'age-time'].")

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
    serosurvey <- simulate_serosurvey_age_and_time_model( #TODO add age-time
      foi,
      survey_features,
      seroreversion_rate
    )
  }

  return(serosurvey)
}
