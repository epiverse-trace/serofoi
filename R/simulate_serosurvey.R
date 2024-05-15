create_exposure_matrix <- function(ages) {

  n_ages <- length(ages)
  exposure_matrix <- matrix(NA, n_ages, n_ages)
  exposure_matrix <- lower.tri(exposure_matrix, diag = TRUE)
  exposure_matrix[exposure_matrix==TRUE] <- 1

  return(exposure_matrix)
}

#' Computes the probability of being seropositive for age-varying
#' force-of-infection including seroreversion
#'
#' @param ages Integer indicating the ages of the exposed cohorts
#' @param foi Numeric atomic vector corresponding to the age-varying
#' force-of-infection to simulate from
#' @param mu Seroreversion rate
#' @return vector of probabilities of being seropositive for age-varying FoI
#' including seroreversion (ordered from youngest to oldest individuals)
probability_exact_time_varying <- function(
    ages,
    foi,
    mu = 0
) {

  exposure_matrix <- create_exposure_matrix(ages)
  probabilities <-
    (foi / (foi + mu)) * (1 - exp(-drop(exposure_matrix %*% (foi + mu))))
  return(probabilities)
}

#' Returns the probability of being seropositive for age-varying
#' force-of-infection including seroreversion
#'
#' @param age Integer corresponding to the age of the exposed cohort
#' @param foi Numeric atomic vector corresponding to the age-varying FoI to
#' simulate from
#' @param mu Seroreversion rate
#' @return probability of being seropositive for age-varying FoI
#' including seroreversion
probability_exact_age_varying <- function(
    age,
    foi,
    mu = 0
) {
  probability <- 0
  # solves ODE exactly within pieces
  for (i in 1:age) {
    probability <-
      (1 / (foi[i] + mu)) *
      (foi[i] + exp(-(foi[i] + mu)) * (probability * (foi[i] + mu) - foi[i]))
  }
  return(probability)
}

probability_seropositive_time_model_by_age <- function(
    foi,
    seroreversion) {

  ages <- seq_along(foi$year)

  probabilities <- probability_exact_time_varying(
    age = ages,
    foi = foi$foi,
    mu = seroreversion
  )

  df <- data.frame(
    age = ages,
    seropositivity = probabilities
  )

  return(df)
}

create_group_interval <- function(age_min, age_max, is_first_row=FALSE) {

  first_element <- dplyr::if_else(is_first_row, "[", "(")
  last_element <- "]"
  age_min <- dplyr::if_else(is_first_row, age_min, age_min - 1)
  interval <- paste0(first_element, age_min, ",", age_max, last_element)

  return(interval)
}

create_features_with_bins <- function(survey_features) {
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

overall_sample_size_by_group <- function(survey_features, age_df) {
  age_df %>%
    left_join(survey_features, by = "group") %>%
    rename(overall_sample_size = sample_size)
}

multinomial_sampling_group <- function(sample_size, n_ages) {
  prob_value <- 1 / n_ages
  probs <- rep(prob_value, n_ages)
  sample_size_by_age <- as.vector(
    rmultinom(1, sample_size, prob=probs)
  )
  return(sample_size_by_age)
}

generate_random_sample_sizes <- function(survey_df) {
  df_new <- NULL
  intervals <- unique(survey_df$group)
  for (interval_aux in intervals) {
    df_tmp <- survey_df %>%
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

sample_size_by_individual_age_random <- function(survey_features) {

  age_bins <- seq(1, max(survey_features$age_max), 1)
  age_df <- data.frame(
    age = age_bins,
    group = cut(age_bins,
                breaks = c(1, survey_features$age_max),
                include.lowest = TRUE)
    )

  survey_features <- create_features_with_bins(survey_features)

  survey_features <- overall_sample_size_by_group(
    survey_features,
    age_df)

  df_new <- generate_random_sample_sizes(survey_features)

  return(df_new)
}

simulate_serosurvey_time_model <- function(
    foi,
    survey_features,
    seroreversion=FALSE
) {

  probability_serop_by_age <- probability_seropositive_time_model_by_age(
    foi = foi,
    seroreversion = seroreversion
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
      n_seropositive=sum(n_seropositive)
    ) %>%
    left_join(
      survey_features,
      by = c("age_min", "age_max", "sample_size")
    )

  return(grouped_df)
}
