#' Computes the probability of being seropositive for age-varying
#' force-of-infection including seroreversion
#'
#' @param ages Integer indicating the ages of the exposed cohorts
#' @param foi Numeric atomic vector corresponding to the age-varying
#' force-of-infection to simulate from
#' @param mu Seroreversion rate
#' @return probability of being seropositive for age-varying FoI
#' including seroreversion
probability_exact_time_varying <- function(
    ages,
    foi,
    mu = 0
) {
  n_ages <- length(ages)
  exposure_matrix <- matrix(0, nrow = n_ages, ncol = n_ages)
  for (k in 1:n_ages) {
    exposure_matrix[k, (n_ages - ages[k] + 1):n_ages] <- 1
  }
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
