# TODO Complete @param documentation


#' Prepare data from a serological survey for modelling
#'
#' This function adds the necessary additional variables to the given dataset
#' `serodata` corresponding to a serological survey.
#' @param serodata A data frame containing the data from a serological survey.
#'  This data frame must contain the following columns:
#' \describe{
#'   \item{`survey`}{survey Label of the current survey}
#'   \item{`total`}{Number of samples for each age group}
#'   \item{`counts`}{Number of positive samples for each age group}
#'   \item{`age_min`}{age_min}
#'   \item{`age_max`}{age_max}
#'   \item{`tsur`}{Year in which the survey took place}
#'   \item{`country`}{The country where the survey took place}
#'   \item{`test`}{The type of test taken}
#'   \item{`antibody`}{antibody}
#' }
#' Alternatively to `age_min` and `age_max`, the dataset could already include
#' the age group marker `age_mean_f`, representing the middle point between
#' `age_min` and `age_max`. If `afe_mean_f` is missing, it will be generated
#' by the function.
#' @param alpha probability of a type I error. For further details refer to
#'   [binconf][Hmisc::binconf].
#' @return serodata with additional columns necessary for the analysis. These
#'   columns are:
#' \describe{
#'   \item{`age_mean_f`}{Floor value of the average between age_min and age_max
#'     for the age groups delimited by `age_min` and `age_max`}
#'   \item{`sample_size`}{The size of the sample}
#'   \item{`birth_year`}{Years in which the subject was born according to the
#'     age group marker `age_mean_f`}
#'   \item{`prev_obs`}{Observed prevalence}
#'   \item{`prev_obs_lower`}{Lower limit of the confidence interval for the
#'     observed prevalence}
#'   \item{`prev_obs_upper`}{Upper limit of the confidence interval for the
#'     observed prevalence}
#' }
#' @examples
#' data(chagas2012)
#' serodata <- prepare_serodata(chagas2012)
#' @export
prepare_serodata <- function(serodata = serodata,
                             alpha = 0.05) {
  checkmate::assert_numeric(alpha, lower = 0, upper = 1)
  serodata <- validate_serodata(serodata)

  if (!any(colnames(serodata) == "age_mean_f")) {
    serodata <- serodata %>%
      dplyr::mutate(
        age_mean_f = floor((.data$age_min + .data$age_max) / 2),
        sample_size = sum(.data$total)
      )
  }

  if (!any(colnames(serodata) == "birth_year")) {
    serodata <- dplyr::mutate(
      serodata,
      birth_year = .data$tsur - .data$age_mean_f
    )
  }

  serodata <- serodata %>%
    cbind(
      Hmisc::binconf(
        serodata$counts,
        serodata$total,
        alpha = alpha,
        method = "exact",
        return.df = TRUE
      )
    ) %>%
    dplyr::rename(
      prev_obs = "PointEst",
      prev_obs_lower = "Lower",
      prev_obs_upper = "Upper"
    ) %>%
    dplyr::arrange(.data$age_mean_f)

  return(serodata)
}


#' Prepare pre-processed serological survey dataset to plot the
#' binomial confidence intervals of the seroprevalence by age group
#'
#' This function prepapares a given pre-processed serological dataset (see
#' [prepare_serodata()]) to plot the binomial confidence intervals of its
#' corresponding seroprevalence grouped by age group.
#' @inheritParams run_seromodel
#' @return data set with the binomial confidence intervals
#' @examples
#' data(chagas2012)
#' serodata <- prepare_serodata(chagas2012)
#' prepare_bin_data(serodata)
#' @keywords internal
#' @noRd
prepare_bin_data <- function(serodata,
                             bin_step = 5,
                             alpha = 0.05) {
  if (!any(colnames(serodata) == "age_mean_f")) {
    serodata <- serodata %>%
      dplyr::mutate(
        age_mean_f = floor((.data$age_min + .data$age_max) / 2),
        sample_size = sum(.data$total)
      )
  }
  serodata$age_group <- get_age_group(
    age = serodata$age_mean_f,
    step = bin_step
  )

  serodata_bin <- serodata %>%
    dplyr::group_by(.data$age_group) %>%
    dplyr::summarise(
      total = sum(.data$total),
      counts = sum(.data$counts)
    ) %>%
    dplyr::mutate(
      survey = unique(serodata$survey),
      tsur = unique(serodata$tsur),
      age_min = as.integer(gsub("[(]|\\,.*", "\\1", .data$age_group)) + 1,
      age_max = as.integer(gsub(".*\\,|[]]", "\\1", .data$age_group)),
      age_mean_f = floor((.data$age_min + .data$age_max) / 2)
    )

    serodata_bin <- cbind(
      serodata_bin,
      Hmisc::binconf(
        serodata_bin$counts,
        serodata_bin$total,
        alpha = alpha,
        method = "exact",
        return.df = TRUE
      )
    ) %>%
    dplyr::rename(
      prev_obs = "PointEst",
      prev_obs_lower = "Lower",
      prev_obs_upper = "Upper"
    ) %>%
    dplyr::arrange(.data$age_mean_f)

  return(serodata_bin)
}

#' Computes the probability of being seropositive for age-varying
#' FoI including seroreversion
#'
#' @param ages Integer indicating the ages of the exposed cohorts
#' @param foi Numeric atomic vector corresponding to the age-varying FoI to
#' simulate from
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

#' Generate probabilities of being previously exposed to a
#' pathogen given a historical force-of-infection.
#'
#' @param sim_data A dataframe object containing the following columns:
#' \describe{
#'   \item{`age`}{Age group markers}
#'   \item{`tsur`}{Year of the survey}
#' }
#' @param foi Numeric atomic vector corresponding to the desired
#' time-varying or age-varying Force-of-Infection to simulate from
#' @param mu Seroreversion rate
#' @param model_type String specifying the type of model to be used.
#' Current valid options are 'time-varying' and 'age-varying'
#' @return A dataframe containing the following columns:
#' \describe{
#'   \item{`age`}{Exposure ages}
#'   \item{`probability`}{Probability to obtain a seropositive case for each age
#'     according to the provided FoI}
#' }
#' @examples
#' n_years <- 50
#' sim_data <- data.frame(
#'   age = seq(1,n_years),
#'   tsur = 2050
#' )
#' foi <- rep(0.02, n_years)
#' sim_probability <- get_sim_probability(sim_data = sim_data, foi=foi)
#' @export
get_sim_probability <- function(
  sim_data,
  foi,
  mu = 0,
  model_type = "time-varying"
  ) {
  # Checks valid model specification
stopifnot(
  "model_type must be either 'time-varying' or 'age-varying'" =
  model_type %in% c("time-varying", "age-varying")
)

  sim_data <- mutate(
    sim_data,
    birth_year = .data$tsur - .data$age
  )
  cohort_ages <- get_cohort_ages(sim_data)
  ages <- rev(cohort_ages$age)

  if (model_type == "time-varying") {
    probabilities <- probability_exact_time_varying(
      ages = ages,
      foi = foi,
      mu = mu
      )
  }
  if (model_type == "age-varying") {
    probabilities <- purrr::map_dbl(
      ages,
      ~probability_exact_age_varying(., foi, mu) # nolint: object_usage_linter
      )
  }

  sim_probability <- data.frame(
    age = ages,
    probability = probabilities
  )
  return(sim_probability)
}

#' Generate sample of counts of seropositive individuals by
#' sampling from a binomial distribution
#'
#' @inheritParams get_sim_probability
#' @param sample_size_by_age Integer indicating the sample size by age group.
#' This corresponds to the number of trials `size` in [rbinom][stats::rbinom].
#' @param seed Seed for random number generation.
#' @return A dataframe containing the following columns:
#' \describe{
#'   \item{`age`}{Age by the time of the survey}
#'   \item{`n_seropositive`}{Number of positive cases sampled according to the
#'     provided FoI}
#' }
#'   simulated list of counts following a binomial distribution in accordance
#'   with a given force of infection and age group sizes.
#' @examples
#' n_years <- 50
#' sim_data <- data.frame(
#'   age = seq(1,n_years),
#'   tsur = 2050
#' )
#' foi <- rep(0.02, n_years)
#' sample_size_by_age <- as.integer(runif(n = n_years, 5, 10))
#' sim_n_seropositive <- get_sim_n_seropositive(
#'   sim_data = sim_data,
#'   foi = foi,
#'   sample_size_by_age = sample_size_by_age
#' )
#' @export
get_sim_n_seropositive <- function(sim_data,
                                   foi,
                                   sample_size_by_age,
                                   mu = 0,
                                   model_type = "time-varying",
                                   seed = 1234) {
  stopifnot(
    "`foi` and `sample_size_by_age` must have the same number of rows" =
      nrow(foi) == nrow(sample_size_by_age)
  )


  sim_probability <- get_sim_probability(
    sim_data = sim_data,
    foi = foi,
    mu = mu,
    model_type = model_type
    )

  set.seed(seed = seed)
  n_seropositive <- purrr::map2_int(
    sample_size_by_age,
    sim_probability$probability,
    ~rbinom(1, size = .x, prob = .y
            )
    )

  sim_n_seropositive <- data.frame(
    age = sim_probability$age,
    n_seropositive = n_seropositive
  )
  return(sim_n_seropositive)
}

#' Generate simulated serosurvey according to the specified FoI
#'
#' @inheritParams get_sim_n_seropositive
#' @param survey_label Label for the resulting simulated serosurvey.
#' @return Dataframe containing the simulated serosurvey.
#' @examples
#' n_years <- 50
#' sim_data <- data.frame(
#'   age = seq(1,n_years),
#'   tsur = 2050
#' )
#' foi <- rep(0.02, n_years)
#' sample_size_by_age <- as.integer(runif(n = n_years, 5, 10))
#' sim_data <- generate_sim_data(
#'   sim_data = sim_data,
#'   foi = foi,
#'   sample_size_by_age = sample_size_by_age,
#'   survey_label = "sim_constant_foi"
#' )
#' @export
generate_sim_data <- function(sim_data,
                              foi,
                              sample_size_by_age,
                              mu = 0,
                              model_type = "time-varying",
                              survey_label = "sim_data",
                              seed = 1234
) {
  sim_n_seropositive <- get_sim_n_seropositive(
    sim_data = sim_data,
    foi = foi,
    sample_size_by_age = sample_size_by_age,
    mu = mu,
    model_type = model_type,
    seed = seed
    )

  # TODO Improve simulation of age_min and age_max
  sim_data <- sim_data %>%
    mutate(
      age_min = .data$age,
      age_max = .data$age,
      counts = sim_n_seropositive$n_seropositive,
      total = sample_size_by_age,
      survey = survey_label
      )

  return(sim_data)
}

#' Construct age-group variable from age column
#'
#' Simplified version of [get_age_group][vaccineff::get_age_group].
#' This function splits an age interval from age_min to age_max by
#' steps of length `step`.
#' `age_min` and `age_max` are calculated from `age`.
#' In cases that `age_max%%(step+1)!=0`, the last age interval is
#' truncated and will have a different length than the others.
#' @param age vector containing age information
#' @param  step step used to split the age interval
#' @return age_group factor variable grouping `age` by the age intervals
#' specified by `min(age)`, `max(age)`.
get_age_group <- function(age, step) {
  age_min <- min(age)
  age_max <- max(age)

  checkmate::assert_int(age_min, lower = 0)
  checkmate::assert_int(age_max, lower = age_min)
  checkmate::assert_int(step, lower = 2, upper = age_max)

  limits_low <- as.integer(
    seq(
      age_min, age_max,
      by = step
      )
    ) - 1

  if ((age_max - age_min) %% step != 0) {
    warn_msg <- "(age_min - age_max) is not an integer multiple of step.
    The last age interval will be truncated to "
    warn_msg <- paste0(
      warn_msg, "(", limits_low[length(limits_low)], ",", age_max, "]"
      )
    warning(warn_msg)
  }

  # prepare breaks
  lim_breaks <- c(limits_low, age_max)

  # define age groups open to the left and closed to the right
  age_group <- cut(
    x = age,
    breaks = lim_breaks
    )

  return(age_group)
}

#' Group simulated serological dataset by age
#'
#' @param sim_data Dataframe with the same structure as the output of
#'   [generate_sim_data()].
#' @param col_age name of the column containing the age information
#' @param step step used to split the age interval
#' @return Dataframe object containing grouped simulated data generated from
#'   `foi`
group_sim_data <- function(sim_data,
                           col_age = "age",
                           step = 5) {
  age <- sim_data[[col_age]]
  sim_data$age_group <- get_age_group(age = age, step = step)
  sim_data_grouped <- sim_data %>%
    group_by(.data$age_group) %>%
    dplyr::summarise(total = sum(.data$total), counts = sum(.data$counts)) %>%
    mutate(
      tsur = sim_data$tsur[1],
      country = "None",
      survey = sim_data$survey[1],
      test = sim_data$test[1],
      antibody = sim_data$antibody[1],
      age_min = as.integer(gsub("[(]|\\,.*", "\\1", .data$age_group)) + 1,
      age_max = as.integer(gsub(".*\\,|[]]", "\\1", .data$age_group))
    ) %>%
    prepare_serodata()

  return(sim_data_grouped)
}
