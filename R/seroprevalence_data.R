# TODO Complete @param documentation


#' Function that prepares the data from a serological survey for modelling
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
  validate_serodata(serodata)


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


#' Function that prepares a pre-processed serological survey dataset to plot the
#' binomial confidence intervals of the seroprevalence grouped by age group
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
prepare_bin_data <- function(serodata) {
  if (!any(colnames(serodata) == "age_mean_f")) {
    serodata <- serodata %>%
      dplyr::mutate(
        age_mean_f = floor((.data$age_min + .data$age_max) / 2),
        sample_size = sum(.data$total)
      )
  }
  serodata$cut_ages <-
    cut(as.numeric(serodata$age_mean_f),
      seq(1, 101, by = 5),
      include.lowest = TRUE
    )
  xx <- serodata %>%
    dplyr::group_by(.data$cut_ages) %>%
    dplyr::summarise(
      bin_size = sum(.data$total),
      bin_pos = sum(.data$counts)
    )
  labs <-
    read.table(
      text = gsub("[^.0-9]", " ", levels(xx$cut_ages)),
      col.names = c("lower", "upper")
    ) %>%
    dplyr::mutate(
      lev = levels(xx$cut_ages),
      mid_age = round((.data$lower + .data$upper) / 2)
    ) %>%
    dplyr::select(.data$mid_age, .data$lev)
  xx$mid_age <- labs$mid_age[labs$lev %in% xx$cut_ages]
  conf <-
    data.frame(Hmisc::binconf(xx$bin_pos, xx$bin_size, method = "exact"))
  xx <- cbind(xx, conf) %>% dplyr::rename(
    age = .data$mid_age,
    p_obs_bin = .data$PointEst,
    p_obs_bin_l = .data$Lower,
    p_obs_bin_u = .data$Upper
  )
  return(xx)
}

#' Function that generates the probabilities of being previously exposed to a
#' pathogen given a historical Force-of-Infection.
#'
#' @param sim_data A dataframe object containing the following columns:
#' \describe{
#'   \item{`age_mean_f`}{Age group markers}
#'   \item{`tsur`}{Year of the survey}
#' }
#' @param foi Numeric atomic vector corresponding to the desired
#'   Force-of-Infection ordered from past to present
#' @return A dataframe containing the following columns:
#' \describe{
#'   \item{`age`}{Exposure ages}
#'   \item{`probability`}{Probability to obtain a seropositive case for each age
#'     according to the provided FoI}
#' }
#' @examples
#' n_years <- 50
#' sim_data <- data.frame(
#'   age_mean_f = seq(1,n_years),
#'   tsur = 2050
#' )
#' foi <- rep(0.02, n_years)
#' sim_probability <- get_sim_probability(sim_data = sim_data, foi=foi)
#' @export
get_sim_probability <- function(sim_data, foi) {
  sim_data <- sim_data %>%
    mutate(birth_year = .data$tsur - .data$age_mean_f)
  cohort_ages <- get_cohort_ages(sim_data)
  exposure_ages <- rev(cohort_ages$age)
  exposure_matrix <- get_exposure_matrix(sim_data) # nolint: object_usage_linter
  probabilities <- purrr::map_dbl(
    exposure_ages,
    ~ 1 - exp(-pracma::dot(exposure_matrix[., ], foi))
  )

  sim_probability <- data.frame(
    age = exposure_ages,
    probability = probabilities
  )
  return(sim_probability)
}

#' Function that generates a sample of counts of seropositive individuals by
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
#'   age_mean_f = seq(1,n_years),
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
                                   seed = 1234) {
  sim_probability <- get_sim_probability(sim_data = sim_data, foi = foi)

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

#' Function that generates a simulated serosurvey according to the specified FoI
#'
#' @inheritParams get_sim_n_seropositive
#' @param survey_label Label for the resulting simulated serosurvey.
#' @return Dataframe containing the simulated serosurvey.
#' @examples
#' n_years <- 50
#' sim_data <- data.frame(
#'   age_mean_f = seq(1,n_years),
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
                              survey_label,
                              seed = 1234
) {
  sim_n_seropositive <- get_sim_n_seropositive(
    sim_data,
    foi,
    sample_size_by_age,
    seed = seed
    )

  # TODO Improve simulation of age_min and age_max
  sim_data <- sim_data %>%
    mutate(
      age_min = age_mean_f,
      age_max = age_mean_f,
      counts = sim_n_seropositive$n_seropositive,
      total = sample_size_by_age,
      survey = survey_label
      )

  return(sim_data)
}

#' Method for constructing age-group variable from age column
#'
#' This function was taken from [get_age_group][vaccineff::get_age_group].
#' This method splits an age interval from age_min to age_max into
#' `(age_max-age_min)/step` intervals.
#' By default age_min is set 0, however it can be assigned by
#' convenience.
#' If the method finds ages greater or equal than age_max
#' it assigns the string `">{age_max}"`.
#' To avoid errors it is necessary to set `step<age_max`.
#' It is also suggested to choose the step such
#' that `age_max%%(step+1)=0`.
#' @param age vector containing age information
#' @param  step step used to split the age interval
#' @return age_group factor variable grouping `age` by the age intervals
#' specified by `min(age)`, `max(age)`.
get_age_group <- function(age, step) {
  age_min <- min(age)
  age_max <- max(age)
  n_steps <- as.integer((age_max - age_min) / step) + 1
  limits_low <- c(as.integer(seq(age_min,
    age_max,
    length.out = n_steps
  )))
  limits_hgh <- limits_low + step
  lim_labels <- paste(as.character(limits_low), as.character(limits_hgh),
    sep = "-"
  )
  lim_labels[length(lim_labels)] <- paste0(
    "+",
    limits_low[length(limits_low)]
  )
  lim_breaks <- c(-Inf, limits_low[2:length(limits_low)] - 1, Inf)

  age_group <- cut(age,
    breaks = lim_breaks,
    labels = lim_labels,
    # this is for the intervals to be closed to the left and open to the right
    right = FALSE
  )
  return(age_group)
}

#' Function that groups a simulated serological dataset by age
#'
#' @param sim_data Dataframe with the same structure as the output of
#'   [generate_sim_data()].
#' @param col_age name of the column containing the age information
#' @param step step used to split the age interval
#' @return Dataframe object containing grouped simulated data generated from
#'   `foi`
group_sim_data <- function(sim_data,
                           col_age = "age_mean_f",
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
      age_min = as.integer(sub("\\-.*", "", .data$age_group)),
      age_max = as.integer(sub(".*\\-", "", .data$age_group))
    ) %>%
    prepare_serodata()

  return(sim_data_grouped)
}
