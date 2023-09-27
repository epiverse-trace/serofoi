# TODO Complete @param documentation


#' Function that prepares the data from a serological survey for modelling
#'
#' This function adds the necessary additional variables to the given dataset \code{serodata} corresponding to a serological survey.
#' @param serodata A data frame containing the data from a serological survey.
#' This data frame must contain the following columns:
#' \tabular{ll}{
#' \code{survey} \tab survey Label of the current survey \cr \tab \cr
#' \code{total} \tab Number of samples for each age group\cr \tab \cr
#' \code{counts} \tab Number of positive samples for each age group\cr \tab \cr
#' \code{age_min} \tab age_min \cr \tab \cr
#' \code{age_max} \tab age_max \cr \tab \cr
#' \code{tsur} \tab Year in which the survey took place \cr \tab \cr
#' \code{country} \tab The country where the survey took place \cr \tab \cr
#' \code{test} \tab The type of test taken \cr \tab \cr
#' \code{antibody} \tab antibody \cr \tab \cr
#' }
#' @param alpha probability of a type I error. For further details refer to \link[Hmisc]{binconf}.
#' @param add_age_mean_f TBD
#' @return serodata with additional columns necessary for the analysis. These columns are:
#' \tabular{ll}{
#' \code{age_mean_f} \tab Floor value of the average between age_min and age_max \cr \tab \cr
#' \code{sample_size} \tab The size of the sample \cr \tab \cr
#' \code{birth_year} \tab Years in which the subjects were borned 
#' according to the age group marker \code{age_mean_f}\cr \tab \cr
#' according to the age group marker \code{age_mean_f}\cr \tab \cr
#' \code{prev_obs} \tab Observed prevalence \cr \tab \cr
#' \code{prev_obs_lower} \tab Lower limit of the confidence interval for the observed prevalence \cr \tab \cr
#' \code{prev_obs_upper} \tab Upper limit of the confidence interval for the observed prevalence \cr \tab \cr
#' }
#' @examples
#' data(chagas2012)
#' serodata <- prepare_serodata(chagas2012)
#' @export
prepare_serodata <- function(serodata = serodata,
                            alpha = 0.05, 
                            add_age_mean_f = TRUE) {
  if(add_age_mean_f){
    serodata <- serodata %>%
      dplyr::mutate(age_mean_f = floor((age_min + age_max) / 2), sample_size = sum(total)) %>%
      dplyr::mutate(birth_year = .data$tsur - .data$age_mean_f)
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
      prev_obs = PointEst,
      prev_obs_lower = Lower,
      prev_obs_upper = Upper
    ) %>%
    dplyr::arrange(.data$age_mean_f)

  return(serodata)
}


#' Function that prepares a pre-processed serological survey dataset to plot the binomial confidence intervals of the seroprevalence grouped by
#' age group
#'
#' This function prepapares a given pre-processed serological dataset (see \code{\link{prepare_serodata}}) to plot the binomial confidence intervals 
#' of its corresponding seroprevalence grouped by age group. 
#' @param serodata A data frame containing the data from a seroprevalence survey. For more information see the function \link{run_seromodel}.
#' This data frame must contain the following columns:
#' \tabular{ll}{
#' \code{survey} \tab survey Label of the current survey \cr \tab \cr
#' \code{total} \tab Number of samples for each age group\cr \tab \cr
#' \code{counts} \tab Number of positive samples for each age group\cr \tab \cr
#' \code{age_min} \tab age_min \cr \tab \cr
#' \code{age_max} \tab age_max \cr \tab \cr
#' \code{tsur} \tab Year in which the survey took place \cr \tab \cr
#' \code{country} \tab The country where the survey took place \cr \tab \cr
#' \code{test} \tab The type of test taken \cr \tab \cr
#' \code{antibody} \tab antibody \cr \tab \cr
#' \code{age_mean_f} \tab Floor value of the average between age_min and age_max \cr \tab \cr
#' \code{sample_size} \tab The size of the sample \cr \tab \cr
#' \code{birth_year} \tab Years in which the subjects were borned 
#' according to the age group marker \code{age_mean_f}\cr \tab \cr
#' \code{prev_obs} \tab Observed prevalence \cr \tab \cr
#' \code{prev_obs_lower} \tab Lower limit of the confidence interval for the observed prevalence \cr \tab \cr
#' \code{prev_obs_upper} \tab Upper limit of the confidence interval for the observed prevalence \cr \tab \cr
#' }
#' The last six colums can be added to \code{serodata} by means of the function \code{\link{prepare_serodata}}.
#' @return data set with the binomial confidence intervals
#' @examples
#' data(chagas2012)
#' serodata <- prepare_serodata(chagas2012)
#' prepare_bin_data(serodata)
#' @export
prepare_bin_data <- function(serodata) {
  serodata$cut_ages <-
    cut(as.numeric(serodata$age_mean_f),
        seq(1, 101, by = 5),
        include.lowest = TRUE)
  xx <- serodata %>%
    dplyr::group_by(.data$cut_ages) %>%
    dplyr::summarise(bin_size = sum(.data$total),
                     bin_pos = sum(.data$counts))
  labs <-
    read.table(
      text = gsub("[^.0-9]", " ", levels(xx$cut_ages)),
      col.names = c("lower", "upper")
    ) %>%
    dplyr::mutate(lev = levels(xx$cut_ages), mid_age = round((lower + upper) / 2)) %>%
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

#' Function that generates the probabilities of being previously exposed to a pathogen.
#'
#' @param sim_data A dataframe object containing the following columns:
#' \tabular{ll}{
#' \code{tsur} \tab Year of the survey\cr \tab \cr
#' \code{age_mean_f} \tab Age \cr \tab \cr
#' \code{birth_year} \tab Years in which the subjects were borned 
#' according to the age group marker \code{age_mean_f}\cr \tab \cr
#' }
#' @param foi Numeric atomic vector corresponding to the desired Force-of-Infection ordered from past to present
#' @param seed The seed for random number generation.
#' @return A dataframe containing the following columns:
#' \tabular{ll}{
#' \code{age} \tab Exposure ages \cr \tab \cr
#' \code{probability} \tab Probability to obtain a seropositive case for each age according to the provided FoI\cr \tab \cr
#' }
#' @export
get_sim_probability <- function(sim_data, foi) {
  exposure_ages <- get_exposure_ages(sim_data)
  exposure_matrix <- get_exposure_matrix(sim_data)
  probabilities <- purrr::map_dbl(exposure_ages, ~1-exp(-pracma::dot(exposure_matrix[., ], foi)))

  sim_probability <- data.frame(
    age = exposure_ages,
    probability = probabilities
  )
  return(sim_probability)
}

#' Function that generates a sample of counts of seropositive individuals by sampling from a binomial distribution
#'
#' @param sim_data A dataframe object containing the following columns:
#' \tabular{ll}{
#' \code{tsur} \tab Year of the survey\cr \tab \cr
#' \code{age_mean_f} \tab Age \cr \tab \cr
#' \code{birth_year} \tab Years in which the subjects were borned
#' according to the age group marker \code{age_mean_f}\cr \tab \cr
#' }
#' @param foi Numeric atomic vector corresponding to the desired Force-of-Infection ordered from past to present
#' @param sample_size_by_age Sample size for each age group: either a single integer indicating that the sample size is 
#' the same for all ages or a vector of sample sizes the same length as
#' This corresponds to the number of trials \code{size} in \link[stats]{rbinom}.
#' @param seed The seed for random number generation.
#' @return A dataframe containing the following columns:
#' \tabular{ll}{
#' \code{age} \tab Age by the time of the survey \cr \tab \cr
#' \code{n_seropositive} \tab Number of positive cases sampled according to the provided FoI \cr \tab \cr
#' }
#' simulated list of counts following a binomial distribution in accordance with a given
#' force of infection and age class sizes.
#' @examples
#'\dontrun{
#' foi <- rep(0.02, 50)
#' sim_data <- generate_sim_data(foi = foi,
#'                               sample_size_by_age = 5,
#'                               tsur = 2050,
#'                               birth_year_min = 2000,
#'                               survey_label = 'foi_sim')
#' sim_n_seropositive <- get_sim_n_seropositive(sim_data = sim_data,
#'                                              foi = foi)
#' }
#' @export
get_sim_n_seropositive <- function(sim_data, foi, sample_size_by_age, seed = 1234) {
  sim_probability <- get_sim_probability(sim_data = sim_data, foi = foi)

  set.seed(seed = seed)
  n_seropositive <- purrr::map_int(sim_probability$probability, ~rbinom(1, sample_size_by_age, .))

  sim_n_seropositive <- data.frame(
    age = sim_probability$age,
    n_seropositive = n_seropositive
  )
  return(sim_n_seropositive)
}

#' Function that generates a simulated serosurvey according to the specified FoI
#'
#' @param foi Numeric atomic vector corresponding to the desired Force-of-Infection ordered from past to present
#' @param sample_size_by_age Size of each age group specified by either an atomic
#' vector of the same size as \code{foi} or an integer.
#' This corresponds to the number of trials \code{size} in \link[stats]{rbinom}.
#' @return Dataframe object containing the simulated data generated from \code{foi}
#' @examples
#'\dontrun{
#' sample_size_by_age = 5
#' foi <- rep(0.02, 50)
#' sim_data <- generate_sim_data(foi = foi,
#'                               sample_size_by_age = sample_size_by_age,
#'                               tsur = 2050,
#'                               birth_year_min = 2000,
#'                               survey_label = 'sim_constant_foi')
#' }
#' @export
generate_sim_data <- function(foi,
                              sample_size_by_age,
                              tsur,
                              birth_year_min,
                              survey_label,
                              test = "fake",
                              antibody = "IgG",
                              seed = 1234
                              ){
    sim_data <- data.frame(birth_year = c(birth_year_min:(tsur - 1))) %>%
        mutate(tsur = tsur,
            country = 'None',
            test = test,
            antibody = antibody,
            survey = survey_label,
            age_mean_f = tsur - birth_year)
    sim_n_seropositive <- get_sim_n_seropositive(sim_data, foi, sample_size_by_age, seed = seed)
    sim_data <- sim_data %>%
        mutate(counts = sim_n_seropositive$n_seropositive,
               total = sample_size_by_age) %>%
        prepare_serodata(add_age_mean_f = FALSE)

    return(sim_data)
}

#' Method for constructing age-group variable from age column
#'
#' This function was taken from \link[vaccineff]{get_age_group}.
#' This method splits an age interval from age_min to age_max into
#' (age_max-age_min)/step intervals.
#' By default age_min is set 0, however it can be assigned by
#' convenience.
#' If the method finds ages greater or equal than age_max
#' it assigns the string ">{age_max}".
#' To avoid errors it is necessary to set step < age_max.
#' It is also suggested to choose the step such
#' that age_max%%(step+1) = 0.
#' @param age vector containing age information
#' @param  step step used to split the age interval
#' @return age_group factor variable grouping \code{age} by the age intervals 
#' specified by \code{min(age)}, \code{max(age)}.
get_age_group <- function(age, step) {
  age_min <- min(age)
  age_max <- max(age)
  n_steps <- as.integer((age_max - age_min) / step) + 1
  limits_low <- c(as.integer(seq(age_min,
                                 age_max,
                                 length.out = n_steps)))
  limits_hgh <- limits_low + step
  lim_labels <- paste(as.character(limits_low), as.character(limits_hgh),
                      sep = "-")
  lim_labels[length(lim_labels)] <- paste0("+",
                                           limits_low[length(limits_low)])
  lim_breaks <- c(-Inf, limits_low[2:length(limits_low)] - 1, Inf)

  age_group <- cut(age,
                   breaks = lim_breaks,
                   labels = lim_labels,
                   # this is for the intervals to be closed to the left and open to the right
                   right = FALSE) 
  return(age_group)
}

#' Function that groups a simulated serological dataset by age
#'
#' @param sim_data Dataframe with the same structure as the output of \code{\linl{generate_sim_data}}.
#' @param col_age name of the column containing the age information
#' @param step step used to split the age interval
#' @return Dataframe object containing grouped simulated data generated from \code{foi}
group_sim_data <- function(sim_data,
                          col_age = "age_mean_f",
                          step = 5) {
    age <- sim_data[[col_age]]
    sim_data$age_group <-  get_age_group(age = age, step = step)
    sim_data_grouped <- sim_data %>%
      group_by(age_group) %>%
    dplyr::summarise(total = sum(total), counts = sum(counts)) %>%
      mutate(tsur = sim_data$tsur[1],
              country = "None",
              survey = sim_data$survey[1],
              test = sim_data$test[1],
              antibody = sim_data$antibody[1]) %>%
      mutate(age_min = as.integer(sub("\\-.*", "", age_group)),
              age_max = as.integer(sub(".*\\-", "", age_group))) %>%
      prepare_serodata()

    return(sim_data_grouped)
}