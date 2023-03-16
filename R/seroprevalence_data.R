#' Prepare data
#'
#' Function that prepares the data for modelling
#' @param seroprev_data A data frame containing the data from a seroprevalence survey.
#' This data frame must contain the following columns:
#' \tabular{ll}{
#' \code{survey} \tab survey Label of the current survey \cr \tab \cr
#' \code{total} \tab Number of samples for each age group\cr \tab \cr
#' \code{counts} \tab Number of positive samples for each age group\cr \tab \cr
#' \code{age_min} \tab age_min \cr \tab \cr
#' \code{age_max} \tab age_max \cr \tab \cr
#' \code{year_init} \tab year_init \cr \tab \cr
#' \code{year_end} \tab year_end \cr \tab \cr
#' \code{tsur} \tab Year in which the survey took place \cr \tab \cr
#' \code{country} \tab The country where the survey took place \cr \tab \cr
#' \code{test} \tab The type of test taken \cr \tab \cr
#' \code{antibody} \tab antibody \cr \tab \cr
#' }
#' @param alpha probability of a type I error. For further details refer to \link{Hmisc::binconf}.
#' @return seroprev_data with additional columns necessary for the analysis. These columns are:
#' \tabular{ll}{
#' \code{age_mean_f} \tab Floor value of the average between age_min and age_max \cr \tab \cr
#' \code{sample_size} \tab The size of the sample \cr \tab \cr
#' \code{birth_year} \tab The year in which the individuals of each age group were bornt \cr \tab \cr
#' \code{prev_obs} \tab Observed prevalence \cr \tab \cr
#' \code{prev_obs_lower} \tab Lower limit of the confidence interval for the observed prevalence \cr \tab \cr
#' \code{prev_obs_upper} \tab Upper limit of the confidence interval for the observed prevalence \cr \tab \cr
#' }
#' @examples
#'\dontrun{
#' data_test <- readRDS("data/data.RDS")
#' data_test <- prepare_seroprev_data(seroprev_data, alpha)
#' }
#' @export
prepare_seroprev_data <- function(seroprev_data = serodata,
                                  alpha = 0.05, 
                                  add_age_mean_f = TRUE) {
  if(add_age_mean_f){
    seroprev_data <- seroprev_data %>%
      dplyr::mutate(age_mean_f = floor((age_min + age_max) / 2), sample_size = sum(total)) %>%
      dplyr::mutate(birth_year = .data$tsur - .data$age_mean_f)
  }
  seroprev_data <- seroprev_data %>%
    cbind(
      Hmisc::binconf(
        seroprev_data$counts,
        seroprev_data$total,
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

  return(seroprev_data)
}


#' Prepare data to plot binomial confidence intervals
#'
#' Function that prepares the data for modelling
#' @param seroprev_data A data frame containing the data from a seroprevalence survey. For more information see the function \link{run_seroprev_model}.
#' This data frame must contain the following columns:
#' \tabular{ll}{
#' \code{survey} \tab survey Label of the current survey \cr \tab \cr
#' \code{total} \tab Number of samples for each age group\cr \tab \cr
#' \code{counts} \tab Number of positive samples for each age group\cr \tab \cr
#' \code{age_min} \tab age_min \cr \tab \cr
#' \code{age_max} \tab age_max \cr \tab \cr
#' \code{year_init} \tab year_init \cr \tab \cr
#' \code{year_end} \tab year_end \cr \tab \cr
#' \code{tsur} \tab Year in which the survey took place \cr \tab \cr
#' \code{country} \tab The country where the survey took place \cr \tab \cr
#' \code{test} \tab The type of test taken \cr \tab \cr
#' \code{antibody} \tab antibody \cr \tab \cr
#' \code{age_mean_f} \tab Floor value of the average between age_min and age_max \cr \tab \cr
#' \code{sample_size} \tab The size of the sample \cr \tab \cr
#' \code{birth_year} \tab The year in which the individuals of each age group were bornt \cr \tab \cr
#' \code{prev_obs} \tab Observed prevalence \cr \tab \cr
#' \code{prev_obs_lower} \tab Lower limit of the confidence interval for the observed prevalence \cr \tab \cr
#' \code{prev_obs_upper} \tab Upper limit of the confidence interval for the observed prevalence \cr \tab \cr
#' }
#' The last six colums can be added to \code{seroprev_data} by means of the function \code{\link{prepare_seroprev_data}}.
#' @return data set with the binomial confidence intervals
#' @examples
#'\dontrun{
#' prepare_bin_data (seroprev_data)
#' }
#' @export
prepare_bin_data <- function(seroprev_data) {
  seroprev_data$cut_ages <-
    cut(as.numeric(seroprev_data$age_mean_f),
        seq(1, 101, by = 5),
        include.lowest = TRUE)
  xx <- seroprev_data %>%
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

# TODO: Complete the documentation of get_sim_counts
#' Function that randomly generates a sample of counts for a simulated dataset
#'
#' @param sim_data A dataframe object containing the following columns:
#' \tabular{ll}{
#' \code{birth_year} \tab List of years in which the subjects were borned \cr \tab \cr
#' \code{tsur} \tab Year of the survey\cr \tab \cr
#' \code{country} \tab Default to 'none'.\cr \tab \cr
#' \code{survey} \tab Survey label \cr \tab \cr
#' \code{age_mean_f} \tab Age \cr \tab \cr
#' }
#' @return A simulated list of counts following a binomial distribution in accordance with a given force of infection and age class sizes.
#' @examples
#'\dontrun{
#' 
#' }
#' @export
get_sim_counts <- function(sim_data, foi, size_age_class, seed = 1234){
  exposure_ages <- get_exposure_ages(sim_data)
  exposure_matrix <- get_exposure_matrix(sim_data)

  set.seed(seed = seed)
  sim_probabilities <- purrr::map_dbl(exposure_ages, ~1-exp(-dot(exposure_matrix[., ], foi)))
  sim_counts <- purrr::map_int(sim_probabilities, ~rbinom(1, size_age_class, .))

  return(sim_counts)
}

# TODO: Complete the documentation of generate_sim_data
#' Function that generates simulated data from a given Force-of-Infection
#'
#' @param foi Numeric atomic vector corresponding to the desired Force-of-Infection
#' @return Dataframe object containing the simulated data generated from \code{foi}
#' @examples
#'\dontrun{
#' size_age_class = 5
#' foi <- rep(0.02, 50)
#' sim_data <- generate_sim_data(foi = foi,
#'                               size_age_class = size_age_class,
#'                               tsur = 2050,
#'                               birth_year_min = 2000,
#'                               survey_label = 'sim_constant_foi')
#' }
#' @export
generate_sim_data <- function(foi,
                              size_age_class,
                              tsur,
                              birth_year_min,
                              survey_label,
                              test = "fake",
                              antibody = "IgG"
                              ){
    sim_data <- data.frame(birth_year = c(birth_year_min:(tsur - 1))) %>%
        mutate(tsur = tsur,
            country = 'None',
            test = test,
            antibody = antibody,
            survey = survey_label,
            age_mean_f = tsur - birth_year)
    sim_data <- sim_data %>%
        mutate(counts = get_sim_counts(sim_data, foi, size_age_class),
            total = size_age_class) %>%
        prepare_seroprev_data(add_age_mean_f = FALSE)

    return(sim_data)
}

# TODO: Complete the documentation of generate_sim_data_grouped
#' Function that generates  grouped simulated data from a given Force-of-Infection
#'
#' @param sim_data Dataframe with the structure of the output of \code{\linl{generate_sim_data}}.
#' @return Dataframe object containing grouped simulated data generated from \code{foi}
#' @examples
#'\dontrun{
#' size_age_class = 5
#' foi <- rep(0.02, 50)
#' sim_data <- generate_sim_data(foi = foi,
#'                               size_age_class = size_age_class,
#'                               tsur = 2050,
#'                               birth_year_min = 2000,
#'                               survey_label = 'sim_constant_foi')
#' sim_data_A_grouped <- generate_sim_data_grouped(sim_data = sim_data,
#'                                                 foi = foi,
#'                                                 size_age_class = size_age_class,
#'                                                 tsur = 2050,
#'                                                 birth_year_min = 2000,
#'                                                 survey_label = 'sim_constant_foi_grouped')
#' }
#' @export
generate_sim_data_grouped <- function(sim_data,
                                      foi,
                                      size_age_class,
                                      tsur,
                                      birth_year_min,
                                      survey_label,
                                      test = "fake",
                                      antibody = "IgG") {

    sim_data <- sim_data %>% mutate(age_group = 'NA', age = age_mean_f) %>% arrange(age)
    sim_data$age_group[sim_data$age > 0 & sim_data$age < 5] <-   '01-04'
    sim_data$age_group[sim_data$age > 4 & sim_data$age < 10] <-  '05-09'
    sim_data$age_group[sim_data$age > 9 & sim_data$age < 15] <-  '10-14'
    sim_data$age_group[sim_data$age > 14 & sim_data$age < 20] <- '15-19'
    sim_data$age_group[sim_data$age > 19 & sim_data$age < 25] <- '20-24'
    sim_data$age_group[sim_data$age > 24 & sim_data$age < 30] <- '25-29'
    sim_data$age_group[sim_data$age > 29 & sim_data$age < 35] <- '30-34'
    sim_data$age_group[sim_data$age > 34 & sim_data$age < 40] <- '35-39'
    sim_data$age_group[sim_data$age > 39 & sim_data$age < 45] <- '40-44'
    sim_data$age_group[sim_data$age > 44 & sim_data$age < 51] <- '45-50'


    sim_data_grouped <- sim_data %>% group_by(age_group) %>%
      dplyr::summarise(total = sum(total), counts = sum(counts)) %>%
      mutate(tsur = sim_data$tsur[1],
              country = "None",
              survey = survey_label,
              test = test,
              antibody = antibody) %>%
      mutate(age_min = as.numeric(substr(age_group, 1, 2)),
              age_max = as.numeric(substr(age_group, 4, 5))) %>%
      prepare_seroprev_data()

    return(sim_data_grouped)
}
