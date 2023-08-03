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
#' \code{birth_year} \tab The year in which the individuals of each age group were bornt \cr \tab \cr
#' \code{prev_obs} \tab Observed prevalence \cr \tab \cr
#' \code{prev_obs_lower} \tab Lower limit of the confidence interval for the observed prevalence \cr \tab \cr
#' \code{prev_obs_upper} \tab Upper limit of the confidence interval for the observed prevalence \cr \tab \cr
#' }
#' @examples
#'\dontrun{
#' data("serodata")
#' data_test <- prepare_serodata(serodata)
#' }
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
#' \code{birth_year} \tab The year in which the individuals of each age group were bornt \cr \tab \cr
#' \code{prev_obs} \tab Observed prevalence \cr \tab \cr
#' \code{prev_obs_lower} \tab Lower limit of the confidence interval for the observed prevalence \cr \tab \cr
#' \code{prev_obs_upper} \tab Upper limit of the confidence interval for the observed prevalence \cr \tab \cr
#' }
#' The last six colums can be added to \code{serodata} by means of the function \code{\link{prepare_serodata}}.
#' @return data set with the binomial confidence intervals
#' @examples
#'\dontrun{
#' prepare_bin_data (serodata)
#' }
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
#' \code{birth_year} \tab List of years in which the subjects were borned \cr \tab \cr
#' \code{tsur} \tab Year of the survey\cr \tab \cr
#' \code{country} \tab Default to 'none'.\cr \tab \cr
#' \code{survey} \tab Survey label \cr \tab \cr
#' \code{age_mean_f} \tab Age \cr \tab \cr
#' }
#' @param foi Numeric atomic vector corresponding to the desired Force-of-Infection.
#' @param seed The seed for random number generation.
#' @return A simulated list containing a seropositivity distribution by age for given simulated
#' dataset and desired foi trend.
#' @examples
#'\dontrun{
#'
#' }
#' @export
get_sim_prob <- function(sim_data, foi, seed = 1234) {
  exposure_ages <- get_exposure_ages(sim_data)
  exposure_matrix <- get_exposure_matrix(sim_data)

  set.seed(seed = seed)
  sim_probabilities <- purrr::map_dbl(exposure_ages, ~1-exp(-pracma::dot(exposure_matrix[., ], foi)))

  return(sim_probabilities)
}

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
#' @param foi Numeric atomic vector corresponding to the desired Force-of-Infection.
#' @param size_age_class Size of each age group specified by either an atomic
#' vector of the same size as \code{foi} or an integer.
#' This corresponds to the number of trials \code{size} in \link[stats]{rbinom}.
#' @param seed The seed for random number generation.
#' @return A simulated list of counts following a binomial distribution in accordance with a given 
#' force of infection and age class sizes.
#' @examples
#'\dontrun{
#'
#' }
#' @export
get_sim_counts <- function(sim_data, foi, size_age_class, seed = 1234) {
  sim_probabilities <- get_sim_prob(sim_data = sim_data, foi = foi, seed = 1234)
  sim_counts <- purrr::map_int(sim_probabilities, ~rbinom(1, size_age_class, .))

  return(sim_counts)
}

#' Function that generates simulated positive counts assuming a known historical force of infection
#'
#' @param foi Numeric atomic vector corresponding to the desired Force-of-Infection
#' @param size_age_class Size of each age group specified by either an atomic
#' vector of the same size as \code{foi} or an integer.
#' This corresponds to the number of trials \code{size} in \link[stats]{rbinom}.
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
    sim_data <- sim_data %>%
        mutate(counts = get_sim_counts(sim_data, foi, size_age_class, seed = seed),
            total = size_age_class) %>%
        prepare_serodata(add_age_mean_f = FALSE)

    return(sim_data)
}

#' Method for constructing age-group variable from age column
#'
#' This function was taken from \link[vaccineff]{get_age_group}.
#' This method splits an age interval from min_val to max_val into
#' (max_val-min_val)/step intervals.
#' By default min_val is set 0, however it can be assigned by
#' convenience.
#' If the method finds ages greater or equal than max_val
#' it assigns the string ">{max_val}".
#' To avoid errors it is necessary to set step < max_val.
#' It is also suggested to choose the step such
#' that max_val%%(step+1) = 0.
#' @param data dataset with at least a column containing the age
#' information
#' @param col_age name of the column containing the age
#' information
#' @param  max_val maximum value of age interval to split
#' @param  min_val minimum value of age interval to split
#' @param  step step used to split the age interval
#' @return age_group
#' @examples
#' \dontrun{
#' sim_data <- generate_sim_data(foi = rep(0.02, 50),
#'                               size_age_class = 5,
#'                               tsur = 2050,
#'                               birth_year_min = 2000,
#'                               survey_label = 'sim_constant_foi')
#' sim_data$age_group <- get_age_group(data = sim_data,
#'                                     col_age = "age_mean_f",
#'                                     max_val = max(sim_data$age_mean_f),
#'                                     min_val = min(sim_data$age_mean_f),
#'                                     step = 5)
#' }
#' @export
get_age_group <- function(data, col_age, max_val, min_val, step) {
  n_steps <- as.integer((max_val - min_val) / step) + 1
  limits_low <- c(as.integer(seq(min_val,
                                 max_val,
                                 length.out = n_steps)))
  limits_hgh <- limits_low + step
  lim_labels <- paste(as.character(limits_low), as.character(limits_hgh),
                      sep = "-")
  lim_labels[length(lim_labels)] <- paste0("+",
                                           limits_low[length(limits_low)])
  lim_breaks <- c(-Inf, limits_low[2:length(limits_low)] - 1, Inf)

  age_group <- cut(data[[col_age]],
                   breaks = lim_breaks,
                   labels = lim_labels)
  return(age_group)
}

# TODO: Complete the documentation of group_sim_data
#' Function that generates grouped simulated data by age from a given Force-of-Infection
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
#' sim_data_grouped <- group_sim_data(sim_data = sim_data,
#'                                    foi = foi,
#'                                    size_age_class = size_age_class,
#'                                    tsur = 2050,
#'                                    birth_year_min = 2000,
#'                                    survey_label = 'sim_constant_foi_grouped')
#' }
#' @export
group_sim_data <- function(sim_data,
                          foi,
                          size_age_class,
                          tsur,
                          birth_year_min,
                          survey_label,
                          test = "fake",
                          antibody = "IgG",
                          seed = 1234) {

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
      prepare_serodata()

    return(sim_data_grouped)
}