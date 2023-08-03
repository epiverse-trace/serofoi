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
  checkmate::assert_numeric(alpha, lower = 0, upper = 1)
  checkmate::assert_logical(add_age_mean_f)
  #Check that serodata has the right columns
  stopifnot("serodata must contain the right columns" =
            all(c("survey", "total", "counts", "age_min", "age_max", "tsur",
                  "country","test","antibody"
                  ) %in%
                  colnames(serodata)
                )
            )
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