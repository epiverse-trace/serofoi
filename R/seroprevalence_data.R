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
