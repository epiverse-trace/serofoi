#' Prepare data
#'
#' Function that prepares the data for modelling
#' @param model_data dataset to be processed
#' @param alpha probability of a type I error (Hmisc::binconf)
#' @return model_data with additional columns necessary for the analysis
#' @examples
#' data_test <- readRDS("data/data.RDS")
#' data_test <- prepare_data(data_test, alpha = 0.05)
#' @export
prepare_data <- function(model_data,
                         alpha = 0.05) {
  model_data <- model_data %>%
    dplyr::mutate(age_mean_f = floor((age_min + age_max) / 2), sample_size = sum(total)) %>%
    cbind(
      Hmisc::binconf(
        model_data$counts,
        model_data$total,
        alpha = alpha,
        method = "exact",
        return.df = TRUE
      )
    ) %>%
    dplyr::rename(
      prev_obs = PointEst,
      prev_obs_lower = Lower,
      prev_obs_upper = Upper
    )

  return(model_data)
}

#' Prepare data to plot binomial confidence intervals
#'
#' Function that prepares the data for modelling
#' @param model_data dataset to be processed
#' @return data set with the binomial confidence intervals
#' @examples
#' prepare_bin_data <- function(model_data)
#' @export
prepare_bin_data <- function(model_data) {
  model_data$cut_ages <-
    cut(as.numeric(model_data$age_mean_f),
        seq(1, 101, by = 5),
        include.lowest = TRUE)
  xx <- model_data %>%
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

