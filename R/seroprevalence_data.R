#' Prepare data
#'
#' Function that prepares the data for modelling
#' @param model_data dataset to be processed
#' @param alpha probability of a type I error (Hmisc::binconf)
#' @return model_data with additional columns necessary for the analysis
#' @export
prepare_data <- function(model_data,
                         alpha = 0.05) {
  model_data <- model_data %>%
    dplyr::mutate(age_mean_f = floor((age_min+age_max)/2), sample_size = sum(total)) %>%
    cbind(Hmisc::binconf(model_data$counts,
                         model_data$total,
                         alpha = alpha,
                         method="exact",
                         return.df = TRUE)) %>%
    dplyr::rename(prev_obs = PointEst, prev_obs_lower = Lower, prev_obs_upper = Upper)

  return(model_data)
}


