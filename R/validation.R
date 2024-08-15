validate_serosurvey <- function(serosurvey) {
  # Check that necessary columns are present
  col_types <- list(
    age_min = "numeric",
    age_max = "numeric",
    sample_size = "numeric",
    n_seropositive = "numeric",
    survey_year = "numeric"
  )

  checkmate::assert_names(names(serosurvey), must.include = names(col_types))

  # Validates column types
  error_messages <- list()
  for (col in names(col_types)) {
    valid_col_types <- as.list(col_types[[col]])

    # Only validates column type if the column exists in the dataframe
    if (col %in% colnames(serosurvey) &&
        !any(vapply(valid_col_types, function(type) {
          do.call(sprintf("is.%s", type), list(serosurvey[[col]]))
        }, logical(1)))) {
      error_messages <- c(
        error_messages,
        sprintf(
          "`%s` must be of any of these types: `%s`",
          col, toString(col_types[[col]])
        )
      )
    }
  }
  if (length(error_messages) > 0) {
    stop(
      "The following columns in `serosurvey` have wrong types:\n",
      toString(error_messages)
    )
  }

  return(serosurvey)
}

validate_survey_features <- function(survey_features) {

  if (!is.data.frame(survey_features) ||
      !all(
        c("age_min", "age_max", "sample_size") %in% names(survey_features))
      ) {
    stop(
      "survey_features must be a dataframe with columns ",
      "'age_min', 'age_max', and 'sample_size'."
      )
  }

  # check that the age_max of a bin does not coincide with
  # the age min of a different bin
  is_age_ok <- check_age_constraints(survey_features)
  if (!is_age_ok)
    stop(
      "Age bins in a survey are inclusive of both bounds, ",
      "so the age_max of one bin cannot equal the age_min of another."
      )
}

validate_foi_df <- function(foi_df, cnames_additional) {
  if (
    !is.data.frame(foi_df) ||
    !all(cnames_additional %in% names(foi_df)) ||
    ncol(foi_df) != (1 + length(cnames_additional))
    ) {
    if (length(cnames_additional) == 1) {
      message_end <- paste0(" and ", cnames_additional, ".")
    } else {
      message_end <- paste0(
        ", ", paste(cnames_additional, collapse = " and "), "."
        )
      message_beginning <- "foi must be a dataframe with columns foi"
      stop(glue::glue("{message_beginning}", "{message_end}"))
    }
  }
}

validate_seroreversion_rate <- function(seroreversion_rate) {
  if (!is.numeric(seroreversion_rate) || seroreversion_rate < 0) {
    stop("seroreversion_rate must be a non-negative numeric value.")
  }
}

validate_survey_and_foi_consistency <- function( #nolint
    survey_features,
    foi_df
) {
  max_age_foi_df <- nrow(foi_df)
  if (max_age_foi_df > max(survey_features$age_max))
    stop(
      "maximum age implicit in foi_df should ",
      "not exceed max age in survey_features."
      )
}

validate_survey_and_foi_consistency_age_time <- function( #nolint
    survey_features,
    foi_df
) {
  max_age_foi_df <- max(foi_df$year) - min(foi_df$year) + 1
  if (max_age_foi_df > max(survey_features$age_max))
    stop(
      "maximum age implicit in foi_df should ",
      "not exceed max age in survey_features."
      )
}
