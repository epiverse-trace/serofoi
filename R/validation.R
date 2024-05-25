stop_if_missing <- function(serosurvey, must_have_cols) {
  if (
    !all(
      must_have_cols
      %in% colnames(serosurvey)
    )
  ) {
    missing <- must_have_cols[which(!(must_have_cols %in% colnames(serosurvey)))]
    stop(
      "The following mandatory columns in `serosurvey` are missing.\n",
      toString(missing)
    )
  }
}

stop_if_wrong_type <- function(serosurvey, col_types) {
  error_messages <- list()
  for (col in names(col_types)) {
    valid_col_types <- as.list(col_types[[col]])

    # Only validates column type if the column exists in the dataframe
    if (col %in% colnames(serosurvey) &&
        !any(vapply(valid_col_types, function(type) {
          do.call(sprintf("is.%s", type), list(serosurvey[[col]]))
        }, logical(1)))) {
      error_messages <- append(
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
}

validate_serosurvey <- function(serosurvey) {
  col_types <- list(
    age_min = "numeric",
    age_max = "numeric",
    sample_size = "numeric",
    n_seropositive = "numeric",
    tsur = "numeric"
  )

  stop_if_missing(serosurvey, must_have_cols = names(col_types))

  stop_if_wrong_type(serosurvey, col_types)

  return(serosurvey)
}
