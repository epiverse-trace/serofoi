# TODO Move to separate package ###
library(testthat)
equal_with_tolerance <- function(tolerance = 1e-2) {
    function(a, b) {
        c <- base::mapply(function(x, y) {
            if (is.na(x) && is.na(y)) {
                return(0)
            } else if (is.na(x) || is.na(y)) {
                return(tolerance + 1)
            } else {
                return(abs(x - y))
            }
        }, a, b)
        return(base::all(c < tolerance))
    }
}
equal_exact <- function() {
    function(a, b) {
        x <- base::mapply(function(x, y) x == y || (is.na(x) && is.na(y)), a, b)
        return(base::all(x == TRUE))
    }
}



compare_dataframes <- function(expected_df_name, actual_df, column_comparation_functions) {
    base_path <- test_path("extdata", "dataframes") # TODO move to config.yml
    actual_df_filename <- paste(file.path(base_path, "actual", expected_df_name), "csv", sep = ".")
    expected_df_filename <- paste(file.path(base_path, "expected", expected_df_name), "csv", sep = ".")

    if (file.exists(expected_df_filename)) {
        expected_df <- read.csv(expected_df_filename)
        all_columns_ok <- TRUE
        for (col in base::names(column_comparation_functions)) {
            if (col %in% colnames(expected_df) && col %in% colnames(actual_df)) {
                compare_function <- column_comparation_functions[[col]]
                col_ok <- compare_function(expected_df[[col]], actual_df[[col]])
                if (!col_ok) {
                    cat("Column", col, "differs ", expected_df[[col]], "!=", actual_df[[col]], "\n")
                }
                all_columns_ok <- all_columns_ok && col_ok
            } else {
                if (!(col %in% colnames(expected_df))) {
                    cat("Column", col, "not in first dataframe")
                }
                if (!(col %in% colnames(expected_df))) {
                    cat("Column", col, "not in second dataframe")
                }
            }
        }
        return(all_columns_ok)
    } else {
        # If expected dataframe file does not exist, creates one based on the actual dataframe
        # This behavior is for convenience and should only be used the first time the test is executed,
        # since it will create the data that will be used as the expected dataframe in further runs
        warning(
            "No expected dataframe found: ", expected_df_filename, "\n",
            "Expected dataframe will be created from the actual dataframe"
        )
        dir.create(dirname(expected_df_filename), recursive = TRUE)
        write.csv(actual_df, file = expected_df_filename)
        return(TRUE)
    }
}
