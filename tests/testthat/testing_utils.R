# TODO Move to separate package ###
# TODO Document all functions and provide examples
library(testthat)
library(vdiffr)
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

# TODO use testthat snapshots
expect_similar_dataframes <- function(name, actual_df, column_comparation_functions) {
    actual_df_filename <- file.path(tempdir(), paste(name, "csv", sep = "."))
    write.csv(actual_df, actual_df_filename)
    compare_fun <- function(expected_df_filename, actual_df_filename) {
        return(compare_dataframes(expected_df_filename, actual_df_filename, column_comparation_functions))
    }
    expect_snapshot_file(actual_df_filename, compare = compare_fun)
}


compare_dataframes <- function(expected_df_filename, actual_df_filename, column_comparation_functions) {
    expected_df <- read.csv(expected_df_filename)
    actual_df <- read.csv(actual_df_filename)

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
}

expect_same_plot <- function(plot_name, actual_plot) {
    if (per_platform_snapshots) {
        title <- file.path(r_version_id(), plot_name)
    } else {
        title <- plot_name
    }
    return(vdiffr::expect_doppelganger(title, actual_plot))
}
