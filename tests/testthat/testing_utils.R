# TODO Move to separate package ###
library(testthat)
library(ggplot2)
library(png)
equal_with_tolerance <- function(tolerance = 1e-4) {
    function(a, b, tolerance = 1e-4) {
        x <- base::mapply(function(x, y) abs(x - y), a, b)
        return(base::all(x < tolerance))
    }
}
equal_exact <- function() {
    function(a, b) {
        x <- base::mapply(function(x, y) x == y, a, b)
        return(base::all(x == TRUE))
    }
}


compare_dataframes <- function(df1, df2, column_comparation_functions) {
    all_columns_ok <- TRUE
    for (col in base::names(column_comparation_functions)) {
        if (col %in% colnames(df1) && col %in% colnames(df2)) {
            compare_function <- column_comparation_functions[[col]]
            all_columns_ok <- all_columns_ok &&
                compare_function(df1[col], df2[col])
        } else {
            if (!(col %in% colnames(df1))) {
                cat("Column", col, "not in first dataframe")
            }
            if (!(col %in% colnames(df1))) {
                cat("Column", col, "not in second dataframe")
            }
        }
    }
    return(all_columns_ok)
}

# withr
compare_plots <- function(expected_plot_name, actual_plot) {
    base_path <- test_path("extdata", "plots")
    actual_plot_filename <- paste(file.path(base_path, "actual", expected_plot_name), "png", sep = ".")
    print(paste("Saving actual plot image:", actual_plot_filename))
    ggsave(filename = actual_plot_filename, plot = actual_plot)
    actual_plot_image <- readPNG(actual_plot_filename)

    expected_plot_filename <- paste(file.path(base_path, "expected", expected_plot_name), "png", sep = ".")
    if (file.exists(expected_plot_filename)) {
        expected_plot_image <- readPNG(expected_plot_filename)
        return(identical(expected_plot_image, actual_plot_image))
    } else {
        # If expected plot file does not exist, creates one based on the actual plot image
        # This behavior is for convenience and should only be used the first time the test is executed,
        # since it will create the data that will be used as the expected plot in further runs
        ggsave(filename = expected_plot_filename, plot = actual_plot)
        warning(
            "No expected plot image found: ", expected_plot_filename, "\n",
            "Expected image will be created from the actual plot"
        )
        return(TRUE)
    }
}
