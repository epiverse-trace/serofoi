# TODO Move to separate package ###
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
