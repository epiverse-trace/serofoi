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


equal_dataframes <- function(df1, df2, column_comparation_functions) {
    all_columns_ok <- TRUE
    for (col in base::names(column_comparation_functions)) {
        compare_function <- column_comparation_functions[[col]]
        all_columns_ok <- all_columns_ok &&
            compare_function(df1[col], df2[col])
    }
    return(all_columns_ok)
}
