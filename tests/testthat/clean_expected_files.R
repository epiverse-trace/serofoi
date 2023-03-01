# This script removes the svg and csv files used in the automatic tests
# This is usually required when a new version of rstan is released
# Random number generation changes for each new version of
# Rstan, https://mc-stan.org/docs/2_18/reference-manual/reproducibility-chapter.html
# thus all files with expected results (both tables and graphs) become
# obsolete and need to be regenerated

library(testthat)
paths <- c(test_path("_snaps", "*"), test_path("extdata", "dataframes", "expected", "*"))

for (path in paths) {
    cat("Deleting", path, "\n")
    unlink(path, recursive = TRUE)
}
