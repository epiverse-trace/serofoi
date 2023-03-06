test_that("compilation", {
    # So far we are skipping tests on these platforms until
    # we find an efficient way to update rstan testthat snapshots on all of them
    skip_on_os(c("windows", "mac"))
    source("testing_utils.R")

    set.seed(1234) # For reproducibility


    # library(devtools)
    library(dplyr)

    mydata <- readRDS(testthat::test_path("extdata", "data.RDS"))

    # Modelling module functions
    seroprev_data <- preprare_seroprev_data(seroprev_data = mydata, alpha = 0.05)

    exposure_years <- get_exposure_years(seroprev_data)

    exposure_matrix <- get_exposure_matrix(
        seroprev_data = seroprev_data,
        exposure_years = exposure_years
    )

    stan_model <- save_or_load_model(seroprev_model_name = "constant_foi_bi")

    fit_seroprev_model_test <- fit_seroprev_model(
        seroprev_data = seroprev_data,
        seroprev_model_name = "constant_foi_bi",
        n_iters = 1000,
        n_thin = 2,
        delta = 0.90,
        m_treed = 10,
        decades = 0
    )


    model_object <- run_seroprev_model(seroprev_data = seroprev_data, seroprev_model_name = "constant_foi_bi")

    model_summary <- extract_seroprev_model_summary(model_object)

    foi <- rstan::extract(model_object$fit, "foi", inc_warmup = FALSE)[[1]]
    expanded_prevalence <- get_prev_expanded(foi, seroprev_data)

    column_comparation_functions <- list(
        age = equal_exact(),
        predicted_prev = equal_with_tolerance(),
        predicted_prev_lower = equal_with_tolerance(),
        predicted_prev_upper = equal_with_tolerance(),
        prev_obs = equal_with_tolerance(),
        prev_obs_lower = equal_with_tolerance(),
        prev_obs_upper = equal_with_tolerance(),
        sample_by_age = equal_exact(), #
        positives = equal_with_tolerance(),
        survey = equal_exact(),
        cut_ages = equal_exact(), #
        bin_size = equal_with_tolerance(),
        bin_pos = equal_with_tolerance(),
        p_obs_bin = equal_with_tolerance(),
        p_obs_bin_l = equal_with_tolerance(),
        p_obs_bin_u = equal_with_tolerance()
    )



    expect_similar_dataframes(
        "expanded_prevalence", expanded_prevalence, column_comparation_functions
    )
})
