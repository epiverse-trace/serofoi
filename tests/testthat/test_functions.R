test_that("compilation", {
    source("testing_utils.R")

    set.seed(1234) # For reproducibility


    library(devtools)
    library(dplyr)

    mydata <- readRDS(test_path("extdata", "data.RDS"))

    # Modelling module functions
    model_data <- prepare_data(model_data = mydata, alpha = 0.05)

    exposure_years <- get_exposure_years(model_data)

    exposure_matrix <- get_exposure_matrix(
        model_data = model_data,
        exposure_years = exposure_years
    )

    stan_model <- save_or_load_model(model_name = "constant_foi_bi")

    fit_model_test <- fit_model(
        model_data = model_data,
        model_name = "constant_foi_bi",
        n_iters = 1000,
        n_thin = 2,
        delta = 0.90,
        m_treed = 10,
        decades = 0
    )


    model_object <- run_model(model_data = model_data, model_name = "constant_foi_bi")

    model_summary <- extract_model_summary(model_object)

    foi <- rstan::extract(model_object$fit, "foi", inc_warmup = FALSE)[[1]]
    expanded_prevalence <- get_prev_expanded(foi, model_data)

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



    expect_true(
        compare_dataframes(
            "expanded_prevalence", expanded_prevalence, column_comparation_functions
        )
    )
})
