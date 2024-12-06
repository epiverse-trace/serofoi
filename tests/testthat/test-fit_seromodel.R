# Setup for testing ----
tol_min = 1e-3
stan_seed = "123"

survey_features <- data.frame(
  age_min = c(1, 5, 15),
  age_max = c(4, 14, 20),
  n_sample = c(500, 500, 500))

mu <- 0.1 # seroreversion_rate

test_foi_estimation <- function(seromodel, serosurvey, foi) {
  foi_estimates <- extract_central_estimates(
    seromodel = seromodel,
    serosurvey = serosurvey,
    par_name = "foi_expanded"
    ) |>
    dplyr::mutate(tol = pmax((upper - lower)/2, tol_min))
  expect_true(
    all(
      dplyr::near(
        foi$foi,
        foi_estimates$median,
        tol = foi_estimates$tol
      )
    )
  )}

test_serorev_estimation <- function(seromodel, serosurvey, mu) {
  mu_estimate <- extract_central_estimates(
    seromodel = seromodel,
    serosurvey = serosurvey,
    par_name = "seroreversion_rate"
  ) |>
    dplyr::mutate(tol = pmax((upper - lower)/2, tol_min))

  expect_true(
    dplyr::near(
      mu,
      mu_estimate$median,
      tol = mu_estimate$tol
    )
  )}

# Test for add_age_group_to_serosurvey ----
test_that("add_age_group_to_serosurvey handles existing age_group column", {
  # Case where serosurvey already has an age_group column
  serosurvey_with_age_group <- dplyr::mutate(
    survey_features,
    age_group = c(2, 10, 18)
    )

  expect_message(
    result <- add_age_group_to_serosurvey(serosurvey_with_age_group),
    "Using `age_group` already present in serosurvey"
  )

  # Check that the existing age_group column is retained
  expect_equal(result$age_group, serosurvey_with_age_group$age_group)
})

test_that("add_age_group_to_serosurvey creates age_group column if missing", {
  result <- add_age_group_to_serosurvey(survey_features)

  # Check that age_group column was created correctly
  expected_age_group <- floor((survey_features$age_min + survey_features$age_max) / 2)
  expect_equal(result$age_group, expected_age_group)
})

# Test fit_seromodel ----
test_that("fit_seromodel correctly estimates constant foi using default settings", {
  skip_on_cran()

  # constant FOI
  foi <- data.frame(
    year = seq(1990, 2009, 1),
    foi = rep(0.01, 20)
  )

  # no seroreversion
  set.seed(123)
  serosurvey <- simulate_serosurvey(
    model = "time",
    foi = foi,
    survey_features = survey_features
  )
  set.seed(Sys.time())
  seromodel <- suppressWarnings(
    fit_seromodel(serosurvey, model_type = "constant", seed = stan_seed)
    )

  test_foi_estimation(seromodel, serosurvey, foi)

  # seroreversion
  set.seed(123)
  serosurvey <- simulate_serosurvey(
    model = "time",
    foi = foi,
    survey_features = survey_features,
    seroreversion_rate = mu
  ) |>
  dplyr::mutate(survey_year = 2050)
  set.seed(Sys.time())

  seromodel <- suppressWarnings(
    fit_seromodel(
      serosurvey,
      model_type = "constant",
      is_seroreversion = TRUE,
      seed = stan_seed)
  )

  test_foi_estimation(seromodel, serosurvey, foi)
  test_serorev_estimation(seromodel, serosurvey, mu)
  })

test_that("fit_seromodel correctly estimates time-varying foi using default priors", {
  skip_on_cran()

  foi <- data.frame(
    year = seq(1990, 2009, 1),
    foi = c(rep(0.01, 10), rep(0.005, 10))
  )

  # no seroreversion
  set.seed(123)
  serosurvey <- simulate_serosurvey(
    model = "time",
    foi = foi,
    survey_features = survey_features
  ) |>
  dplyr::mutate(survey_year = 2050)
  set.seed(Sys.time())

  foi_index <- get_foi_index(serosurvey, group_size = 10, model_type = "time")
  seromodel <- suppressWarnings(
    fit_seromodel(
      serosurvey,
      model_type = "time",
      foi_index = foi_index,
      seed = stan_seed
      )
    )

  test_foi_estimation(seromodel, serosurvey, foi)

  # seroreversion
  set.seed(123)
  serosurvey <- simulate_serosurvey(
    model = "time",
    foi = foi,
    survey_features = survey_features,
    seroreversion_rate = mu
  ) |>
  dplyr::mutate(survey_year = 2050)
  set.seed(Sys.time())

  foi_index <- get_foi_index(serosurvey, group_size = 10, model_type = "time")
  seromodel <- suppressWarnings(
    fit_seromodel(
      serosurvey,
      model_type = "time",
      foi_index = foi_index,
      foi_prior = sf_uniform(),
      is_seroreversion = TRUE,
      seroreversion_prior = sf_normal(mu, mu/10),
      seed = stan_seed
    )
  )

  test_foi_estimation(seromodel, serosurvey, foi)
  test_serorev_estimation(seromodel, serosurvey, mu)
  })

test_that("fit_seromodel correctly estimates age-varying foi", {
  skip_on_cran()

  # age-varying FOI
  foi <- data.frame(
    age = seq(1, 20, 1),
    foi = c(rep(0.0, 10), rep(0.01, 10))
  )

  # no seroreversion
  set.seed(123)
  serosurvey <- simulate_serosurvey(
    model = "age",
    foi = foi,
    survey_features = survey_features
  )
  set.seed(Sys.time())

  foi_index <- get_foi_index(serosurvey, group_size = 10, model_type = "age")
  seromodel <- suppressWarnings(
    fit_seromodel(
      serosurvey,
      model_type = "age",
      foi_index = foi_index,
      seed = stan_seed
    )
  )

  test_foi_estimation(seromodel, serosurvey, foi)

  # seroreversion
  set.seed(123)
  serosurvey <- simulate_serosurvey(
    model = "age",
    foi = foi,
    survey_features = survey_features,
    seroreversion_rate = mu
  )
  set.seed(Sys.time())

  foi_index <- get_foi_index(serosurvey, group_size = 10, model_type = "age")
  seromodel <- suppressWarnings(
    fit_seromodel(
      serosurvey,
      model_type = "age",
      foi_index = foi_index,
      foi_prior = sf_normal(0, 1e-4),
      is_seroreversion = TRUE,
      seroreversion_prior = sf_normal(mu, mu/10),
      seed = stan_seed
    )
  )

  test_foi_estimation(seromodel, serosurvey, foi)
  test_serorev_estimation(seromodel, serosurvey, mu)
  })

test_that("fit_seromodel correctly identifies outbreak using time-log-foi model", {
  skip_on_cran()

  # time-varying FOI
  outbreak_years <- 2
  foi <- data.frame(
    year = seq(1990, 2009, 1),
    foi = c(rep(0, 10), rep(0.5, outbreak_years), rep(0.01, 10 - outbreak_years))
  )

  # no seroreversion
  set.seed(123)
  serosurvey <- simulate_serosurvey(
    model = "time",
    foi = foi,
    survey_features = survey_features
  )
  set.seed(Sys.time())

  foi_index <- data.frame(
    year = foi$year,
    foi_index = c(
      rep(1, 10),
      rep(2, outbreak_years),
      rep(3, 10 - outbreak_years)
    )
  )
  seromodel <- suppressWarnings(
    fit_seromodel(
      serosurvey,
      model_type = "time",
      is_log_foi = TRUE,
      foi_prior = sf_normal(0, 1e-4),
      foi_index = foi_index,
      seed = stan_seed
    )
  )

  test_foi_estimation(seromodel, serosurvey, foi)

  # seroreversion
  set.seed(123)
  serosurvey <- simulate_serosurvey(
    model = "time",
    foi = foi,
    survey_features = survey_features,
    seroreversion_rate = mu
  ) |>
  dplyr::mutate(survey_year = 2050)
  set.seed(Sys.time())

  foi_index <- data.frame(
    year = foi$year,
    foi_index = c(
      rep(1, 10),
      rep(2, outbreak_years),
      rep(3, 10 - outbreak_years)
    )
  )
  seromodel <- suppressWarnings(
    fit_seromodel(
      serosurvey,
      model_type = "time",
      is_log_foi = TRUE,
      foi_prior = sf_normal(0, 1e-4),
      foi_index = foi_index,
      is_seroreversion = TRUE,
      seroreversion_prior = sf_normal(mu, mu/10),
      seed = stan_seed
    )
  )

  test_foi_estimation(seromodel, serosurvey, foi)
  test_serorev_estimation(seromodel, serosurvey, mu)
  })

