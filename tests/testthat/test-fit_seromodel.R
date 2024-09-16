# Setup for testing ----
tol_min = 1e-3
tol_max = 0.1
stan_seed = "123"

survey_features <- data.frame(
  age_min = c(1, 3, 15),
  age_max = c(2, 14, 20),
  n_sample = c(1000, 2000, 1500))

mu <- 0.01 # seroreversion_rate

test_foi_estimation <- function(seromodel, serosurvey, foi) {
  foi_estimates <- extract_central_estimates(
    seromodel = seromodel,
    serosurvey = serosurvey,
    par_name = "foi_expanded"
    ) %>%
    mutate(tol = pmax(pmin((upper - lower)/2, tol_max), tol_min))

  expect_true(
    all(
      dplyr::near(
        foi$foi,
        foi_estimates$median,
        tol = foi_estimates$tol
      )
    )
  )
}

test_serorev_estimation <- function(seromodel, serosurvey, mu) {
  mu_estimate <- extract_central_estimates(
    seromodel = seromodel,
    serosurvey = serosurvey,
    par_name = "seroreversion_rate"
  ) %>%
    mutate(tol = pmax((upper - lower)/2, tol_min))

  expect_true(
    dplyr::near(
      mu,
      mu_estimate$median,
      tol = mu_estimate$tol
    )
  )
}

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
  )
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
  }
)

test_that("fit_seromodel correctly estimates time-varying foi using default priors", {
  skip_on_cran()

  # no seroreversion
  foi <- data.frame(
    year = seq(1990, 2009, 1),
    foi = c(rep(0.01, 10), rep(0, 10))
  )

  set.seed(123)
  serosurvey <- simulate_serosurvey(
    model = "time",
    foi = foi,
    survey_features = survey_features
  )
  set.seed(Sys.time())

  seromodel <- suppressWarnings(
    fit_seromodel(
      serosurvey,
      model_type = "time",
      foi_index = get_foi_index(serosurvey, group_size = 10),
      seed = stan_seed
      )
    )

  test_foi_estimation(seromodel, serosurvey, foi)

  # seroreversion
  foi <- data.frame(
    year = seq(1990, 2009, 1),
    foi = c(rep(0, 10), rep(0.01, 10))
  )

  set.seed(123)
  serosurvey <- simulate_serosurvey(
    model = "time",
    foi = foi,
    survey_features = survey_features,
    seroreversion_rate = mu
  )
  set.seed(Sys.time())

  seromodel <- suppressWarnings(
    fit_seromodel(
      serosurvey,
      model_type = "time",
      foi_index = get_foi_index(serosurvey, group_size = 10),
      foi_prior = sf_normal(0, 1e-4),
      is_seroreversion = TRUE,
      seroreversion_prior = sf_normal(mu, mu/10),
      seed = stan_seed
    )
  )

  test_foi_estimation(seromodel, serosurvey, foi)
  test_serorev_estimation(seromodel, serosurvey, mu)
  }
)

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

  seromodel <- suppressWarnings(
    fit_seromodel(
      serosurvey,
      model_type = "age",
      foi_index = get_foi_index(serosurvey, group_size = 10),
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

  seromodel <- suppressWarnings(
    fit_seromodel(
      serosurvey,
      model_type = "age",
      foi_index = get_foi_index(serosurvey, group_size = 10),
      foi_prior = sf_normal(0, 1e-4),
      is_seroreversion = TRUE,
      seroreversion_prior = sf_normal(mu, mu/10),
      seed = stan_seed
    )
  )

  test_foi_estimation(seromodel, serosurvey, foi)
  test_serorev_estimation(seromodel, serosurvey, mu)
  }
)

test_that("fit_seromodel correctly identifies outbreak using time-log-foi model", {
  skip_on_cran()

  # time-varying FOI
  foi <- data.frame(
    year = seq(1990, 2009, 1),
    foi = c(rep(0, 10), rep(0.1, 5), rep(0,5))
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
    fit_seromodel(
      serosurvey,
      model_type = "time",
      is_log_foi = TRUE,
      foi_prior = sf_normal(0, 1e-4),
      foi_index = c(rep(1, 10), rep(2, 5), rep(3, 5)),
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
  )
  set.seed(Sys.time())

  seromodel <- suppressWarnings(
    fit_seromodel(
      serosurvey,
      model_type = "time",
      is_log_foi = TRUE,
      foi_prior = sf_normal(0, 1e-4),
      foi_index = c(rep(1, 10), rep(2, 5), rep(3, 5)),
      is_seroreversion = TRUE,
      seroreversion_prior = sf_normal(mu, mu/10),
      seed = stan_seed
    )
  )

  test_foi_estimation(seromodel, serosurvey, foi)
  test_serorev_estimation(seromodel, serosurvey, mu)
  }
)

