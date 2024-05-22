test_that("fit_seromodel constant model estimates the right force of infection", {

  foi_df <- data.frame(
    year = 1970:1999,
    foi = rep(0.03, 30)
  )
  survey_features <- data.frame(
    age_min = seq(1, 21, 10),
    age_max = seq(10, 30, 10),
    sample_size = c(1000, 1500, 2000)
  )

  simdata <- simulate_serosurvey(
    model = "time",
    foi = foi_df,
    survey_features = survey_features
  )  %>%
    mutate(
      survey = "constant_foi",
      tsur = max(foi_df$year) + 1
    ) %>%
    rename(
      total = sample_size,
      counts = n_seropositive
  ) %>% prepare_serodata()

  model_object <- fit_seromodel(
    serodata = simdata,
    foi_model = "constant",
    iter = 1000
  )

  foi_central_estimates <- get_foi_central_estimates(
    seromodel_object = model_object,
    serodata = simdata
  ) %>%
    mutate(
      # calculates tolerance as half the confidence interval size
      tol = (upper - lower)/2
    ) %>%
    left_join(foi_df, by = "year")

  expect_true(
    all(
      dplyr::near(
        foi_central_estimates$medianv,
        foi_central_estimates$foi,
        tol = max(foi_central_estimates$tol, 0.05)
      )
    )
  )
})

test_that("fit_seromodel tv_normal model estimates the right force of infection", {

  foi_df <- data.frame(
    year = 1940:1999,
    foi = rep(c(0.06, 0.03, 0.01), c(20, 20, 20))
  )
  survey_features <- data.frame(
    age_min = seq(1, 51, 10),
    age_max = seq(10, 60, 10),
    sample_size = 100
  )

  simdata <- simulate_serosurvey(
    model = "time",
    foi = foi_df,
    survey_features = survey_features
  )  %>%
    mutate(
      survey = "sw_dec_foi",
      tsur = max(foi_df$year) + 1
    ) %>%
    rename(
      total = sample_size,
      counts = n_seropositive
    ) %>% prepare_serodata()

  model_object <- fit_seromodel(
    serodata = simdata,
    foi_model = "tv_normal",
    chunks = rep(c(1, 2, 3), c(15, 20, 20)),
    iter = 800
  )

  foi_central_estimates <- get_foi_central_estimates(
    seromodel_object = model_object,
    serodata = simdata
  ) %>%
    mutate(
      # calculates tolerance as half the confidence interval size
      tol = (upper - lower) / 2
    ) %>%
    left_join(foi_df, by = "year")

  expect_true(
    all(
      dplyr::near(
        foi_central_estimates$medianv,
        foi_central_estimates$foi,
        tol = max(foi_central_estimates$tol, 0.05)
      )
    )
  )
})

test_that("fit_seromodel tv_normal_log model estimates the right force of infection", {

  foi_df <- data.frame(
    year = 1950:1999,
    foi = rep(
      c(0.001, 0.4, 0.001),
      c(30, 5, 15)
    )
  )

  survey_features <- data.frame(
    age_min = seq(1, 41, 10),
    age_max = seq(10, 50, 10),
    sample_size = 100
  )

  simdata <- simulate_serosurvey(
    model = "time",
    foi = foi_df,
    survey_features = survey_features
  )  %>%
    mutate(
      survey = "big_outbreak",
      tsur = max(foi_df$year) + 1
    ) %>%
    rename(
      total = sample_size,
      counts = n_seropositive
    ) %>% prepare_serodata()

  model_object <- fit_seromodel(
    serodata = simdata,
    foi_model = "tv_normal_log",
    chunks = rep(c(1, 2, 3), c(25, 5, 15)),
    iter = 700
  )

  foi_central_estimates <- get_foi_central_estimates(
    seromodel_object = model_object,
    serodata = simdata
  ) %>%
    mutate(
      # calculates tolerance as half the confidence interval size
      tol = (upper - lower) / 2
    ) %>%
    left_join(foi_df, by = "year")

  expect_true(
    all(
      dplyr::near(
        foi_central_estimates$medianv,
        foi_central_estimates$foi,
        tol = max(foi_central_estimates$tol, 0.05)
      )
    )
  )
})
