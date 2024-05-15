test_that("probability_seropositive_time_model_by_age works", {

  foi <- data.frame(
    year=seq(1990, 2009, 1)
  ) %>%
    mutate(foi=rnorm(20, 0.2, 0.01))

  seroreversion <- 0.0
  prob_df <- probability_seropositive_time_model_by_age(
    foi = foi,
    seroreversion = seroreversion
  )

  # check output dimensions
  expect_equal(nrow(prob_df), nrow(foi))
  ages <- seq(1, nrow(foi), 1)
  expect_equal(ages, prob_df$age)

  # checking monotonicity
  derivative_foi <- diff(prob_df$seropositivity)
  expect_true(all(derivative_foi > 0))

  seroreversion <- 0.1
  prob_df_1 <- probability_seropositive_time_model_by_age(
    foi = foi,
    seroreversion = seroreversion
  )

  # check output dimensions
  expect_equal(nrow(prob_df_1), nrow(foi))
  expect_equal(ages, prob_df_1$age)

  # check seropositivities always lower (due to seroreversion)
  expect_true(all(prob_df_1$seropositivity < prob_df$seropositivity))
})
