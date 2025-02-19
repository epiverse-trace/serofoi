# Setup for testing ----
serosurvey <- data.frame(
  age_min = c(1, 6, 11, 16, 21),
  age_max = c(5, 10, 15, 20, 25),
  n_sample = c(100, 200, 300, 100, 200),
  n_seropositive = c(10, 70, 200, 80, 170)
)

# Test for sf_normal ----
test_that("sf_normal handles invalid inputs correctly", {

  # Test that function errors with negative mean
  expect_error(
    sf_normal(mean = -1, sd = 1),
    "Normal distribution only accepts `mean>=0` and `sd>0"
  )

  # Test that function errors with non-positive sd
  expect_error(
    sf_normal(mean = 0, sd = 0),
    "Normal distribution only accepts `mean>=0` and `sd>0`"
  )
  expect_error(
    sf_normal(mean = 0, sd = -1),
    "Normal distribution only accepts `mean>=0` and `sd>0`"
  )

})

# Test for sf_uniform ----
test_that("sf_uniform handles invalid inputs correctly", {

  # Test that function errors with negative min
  expect_error(
    sf_uniform(min = -1, max = 5),
    "Uniform distribution only accepts 0<=min<max"
  )

  # Test that function errors if min is equal to max
  expect_error(
    sf_uniform(min = 5, max = 5),
    "Uniform distribution only accepts 0<=min<max"
  )

  # Test that function errors if min is greater than max
  expect_error(
    sf_uniform(min = 10, max = 5),
    "Uniform distribution only accepts 0<=min<max"
  )

})

# Test for sf_cauchy ----
test_that("sf_cauchy handles valid and invalid inputs correctly", {

  # Test the structure of the output for valid inputs
  result <- sf_cauchy(location = 1, scale = 2)
  expect_type(result, "list")
  expect_equal(result$location, 1)
  expect_equal(result$scale, 2)
  expect_equal(result$name, "cauchy")

  # Test that function errors with negative location
  expect_error(
    sf_cauchy(location = -1, scale = 5),
    "Cauchy distribution only accepts `location>=0` and `scale>=0` for median and median absolute deviation"
  )

  # Test that function errors with negative scale
  expect_error(
    sf_cauchy(location = 1, scale = -5),
    "Cauchy distribution only accepts `location>=0` and `scale>=0` for median and median absolute deviation"
  )

})

# Test build_stan_data ----
test_that("build_stan_data works with basic inputs", {

  stan_data <- build_stan_data(serosurvey)
  expect_equal(stan_data$n_observations, 5)
  expect_equal(stan_data$age_max, 25)
  expect_equal(stan_data$ages, 1:25)
  expect_equal(stan_data$n_seropositive, serosurvey$n_seropositive)
  expect_equal(stan_data$n_sample, serosurvey$n_sample)
})

test_that("build_stan_data sets right FoI index with constant model_type", {
  stan_data <- build_stan_data(serosurvey, model_type = "constant")
  expect_equal(stan_data$foi_index, rep(1, max(serosurvey$age_max)))
})

test_that("build_stan_data handles uniform prior for foi", {
  foi_prior <-
  stan_data <- build_stan_data(
    serosurvey,
    foi_prior = sf_uniform(min = 0, max = 10)
    )

  expect_equal(stan_data$foi_prior_index, 0) # assuming uniform is indexed as 0
  expect_equal(stan_data$foi_min, 0)
  expect_equal(stan_data$foi_max, 10)
})

test_that("build_stan_data handles normal prior for foi", {
  stan_data <- build_stan_data(
    serosurvey,
    model_type = "age",
    foi_prior = sf_normal(mean = 1, sd = 0.5),
    foi_index = get_foi_index(serosurvey, group_size = 5, model_type = "age")
    )

  expect_equal(stan_data$foi_prior_index, 1) # assuming normal is indexed as 1
  expect_equal(stan_data$foi_mean, 1)
  expect_equal(stan_data$foi_sd, 0.5)
})

test_that("build_stan_data handles cauchy prior for foi_sigma_rw", {
  stan_data <- build_stan_data(
    serosurvey,
    model_type = "age",
    foi_sigma_rw = sf_cauchy(location = 0, scale = 1)
    )

  expect_equal(stan_data$foi_sigma_rw_loc, 0)
  expect_equal(stan_data$foi_sigma_rw_sc, 1)
})

test_that("build_stan_data errors on missing seroreversion_prior if is_seroreversion is TRUE", {
  expect_error(
    build_stan_data(serosurvey, is_seroreversion = TRUE, seroreversion_prior = sf_none()),
    "seroreversion_prior not specified"
  )
})

test_that("build_stan_data handles seroreversion_prior with uniform and normal priors", {
  # Uniform prior
  stan_data <- build_stan_data(
    serosurvey,
    is_seroreversion = TRUE,
    seroreversion_prior = sf_uniform(min = 0, max = 10)
  )

  expect_equal(stan_data$seroreversion_prior_index, 0) # assuming uniform is indexed as 0
  expect_equal(stan_data$seroreversion_min, 0)
  expect_equal(stan_data$seroreversion_max, 10)

  # Normal prior
  stan_data <- build_stan_data(
    serosurvey,
    is_seroreversion = TRUE,
    seroreversion_prior = sf_normal(mean = 0, sd = 1e-4)
  )

  expect_equal(stan_data$seroreversion_prior_index, 1) # assuming uniform is indexed as 1
  expect_equal(stan_data$seroreversion_mean, 0)
  expect_equal(stan_data$seroreversion_sd, 1e-4)
})
