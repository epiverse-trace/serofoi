set.seed(123)
data("serodata")
serodata <- prepare_serodata(serodata)
foi_constant <- run_seromodel(serodata, foi_model = "constant")

test_that("run_seromodel output is in right format", {
  # Load the already prepared haiti dataset
  data_path <- test_path("testdata", "haiti_ssa_sample.RDS")
  haiti_serodata <- readRDS(data_path)

  # Pareto k diagnostic values error reproduction
  model_fit <- run_seromodel(haiti_serodata, foi_model = "tv_normal",
                             n_iters = 1500, print_summary = FALSE)
  foi <- rstan::extract(model_fit$fit, "foi", inc_warmup = FALSE)[[1]]
  age_max <- max(haiti_serodata$age_mean_f)
  prev_expanded <- get_prev_expanded(foi, serodata = haiti_serodata)

  # Test
  expect_length(prev_expanded$age, n = age_max)
})

