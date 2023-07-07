test_that("individual models", {
  set.seed(1234) # For reproducibility

  library(devtools)
  library(vdiffr)

  #----- Read and prepare data
  data(chagas2012)
  serodata <- prepare_serodata(chagas2012, alpha = 0.05)

  data_constant_path <- test_path("testdata", "prev_expanded_constant.RDS")
  data_tv_normal_path <- test_path("testdata", "prev_expanded_tv_normal.RDS")
  data_tv_normal_log_path <- test_path("testdata", "prev_expanded_tv_normal_log.RDS")

  prev_expanded_tv_normal_log <- readRDS(data_constant_path)

  #----- Test for the constant model

  model_name <- "constant"
  model_object <- run_seromodel(serodata = serodata,
                                foi_model = model_name,
                                n_iters = 1000,
                                print_summary = FALSE)

  foi <- rstan::extract(model_object$seromodel_fit, "foi", inc_warmup = FALSE)[[1]]
  prev_expanded <- get_prev_expanded(foi, serodata = model_object$serodata)
  prev_expanded_constant <- readRDS(data_constant_path)

  testthat::expect_equal(prev_expanded, prev_expanded_constant, tolerance = TRUE)

  #----- Test for the tv_normal model

  model_name <- "tv_normal"
  model_object <- run_seromodel(serodata = serodata,
                                foi_model = model_name,
                                n_iters = 1000)

  foi <- rstan::extract(model_object$seromodel_fit, "foi", inc_warmup = FALSE)[[1]]
  prev_expanded <- get_prev_expanded(foi, serodata = model_object$serodata)
  prev_expanded_tv_normal <- readRDS(data_tv_normal_path)
  testthat::expect_equal(prev_expanded, prev_expanded_tv_normal, tolerance = TRUE)

  #----- Test for the tv_normal_log model

  model_name <- "tv_normal_log"
  model_object <- run_seromodel(serodata = serodata,
                                foi_model = model_name,
                                n_iters = 1000)

  foi <- rstan::extract(model_object$seromodel_fit, "foi", inc_warmup = FALSE)[[1]]
  prev_expanded <- get_prev_expanded(foi, serodata = model_object$serodata)
  prev_expanded_tv_normal <- readRDS(data_tv_normal_path)
  testthat::expect_equal(prev_expanded, prev_expanded_tv_normal_log, tolerance = TRUE)

})
