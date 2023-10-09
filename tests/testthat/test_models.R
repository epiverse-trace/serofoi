library(testthat)

#----- Test for "constant" model
test_that("Test model - constant", {
  set.seed(1234) # For reproducibility
  library(devtools)
  library(vdiffr)
  
  data(simdata_constant)
  serodata <- prepare_serodata(simdata_constant)
   
   #----- Test for get_cohort_ages
  cohort_ages <- get_cohort_ages(serodata = serodata)
  expect_equal(nrow(cohort_ages), max(unique(serodata$tsur)) - min(serodata$birth_year))
  
  #----- Test for the constant model
  foi_sim <- rep(0.02, nrow(cohort_ages)+1) 
  foi_model <- "constant"
  seromodel_object <- run_seromodel(serodata = serodata,
                                    foi_model = foi_model,
                                    n_iters = 1000,
                                    print_summary = FALSE)

  #----- Results visualisation
  size_text <- 6
  max_lambda <- 0.025
  seromodel_plot <- plot_seromodel(seromodel_object = seromodel_object,
                                    serodata = serodata,
                                    size_text = size_text,
                                    max_lambda = max_lambda,
                                    foi_sim = foi_sim)
  
  vdiffr::expect_doppelganger(paste0("plot_", foi_model), seromodel_plot)
})

#----- Test for "tv_normal" model
test_that("Test model - tv_normal", {
  set.seed(1234) # For reproducibility
  library(devtools)
  library(vdiffr)
  
  data(chagas2012)
  serodata <- prepare_serodata(chagas2012)

   #----- Test for get_cohort_ages
  cohort_ages <- get_cohort_ages(serodata = serodata)
  expect_equal(nrow(cohort_ages), max(unique(serodata$tsur)) - min(serodata$birth_year))
  
  #----- Test for the constant model
  foi_model <- "tv_normal"
  seromodel_object <- run_seromodel(serodata = serodata,
                                    foi_model = foi_model,
                                    n_iters = 1000,
                                    print_summary = TRUE)
  
  #----- Test for get_prev_expanded
  data_constant_path <- testthat::test_path("extdata", paste0("prev_expanded_", foi_model, ".RDS"))
  foi <- rstan::extract(seromodel_object, "foi", inc_warmup = FALSE)[[1]]

  prev_expanded <- get_prev_expanded(foi, serodata = serodata)
  prev_expanded_constant <- readRDS(data_constant_path)
  testthat::expect_equal(prev_expanded, prev_expanded_constant, tolerance = TRUE)

  #----- Results visualisation
  size_text <- 6
  max_lambda <- 0.025
  seromodel_plot <- plot_seromodel(seromodel_object = seromodel_object,
                                    serodata = serodata,
                                    size_text = size_text,
                                    max_lambda = max_lambda)
  
  vdiffr::expect_doppelganger(paste0("plot_", foi_model), seromodel_plot)
})

#----- Test for "tv_normal_log" model
test_that("Test model - tv_normal_log", {
  set.seed(1234) # For reproducibility
  library(devtools)
  library(vdiffr)
  
  data(veev2012)
  serodata <- prepare_serodata(veev2012)

   #----- Test for get_cohort_ages
  cohort_ages <- get_cohort_ages(serodata = serodata)
  expect_equal(nrow(cohort_ages), max(unique(serodata$tsur)) - min(serodata$birth_year))
  
  #----- Test for the constant model
  foi_model <- "tv_normal_log"
  seromodel_object <- run_seromodel(serodata = serodata,
                                    foi_model = foi_model,
                                    n_iters = 1000,
                                    print_summary = TRUE)

  #----- Results visualisation
  size_text <- 6
  max_lambda <- 0.6
  seromodel_plot <- plot_seromodel(seromodel_object = seromodel_object,
                                    serodata = serodata,
                                    size_text = size_text,
                                    max_lambda = max_lambda)
  
  vdiffr::expect_doppelganger(paste0("plot_", foi_model), seromodel_plot)
})
