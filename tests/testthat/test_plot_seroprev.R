library(testthat)

# Define la función de prueba
test_that("plot_seroprev_fitted() maneja las excepciones correctamente", {
  #skip_on_os(c("windows", "mac"))
  skip_on_ci()
  source("testing_utils.R")
  set.seed(1234) # For reproducibility

  library(devtools)
  library(dplyr)
  library(vdiffr)

  #----- Read and prepare data
  data("serodata")
  data_test <- serodata %>% prepare_serodata(alpha = 0.05)

  #----- Plot raw data
  data_test_plot <- plot_seroprev(data_test, size_text = 15)
  vdiffr::expect_doppelganger("serodata_plot", data_test_plot)


  # Prueba para cuando se proporciona un objeto de modelo inválido
  expect_warning(plot_seroprev_fitted(list(), size_text = 6), "model did not run")

  # Prueba para cuando se proporciona un objeto de modelo sin ajustar
  seromodel_object <- list(fit = "unfitted", serodata = data.frame())
  expect_warning(plot_seroprev_fitted(seromodel_object, size_text = 6), "model did not run")

  # Prueba para cuando se proporciona un objeto de modelo ajustado pero no se encuentra la clase 'stanfit'
  seromodel_object <- list(fit = "fitted", serodata = data.frame())
  expect_warning(plot_seroprev_fitted(seromodel_object, size_text = 6), "model did not run")

  # Prueba para cuando se proporciona un objeto de modelo ajustado con una clase 'stanfit'
  fit <- list(sim = list(samples = list(foi = matrix(rnorm(100), ncol = 1))))
  seromodel_object <- list(fit = fit, serodata = data.frame())
  expect_silent(plot_seroprev_fitted(seromodel_object, size_text = 6))

 })

