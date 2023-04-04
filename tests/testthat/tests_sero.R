library(testthat)

# Creamos una función auxiliar que devuelve un objeto seromodel sin ajustar
create_dummy_seromodel <- function() {
  # Objeto vacío
  seromodel <- list()
  seromodel$fit <- "dummy"
  seromodel$serodata <- data.frame(age = c(0, 10, 20, 30, 40, 50, 60),
                                   p_obs_bin = c(0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99),
                                   bin_size = c(5, 10, 20, 30, 40, 50, 60))
  return(seromodel)
}

# Pruebas unitarias para plot_seroprev_fitted()
test_that("plot_seroprev_fitted() devuelve una gráfica vacía para un objeto seromodel sin ajustar", {
  # Creamos un objeto seromodel sin ajustar
  seromodel <- create_dummy_seromodel()
  # Llamamos a la función plot_seroprev_fitted()
  plot <- plot_seroprev_fitted(seromodel)
  # Verificamos que el objeto plot sea una ggplot vacía
  expect_true(class(plot) == "ggplot")
  expect_true(length(plot$layers) == 0)
})

test_that("plot_seroprev_fitted() devuelve un mensaje de error y una gráfica vacía cuando seromodel$fit es un string", {
  # Creamos un objeto seromodel con seromodel$fit como string
  seromodel <- create_dummy_seromodel()
  seromodel$fit <- "dummy"
  # Llamamos a la función plot_seroprev_fitted()
  plot <- plot_seroprev_fitted(seromodel)
  # Verificamos que el mensaje de error se imprima correctamente
  expect_output(plot_seroprev_fitted(seromodel), "model did not run")
  # Verificamos que el objeto plot sea una ggplot vacía
  expect_true(class(plot) == "ggplot")
  expect_true(length(plot$layers) == 0)
})

