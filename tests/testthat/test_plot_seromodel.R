library(testthat)

# Prueba para excepción en else
test_that("plot_seromodel imprime un mensaje de error y devuelve un objeto de trazado vacío cuando no se puede ajustar un modelo", {
  # Crea un objeto de seromodel con el ajuste como característica
  seromodel_object <- list(fit = "no_fit", model = "mi_modelo")
  # Ejecuta la función y comprueba que devuelve un objeto de trazado vacío
  expect_silent(plot_seromodel(seromodel_object))
  expect_equal(length(plot_seromodel(seromodel_object)$grobs), 5)
})
