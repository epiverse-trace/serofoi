library(testthat)

# Definir el contexto de la prueba
context("Test de la función plot_foi")

# Crear una prueba para verificar el bloque else
test_that("La función plot_foi debe generar un gráfico vacío cuando el modelo no se ejecuta", {

  # Crear un objeto vacío que simule la salida del modelo que falló
  empty_model <- list(fit = "failure")

  # Ejecutar la función plot_foi con el modelo vacío
  plot <- plot_foi(empty_model)

  # Verificar si la salida es un gráfico vacío
  expect_identical(ggplot2::ggplot(), plot)
})
