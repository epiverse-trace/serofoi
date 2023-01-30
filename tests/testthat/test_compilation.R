test_that("hello works", {
  print("aaa")
  testthat::expect_equal(
    "hello",
    "hello"
  )
})
