set.seed(123)
data("serodata")
data_test <- prepare_serodata(serodata = serodata)
model_constant <- run_seromodel(serodata = data_test,
                                foi_model = "constant",
                                n_iters = 800)

rhats_df <- get_table_rhats(seromodel_object = model_constant)

test_that("The structure of rhats dataframe is correct", {
  expect_s3_class(rhats_df, "data.frame")
  expect_equal(ncol(rhats_df), 2)
  expect_false(any(is.nan(rhats_df$rhat)))
  expect_equal(names(rhats_df), c("year", "rhat"))
})
