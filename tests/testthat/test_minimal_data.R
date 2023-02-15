# TODO Complete
test_that("minimal data", {
  library(devtools)
  library(dplyr)

  #----- Read Data
  data_test_path <- test_path(
    "extdata", "data.RDS"
  )
  data_test0 <- readRDS(data_test_path)  

  age_mean_f <- floor((data_test0$age_min + data_test0$age_max) / 2)
  age_mean_f == data_test0$age_mean_f

  sample_size <- sum(data_test0$total)
  sample_size == data_test0$sample_size

  bin_int <- Hmisc::binconf(data_test0$counts,
    data_test0$total,
    alpha = 0.05,
    method = "exact",
    return.df = TRUE
  ) %>%
    rename(prev_obs = PointEst, prev_obs_lower = Lower, prev_obs_upper = Upper)

bin_int$prev_obs == data_test0$prev_obs
bin_int$prev_obs_lower == data_test0$prev_obs_lower
bin_int$prev_obs_upper == data_test0$prev_obs_upper

# data_test <- data_test0 %>% select('survey',
#                                    'total',
#                                    'counts',
#                                    'age_min',
#                                    'age_max',
#                                    'year_init',
#                                    'year_end',
#                                    'tsur',
#                                    'country',
#                                    'test',
#                                    'antibody')
# saveRDS(data_test, file="data/data.RDS")

  data_test <- readRDS(test_path("extdata", "data.RDS"))
  data_test <- prepare_data(data_test, alpha = 0.05)

  model_0_object <- run_model(
    model_data = data_test,
    model_name = "constant_foi_bi",
    n_iters = 1000
  )
# TODO Complete test ###
})
