test_that("simulated data", {

    library(dplyr)
    library(serofoi)

    seed <- 1234
    sample_size_by_age <- 5
    #----- Test for constant FoI
    # foi_model <- "constant"
    # foi_sim <- rep(0.02, 50)
    # case_label <- "constant_foi_"
    # max_lambda <- 0.035

    #----- Test for stepwise decreasing FoI
    foi_model = "tv_normal"
    no_transm <- 0.0000000001
    foi_sim <- c(rep(0.2, 25), rep(0.1, 10), rep(no_transm, 15))
    case_label <- "sw_dec_foi_"
    max_lambda <- 0.3

    #----- Test function generate_sim_data
    sim_data <- generate_sim_data(foi = foi_sim,
                                  sample_size_by_age = sample_size_by_age,
                                  tsur = 2050,
                                  birth_year_min = 2000,
                                  survey_label = 'foi_sim',
                                  seed = seed)

    expect_s3_class(sim_data, "data.frame")

    #----- Test function group_sim_data
    sim_data <- sim_data %>% mutate(age_min = age_mean_f, age_max = age_mean_f)
    sim_data_grouped <- group_sim_data(sim_data = sim_data)
    expect_s3_class(sim_data_grouped, "data.frame")
    expect_s3_class(sim_data_grouped$age_group, "factor")
})
