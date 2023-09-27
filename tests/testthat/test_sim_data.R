test_that("simulated data", {

    library(dplyr)
    library(serofoi)

    seed <- 1234
    n_iters <- 1000
    sample_size_by_age <- 5
    tsur <- 2050
    birth_year_min <- 2000

    #----- Test for constant FoI
    # case_label <- "constant_foi_"
    # foi_model <- "constant"
    # foi_sim <- rep(0.02, tsur - birth_year_min)
    # max_lambda <- 0.035

    #----- Test for stepwise decreasing FoI
    case_label <- "smth_dec_foi_" #Smooth-decendent FoI
    foi_model = "tv_normal"
    foi_max = 0.2
    stretch = 0.15
    x <- 1:(tsur - birth_year_min)
    foi_sim <- (-foi_max * (atan(stretch * (x - 25))) / (0.5 * pi) + foi_max) / 2
    max_lambda <- 0.3

    #----- Test function generate_sim_data
    sim_data <- generate_sim_data(foi = foi_sim,
                                  sample_size_by_age = sample_size_by_age,
                                  tsur = tsur,
                                  birth_year_min = birth_year_min,
                                  survey_label = 'foi_sim',
                                  seed = seed)
    # Check sim_data structure
    expect_s3_class(sim_data, "data.frame")
    expect_length(sim_data$birth_year, tsur - birth_year_min)

    # Check consistency between sim_foi and the fitted foi
    sim_seromodel <- run_seromodel(sim_data, foi_model = foi_model, n_iters = n_iters)
    foi <- rstan::extract(sim_seromodel$seromodel_fit, "foi", inc_warmup = FALSE)[[1]]
    foi_lower <- apply(foi, 2, function(x) quantile(x, 0.05))
    foi_upper <- apply(foi, 2, function(x) quantile(x, 0.95))

    expect_true(all((foi_sim >= foi_lower) & (foi_sim <= foi_upper)))

    #----- Test function group_sim_data
    sim_data <- sim_data %>% mutate(age_min = age_mean_f, age_max = age_mean_f)
    sim_data_grouped <- group_sim_data(sim_data = sim_data)
    expect_s3_class(sim_data_grouped, "data.frame")
    expect_s3_class(sim_data_grouped$age_group, "factor")
})
