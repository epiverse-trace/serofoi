test_that("simulated data", {
    # So far we are skipping tests on these platforms until
    # we find an efficient way to update rstan testthat snapshots on all of them
    skip_on_os(c("windows", "mac"))
    skip_on_ci()

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

    sim_data <- sim_data %>% mutate(age_min = age_mean_f, age_max = age_mean_f)
    model_object <- run_seromodel(serodata = sim_data,
                                  foi_model = foi_model)

    model_plot <- plot_seromodel(model_object, size_text = 6, , max_lambda = max_lambda)
    vdiffr::expect_doppelganger(paste0(case_label, foi_model, "_no_group"), model_plot, foi_plot = foi_plot)

    foi_plot <- plot_foi(model_object, size_text = 10, max_lambda = max_lambda, foi_sim = foi_sim) +
      ggplot2::ggtitle(paste0(case_label, "no_group"))


    vdiffr::expect_doppelganger(paste0(case_label, foi_model, "_no_group_foi"), foi_plot)

    #----- Test function sim_data_grouped
    sim_data_grouped <- group_sim_data(sim_data = sim_data)
    model_object_grouped <- run_seromodel(serodata = sim_data_grouped,
                                          foi_model = foi_model)

    model_grouped_plot <- plot_seromodel(model_object_grouped, size_text = 6, max_lambda = max_lambda)

    foi_grouped_plot <- plot_foi(model_object_grouped, size_text = 10, max_lambda = max_lambda, foi_sim = foi_sim) +
      ggplot2::ggtitle(paste0(case_label, foi_model, "_group_"))

    vdiffr::expect_doppelganger(paste0(case_label, "group_foi"), foi_grouped_plot)
})
