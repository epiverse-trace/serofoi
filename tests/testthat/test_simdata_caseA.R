test_that("simulated data", {
    # So far we are skipping tests on these platforms until
    # we find an efficient way to update rstan testthat snapshots on all of them
    skip_on_os(c("windows", "mac"))
    skip_on_ci()

    library(dplyr)
    library(serofoi)
    library(pracma)

#----- Test foi case A
    size_age_class_A = 5
    foi_A <- rep(0.02, 50)

    sim_data_A <- generate_sim_data(foi = foi_A,
                                    size_age_class = size_age_class_A,
                                    tsur = 2050,
                                    birth_year_min = 2000,
                                    survey_label = 'foi_A')

    sim_data_A <- sim_data_A %>% mutate(age_min = age_mean_f, age_max = age_mean_f)

    model_A_constant <- run_seroprev_model(seroprev_data = sim_data_A,
                                          seroprev_model_name = "constant_foi_bi")

    model_A_normal <- run_seroprev_model(seroprev_data = sim_data_A,
                                         seroprev_model_name = "continuous_foi_normal_bi")

    model_A_log <- run_seroprev_model(seroprev_data = sim_data_A,
                                      seroprev_model_name = "continuous_foi_normal_log")

    model_A_plot_constant <- plot_seroprev_model(model_A_constant, size_text = 6, , max_lambda = 0.05)
    model_A_plot_normal <- plot_seroprev_model(model_A_normal, size_text = 6, , max_lambda = 0.05)
    model_A_plot_log <- plot_seroprev_model(model_A_log, size_text = 6, max_lambda = 0.05)
    model_A_plot_arrange <- plot_seroprev_models_grid(model_A_plot_constant,
                                                      model_A_plot_normal,
                                                      model_A_plot_log,
                                                      size_text = 6,
                                                      n_col = 3,
                                                      n_row = 1)
    vdiffr::expect_doppelganger("plot_A_no_group", model_A_plot_arrange)

    foi_A_plot_constant <- plot_foi(model_A_constant, size_text = 10) +
      ggplot2::geom_point(data = data.frame(year = model_A_constant$exposure_years,
                                            foi = foi_A),
                          ggplot2::aes(year, foi))
    foi_A_plot_normal <- plot_foi(model_A_normal, size_text = 10) +
      ggplot2::geom_point(data = data.frame(year = model_A_normal$exposure_years,
                                            foi = foi_A),
                          ggplot2::aes(year, foi))
    foi_A_plot_log <- plot_foi(model_A_log, size_text = 10) +
      ggplot2::geom_point(data = data.frame(year = model_A_log$exposure_years,
                                            foi = foi_A),
                          ggplot2::aes(year, foi))
    foi_A_plot_arrange <- plot_seroprev_models_grid(foi_A_plot_constant,
                                                    foi_A_plot_normal,
                                                    foi_A_plot_log,
                                                    n_col = 1, n_row = 3)
    vdiffr::expect_doppelganger("plot_foi_A_no_group", foi_A_plot_arrange)


    sim_data_A_grouped <- generate_sim_data_grouped(sim_data = sim_data_A,
                                                    foi = foi_A,
                                                    size_age_class = size_age_class_A,
                                                    tsur = 2050,
                                                    birth_year_min = 2000,
                                                    survey_label = 'foi_A')

    model_A_grouped_constant <- run_seroprev_model(seroprev_data = sim_data_A_grouped,
                                          seroprev_model_name = "constant_foi_bi")
    model_A_grouped_normal <- run_seroprev_model(seroprev_data = sim_data_A_grouped,
                                                   seroprev_model_name = "continuous_foi_normal_bi")
    model_A_grouped_log <- run_seroprev_model(seroprev_data = sim_data_A_grouped,
                                                   seroprev_model_name = "continuous_foi_normal_log")

    model_A_grouped_constant_plot <- plot_seroprev_model(model_A_grouped_constant, size_text = 6)
    model_A_grouped_normal_plot <- plot_seroprev_model(model_A_grouped_normal, size_text = 6)
    model_A_grouped_log_plot <- plot_seroprev_model(model_A_grouped_log, size_text = 6)
    model_A_grouped_plot_arrange <- plot_seroprev_models_grid(model_A_grouped_constant_plot,
                                                              model_A_grouped_normal_plot,
                                                              model_A_grouped_log_plot,
                                                              n_col = 3, n_row = 1)
    vdiffr::expect_doppelganger("plot_A_group", model_A_grouped_plot_arrange)
    rm(list = ls(pat ="*_plot"))

    foi_A_grouped_constant_plot <- plot_foi(model_A_grouped_constant, size_text = 10) +
      ggplot2::geom_point(data = data.frame(year = model_A_constant$exposure_years,
                                            foi = foi_A),
                          ggplot2::aes(year, foi))
    foi_A_grouped_normal_plot <- plot_foi(model_A_grouped_normal, size_text = 10) +
      ggplot2::geom_point(data = data.frame(year = model_A_normal$exposure_years,
                                            foi = foi_A),
                          ggplot2::aes(year, foi))
    foi_A_grouped_log_plot <- plot_foi(model_A_grouped_log, size_text = 10) +
      ggplot2::geom_point(data = data.frame(year = model_A_log$exposure_years,
                                            foi = foi_A),
                          ggplot2::aes(year, foi))
    foi_A_grouped_plot_arrange <- plot_seroprev_models_grid(foi_A_grouped_constant_plot,
                                                            foi_A_grouped_normal_plot,
                                                            foi_A_grouped_log_plot,
                                                            n_col = 1, n_row = 3)
    vdiffr::expect_doppelganger("plot_foi_A_group", foi_A_grouped_plot_arrange)
})
# TODO: solve error 'address 0x18, cause 'memory not mapped' when testing. 
