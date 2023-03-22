test_that("simulated data", {
    # So far we are skipping tests on these platforms until
    # we find an efficient way to update rstan testthat snapshots on all of them
    skip_on_os(c("windows", "mac"))
    skip_on_ci()

    library(dplyr)
    library(serofoi)
    library(pracma)

    seed = 1234
    size_age_class = 5
    #----- Test foi case A
    sim_foi <- rep(0.02, 50)
    case_label <- "case_A_"
    max_lambda <- 0.035

    #----- Test foi case B
    # no_transm <- 0.0000000001
    # sim_foi <- c(rep(0.2, 25), rep(0.1, 10), rep(no_transm, 15))
    # case_label <- "case_B_"
    # max_lambda <- 0.3

    #----- Results paths
    data_path_no_grouped <- testthat::test_path(
      "extdata", paste0(case_label, "sim_data_no_grouped.csv")
    )

    data_path_grouped <- testthat::test_path(
      "extdata", paste0(case_label, "sim_data_grouped.csv")
    )

    #----- Data simulation
    sim_data <- generate_sim_data(foi = sim_foi,
                                  size_age_class = size_age_class,
                                  tsur = 2050,
                                  birth_year_min = 2000,
                                  survey_label = 'sim_foi',
                                  seed = seed)

    sim_data <- sim_data %>% mutate(age_min = age_mean_f, age_max = age_mean_f)
    write.csv(sim_data, data_path_no_grouped, row.names = FALSE)

    model_constant <- run_seroprev_model(seroprev_data = sim_data,
                                          seroprev_model_name = "constant_foi_bi")

    model_normal <- run_seroprev_model(seroprev_data = sim_data,
                                      seroprev_model_name = "continuous_foi_normal_bi")

    model_log <- run_seroprev_model(seroprev_data = sim_data,
                                    seroprev_model_name = "continuous_foi_normal_log")

    model_plot_constant <- plot_seroprev_model(model_constant, size_text = 6, , max_lambda = max_lambda)
    model_plot_normal <- plot_seroprev_model(model_normal, size_text = 6, , max_lambda = max_lambda)
    model_plot_log <- plot_seroprev_model(model_log, size_text = 6, max_lambda = max_lambda)
    model_plot_arrange <- plot_seroprev_models_grid(model_plot_constant,
                                                    model_plot_normal,
                                                    model_plot_log,
                                                    size_text = 6,
                                                    n_col = 3,
                                                    n_row = 1)
    vdiffr::expect_doppelganger(paste0(case_label, "no_group_models_", seed), model_plot_arrange)

    sim_foi_plot_constant <- plot_foi(model_constant, size_text = 10, max_lambda = max_lambda) +
      ggplot2::geom_point(data = data.frame(year = model_constant$exposure_years,
                                            foi = sim_foi),
                          ggplot2::aes(year, foi)) +
      ggplot2::ggtitle(paste0(case_label, "no_group_", seed))

    sim_foi_plot_normal <- plot_foi(model_normal, size_text = 10, max_lambda = max_lambda) +
      ggplot2::geom_point(data = data.frame(year = model_normal$exposure_years,
                                            foi = sim_foi),
                          ggplot2::aes(year, foi))
    sim_foi_plot_log <- plot_foi(model_log, size_text = 10, max_lambda = max_lambda) +
      ggplot2::geom_point(data = data.frame(year = model_log$exposure_years,
                                            foi = sim_foi),
                          ggplot2::aes(year, foi))
    sim_foi_plot_arrange <- plot_seroprev_models_grid(sim_foi_plot_constant,
                                                      sim_foi_plot_normal,
                                                      sim_foi_plot_log,
                                                      n_col = 1, n_row = 3)

    vdiffr::expect_doppelganger(paste0(case_label, "no_group_foi_", seed), sim_foi_plot_arrange)


    sim_data_grouped <- group_sim_data(sim_data = sim_data,
                                      foi = sim_foi,
                                      size_age_class = size_age_class,
                                      tsur = 2050,
                                      birth_year_min = 2000,
                                      survey_label = 'sim_foi')
    write.csv(sim_data_grouped, data_path_grouped, row.names = FALSE)

    model_grouped_constant <- run_seroprev_model(seroprev_data = sim_data_grouped,
                                          seroprev_model_name = "constant_foi_bi")
    model_grouped_normal <- run_seroprev_model(seroprev_data = sim_data_grouped,
                                              seroprev_model_name = "continuous_foi_normal_bi")
    model_grouped_log <- run_seroprev_model(seroprev_data = sim_data_grouped,
                                            seroprev_model_name = "continuous_foi_normal_log")

    model_grouped_constant_plot <- plot_seroprev_model(model_grouped_constant, size_text = 6, max_lambda = max_lambda)
    model_grouped_normal_plot <- plot_seroprev_model(model_grouped_normal, size_text = 6, max_lambda = max_lambda)
    model_grouped_log_plot <- plot_seroprev_model(model_grouped_log, size_text = 6, max_lambda = max_lambda)
    model_grouped_plot_arrange <- plot_seroprev_models_grid(model_grouped_constant_plot,
                                                            model_grouped_normal_plot,
                                                            model_grouped_log_plot,
                                                            n_col = 3, n_row = 1)
    vdiffr::expect_doppelganger(paste0(case_label, "group_models_", seed), model_grouped_plot_arrange)
    rm(list = ls(pat = "*_plot"))

    sim_foi_grouped_constant_plot <- plot_foi(model_grouped_constant, size_text = 10, max_lambda = max_lambda) +
      ggplot2::geom_point(data = data.frame(year = model_constant$exposure_years,
                                            foi = sim_foi),
                          ggplot2::aes(year, foi)) +
      ggplot2::ggtitle(paste0(case_label, "group_", seed))
    sim_foi_grouped_normal_plot <- plot_foi(model_grouped_normal, size_text = 10, max_lambda = max_lambda) +
      ggplot2::geom_point(data = data.frame(year = model_normal$exposure_years,
                                            foi = sim_foi),
                          ggplot2::aes(year, foi))
    sim_foi_grouped_log_plot <- plot_foi(model_grouped_log, size_text = 10, max_lambda = max_lambda) +
      ggplot2::geom_point(data = data.frame(year = model_log$exposure_years,
                                            foi = sim_foi),
                          ggplot2::aes(year, foi))
    sim_foi_grouped_plot_arrange <- plot_seroprev_models_grid(sim_foi_grouped_constant_plot,
                                                              sim_foi_grouped_normal_plot,
                                                              sim_foi_grouped_log_plot,
                                                              n_col = 1, n_row = 3)
    vdiffr::expect_doppelganger(paste0(case_label, "group_foi_", seed), sim_foi_grouped_plot_arrange)
})
# TODO: solve error 'address 0x18, cause 'memory not mapped' when testing.
