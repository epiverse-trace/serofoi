rm(list=ls())

library(dplyr)
library(serofoi)
library(tidyverse)
library(pracma)
library(cowplot)
library(Hmisc)

# TODO I would like to take this functions out of generate_sim_data in order to incorporate them to seroprevalence_data module
# For some reason when I take this functions out of generate_sim_data R gives me the Error: C stack usage  7970980 is too close to the limit.

# get_foi_probabilities <- function(exposure_matrix, foi){
#   probabilities <- purrr::map_dbl(1:nrow(exposure_matrix), ~1-exp(-dot(exposure_matrix[., ], foi)))
#   return(probabilities)
# }
#
# generate_sim_data <- function(foi, sample_size, exposure_matrix){
#   probabilities <- get_foi_probabilities(exposure_matrix, foi)
#   counts <- purrr::map_int(probabilities, ~rbinom(1, sample_size, .))
#   return(counts)
# }

generate_sim_data <- function(n, foi, name_example = 'fake', grouping = FALSE) {
  dat <- data.frame(birth_year = c(2000:2049)) %>%
    mutate(tsur = rep(2050, length(birth_year)),
           age_mean_f = 2050 - birth_year) %>%
    mutate(tsur = tsur)
  yexpo <- get_exposure_years(dat)
  yexpo <- yexpo[-length(yexpo)]
  RealYexpo <- (min(dat$birth_year):dat$tsur[1])[-1]
  exposure_matrix <- get_exposure_matrix(dat, yexpo)
  Nobs <- nrow(dat)

  get_foi_probabilities <- function(exposure_matrix, foi){
    probabilities <- purrr::map_dbl(1:nrow(exposure_matrix), ~1-exp(-dot(exposure_matrix[., ], foi)))
    return(probabilities)
  }
  generate_sim_data <- function(foi, sample_size, exposure_matrix){
    probabilities <- get_foi_probabilities(exposure_matrix, foi)
    counts <- purrr::map_int(probabilities, ~rbinom(1, sample_size, .))
    return(counts)
  }

  dat$counts <- generate_sim_data(foi, n, exposure_matrix)
  dat$total <- n
  dat$survey <- name_example
  dat$country <- 'None'

  dat <- dat %>% mutate(age_group = 'NA', age = age_mean_f) %>% arrange(age)
  dat$age_group[dat$age > 0 & dat$age < 5] <-   '01-04'
  dat$age_group[dat$age > 4 & dat$age < 10] <-  '05-09'
  dat$age_group[dat$age > 9 & dat$age < 15] <-  '10-14'
  dat$age_group[dat$age > 14 & dat$age < 20] <- '15-19'
  dat$age_group[dat$age > 19 & dat$age < 25] <- '20-24'
  dat$age_group[dat$age > 24 & dat$age < 30] <- '25-29'
  dat$age_group[dat$age > 29 & dat$age < 35] <- '30-34'
  dat$age_group[dat$age > 34 & dat$age < 40] <- '35-39'
  dat$age_group[dat$age > 39 & dat$age < 45] <- '40-44'
  dat$age_group[dat$age > 44 & dat$age < 51] <- '45-50'

  if (grouping == TRUE) {
    datg <- dat %>% group_by(age_group) %>%
      dplyr::summarise(total = sum(total), counts = sum(counts)) %>%
      mutate(tsur = dat$tsur[1], tsur = dat$tsur[1]) %>%
      mutate(age_min = as.numeric(substr(age_group, 1, 2)),
             age_max = as.numeric(substr(age_group, 4, 5))) %>%
      mutate(age_mean_f = floor((age_min + age_max)/2)) %>%
      mutate(age = age_mean_f, age_mean = age_mean_f,
              survey = name_example,
              country = 'None') %>%
      mutate(birth_year = tsur - age_mean_f)

    dat <- datg[, names(dat)] %>% arrange(age)
    conf <- data.frame(Hmisc::binconf(dat$counts, dat$total, method = "exact"))
    dat <- cbind(dat, conf) %>% dplyr::rename(
      p_obs = PointEst,
      p_obs_l = Lower,
      p_obs_u = Upper
    )
  }
  else{
    dat <- dat %>% mutate(p_obs = counts/total,
                          p_obs_l = 0,
                          p_obs_u = 0)
  }
  dat <- dat %>% mutate(sample_size = sum(total))
  return(dat)
}


no_transm <- 0.0000000001
small_outbreak <- 0.5
big_outbreak   <- 1.5

foiA <- rep(0.02, 50)
foiB <- c(rep(0.2, 25), rep(0.1, 10), rep(no_transm, 15)) # interruption
foiC <- seq(no_transm, 0.05, length.out = 50) # interruption
foiD <- c(rep(no_transm, 32), rep(big_outbreak, 3), rep(no_transm, 15)) # 1 epidemics
foiE <- c(rep(no_transm, 37), rep(small_outbreak, 3), rep(no_transm, 10)) # 1 epidemics
foiF <- c(rep(no_transm, 11), rep(big_outbreak, 2), rep(no_transm, 23), rep(small_outbreak, 3), rep(no_transm, 11)) # 1 epidemics

genetate_4_datasets <- function(foi, name_foi)  {
  foi1 <- generate_sim_data(n = 5, foi,  paste0(name_foi, '_n05_sing')) %>% mutate(simulation = name_foi)
  foi2 <- generate_sim_data(n = 5, foi,  paste0(name_foi, '_n05_group'), grouping = TRUE) %>% mutate(simulation = name_foi)
  foi3 <- generate_sim_data(n = 10, foi, paste0(name_foi, '_n10_sing')) %>% mutate(simulation = name_foi)
  foi4 <- generate_sim_data(n = 10, foi, paste0(name_foi, '_n10_group'), grouping = TRUE) %>% mutate(simulation = name_foi)
  datasets <- rbind(foi1, foi2, foi3, foi4) %>% mutate(antibody = 'IgG', test = 'simulated')
  return(datasets)
}

datA <- genetate_4_datasets(foiA, 'foiA')
datB <- genetate_4_datasets(foiB, 'foiB')
datC <- genetate_4_datasets(foiC, 'foiC')
datD <- genetate_4_datasets(foiD, 'foiD')
datE <- genetate_4_datasets(foiE, 'foiE')
datF <- genetate_4_datasets(foiF, 'foiF')

plot_sim_data <- function(data_sim, n_row=1){
  data_sim_plot <- ggplot(data = data_sim, aes(x = age, y = p_obs)) +
    geom_errorbar(aes(age, ymin = p_obs_l, ymax = p_obs_u), width = 0.1) +
    geom_point(aes(age, p_obs, size = data_sim$total), fill = "#7a0177", colour = "black", shape = 21) +
    facet_wrap(~survey, nrow = n_row) +
    coord_cartesian(xlim = c(0, 50), ylim = c(0,1)) +
    theme(legend.position = "none")
  return(data_sim_plot)
}

plot_datA <- plot_sim_data(data_sim = datA)
plot_datB <- plot_sim_data(data_sim = datB)
plot_datC <- plot_sim_data(data_sim = datC)
plot_datD <- plot_sim_data(data_sim = datD)
plot_datE <- plot_sim_data(data_sim = datE)
plot_datF <- plot_sim_data(data_sim = datF)

grid_plot <- plot_grid(plot_datA, plot_datB, plot_datC, plot_datD, plot_datE, plot_datF,
                       ncol = 1, nrow = 6, rel_heights = c(1, 1, 1, 1, 1, 1))

save_plot("tests/sim_data_test.png", grid_plot, nrow = 6, ncol = 4)
# dat0 <- dat0 %>% mutate(Antibody = 'IgG', Test = 'Fake')
# Here we assume the test used is IgG
# dat0 %>%
#   ggplot(aes(x = age, y = prev_obs)) +
#   geom_point() +
#   facet_wrap(~survey, nrow = 2)


# fake_data <- list(dat0 = dat0,
#                   foiA = foiA,
#                   foiB = foiB,
#                   foiC = foiC,
#                   foiD = foiD,
#                   foiE = foiE,
#                   foiF = foiF)

#(name_file <- paste0('data/fake_data_', round(as.numeric(Sys.time()), 0), '.RDS'))
#saveRDS(fake_data, name_file)
