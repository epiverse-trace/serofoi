rm(list=ls())
# remotes::install_github("TRACE-LAC/serofoi", ref = "dev")
library(dplyr)
library(serofoi)
library(tidyverse)
library(pracma)

# Generate exposure matrix

generate_fake_data <- function(n, foi, name_example = 'fake', grouping = FALSE) {
  # dat <- readRDS('data/chik_serop_IND 1965(S095).RDS')
  dat <- data.frame(birth_year=c(2000:2049)) %>%
    mutate(tSur=rep(2050, length(birth_year)),
           age_mean_f=2050-birth_year) %>%
    mutate(tsur=tSur)
  yexpo <- get_exposure_years(dat)
  yexpo <- yexpo[-length(yexpo)]
  RealYexpo <- (min(dat$birth_year):dat$tsur[1])[-1]
  ExposureMatrix <- get_exposure_matrix(dat, yexpo)
  Nobs <- nrow(dat)

  fGenerateProbs <- function(ExposureMatrix, foi){
    p <- map_dbl(1:nrow(ExposureMatrix), ~1-exp(-dot(ExposureMatrix[., ], foi)))
    return(p)
  }

  fGenerateFakeData <- function(foi, sample_size, ExposureMatrix){
    p <- fGenerateProbs(ExposureMatrix, foi)
    counts <- map_int(p, ~rbinom(1, sample_size, .))
    return(counts)
  }


  dat$counts <- fGenerateFakeData(foi, n, ExposureMatrix)
  dat$total <- n
  dat$survey <- name_example
  dat$country <- 'None'

  dat <- dat %>% mutate(age_group = 'NA', age = age_mean_f)
  dat$age_group[dat$age >0 & dat$age <5] <-   '01-04'
  dat$age_group[dat$age >4 & dat$age <10] <-  '05-09'
  dat$age_group[dat$age >9 & dat$age <15] <-  '10-14'
  dat$age_group[dat$age >14 & dat$age <20] <- '15-19'
  dat$age_group[dat$age >19 & dat$age <25] <- '20-24'
  dat$age_group[dat$age >24 & dat$age <30] <- '25-29'
  dat$age_group[dat$age >29 & dat$age <35] <- '30-34'
  dat$age_group[dat$age >34 & dat$age <40] <- '35-39'
  dat$age_group[dat$age >39 & dat$age <45] <- '40-44'
  dat$age_group[dat$age >44 & dat$age <51] <- '45-50'

  dat <- dat%>% arrange(age)

  if (grouping == TRUE) {

    datg <- dat %>% group_by(age_group) %>%
      dplyr::summarise(total = sum(total), counts = sum(counts)) %>%
      mutate(tsur = dat$tsur[1], tSur = dat$tSur[1]) %>%
      mutate(age_min = as.numeric(substr(age_group, 1, 2)),
             age_max = as.numeric(substr(age_group, 4, 5))) %>%
      mutate(age_mean_f = floor((age_min + age_max)/2)) %>%
      mutate (age = age_mean_f, age_mean = age_mean_f,
              survey = name_example,
              country = 'None') %>%
      mutate(birth_year = tsur - age_mean_f)

    dat <- datg[, names(dat)]

  }

  dat <- dat%>% arrange(age)
  conf_int <- data.frame(Hmisc::binconf(x=dat$counts,n = dat$total, alpha = 0.05, method ='exact'))
  dat <- cbind(dat, conf_int)
  dat <- dat %>% rename(prev_obs_lower = Lower,
                        prev_obs = PointEst,
                        prev_obs_upper = Upper)


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

genetate_4_datasets <- function (foi, name_foi)  {
  datfoiA1 <- generate_fake_data(n = 5, foi,  paste0(name_foi, '_n05_sing')) %>% mutate(simulation = name_foi)
  datfoiA2 <- generate_fake_data(n = 5, foi,  paste0(name_foi, '_n05_group'), grouping = TRUE) %>% mutate(simulation = name_foi)
  datfoiA3 <- generate_fake_data(n = 10, foi, paste0(name_foi, '_n10_sing')) %>% mutate(simulation = name_foi)
  datfoiA4 <- generate_fake_data(n = 10, foi, paste0(name_foi, '_n10_group'), grouping = TRUE)%>% mutate(simulation = name_foi)
  datasets <- list(foi1 = datfoiA1,
                   foi2 = datfoiA2,
                   foi3 = datfoiA3,
                   foi4 = datfoiA4)
  return(datasets)
}

foiA_ds <- genetate_4_datasets(foiA, 'foiA')
foiB_ds <- genetate_4_datasets(foiB, 'foiB')
foiC_ds <- genetate_4_datasets(foiC, 'foiC')
foiD_ds <- genetate_4_datasets(foiD, 'foiD')
foiE_ds <- genetate_4_datasets(foiE, 'foiE')
foiF_ds <- genetate_4_datasets(foiF, 'foiF')

dat0 <- rbind (foiA_ds$foi1, foiA_ds$foi2, foiA_ds$foi3, foiA_ds$foi4,
               foiB_ds$foi1, foiB_ds$foi2, foiB_ds$foi3, foiB_ds$foi4,
               foiC_ds$foi1, foiC_ds$foi2, foiC_ds$foi3, foiC_ds$foi4,
               foiD_ds$foi1, foiD_ds$foi2, foiD_ds$foi3, foiD_ds$foi4,
               foiE_ds$foi1, foiE_ds$foi2, foiE_ds$foi3, foiE_ds$foi4,
               foiF_ds$foi1, foiF_ds$foi2, foiF_ds$foi3, foiF_ds$foi4)

dat0 <- dat0 %>% mutate(Antibody = 'IgG', Test = 'Fake')
# Here we assume the test used is IgG
dat0 %>%
  ggplot(aes(x=age, y=prev_obs)) +
  geom_point() +
  facet_wrap(~survey, nrow = 2)


fake_data <- list  (dat0 = dat0,
                    foiA = foiA,
                    foiB = foiB,
                    foiC = foiC,
                    foiD = foiD,
                    foiE = foiE,
                    foiF = foiF)


(name_file <- paste0('data/fake_data_', round(as.numeric(Sys.time()), 0), '.RDS'))
saveRDS(fake_data, name_file)
