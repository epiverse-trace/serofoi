% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seroprevalence_data.R
\name{generate_sim_data}
\alias{generate_sim_data}
\title{Generate simulated serosurvey according to the specified FoI}
\usage{
generate_sim_data(
  sim_data,
  foi,
  sample_size_by_age,
  mu = 0,
  model_type = "time-varying",
  survey_label = "sim_data",
  seed = 1234
)
}
\arguments{
\item{sim_data}{A dataframe object containing the following columns:
\describe{
\item{\code{age}}{Age group markers}
\item{\code{tsur}}{Year of the survey}
}}

\item{foi}{Numeric atomic vector corresponding to the desired
time-varying or age-varying Force-of-Infection to simulate from}

\item{sample_size_by_age}{Integer indicating the sample size by age group.
This corresponds to the number of trials \code{size} in \link[stats:Binomial]{rbinom}.}

\item{mu}{Seroreversion rate}

\item{model_type}{String specifying the type of model to be used.
Current valid options are 'time-varying' and 'age-varying'}

\item{survey_label}{Label for the resulting simulated serosurvey.}

\item{seed}{Seed for random number generation.}
}
\value{
Dataframe containing the simulated serosurvey.
}
\description{
Generate simulated serosurvey according to the specified FoI
}
\examples{
n_years <- 50
sim_data <- data.frame(
  age = seq(1,n_years),
  tsur = 2050
)
foi <- rep(0.02, n_years)
sample_size_by_age <- as.integer(runif(n = n_years, 5, 10))
sim_data <- generate_sim_data(
  sim_data = sim_data,
  foi = foi,
  sample_size_by_age = sample_size_by_age,
  survey_label = "sim_constant_foi"
)
}
