% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_stan_data.R
\name{sf_normal}
\alias{sf_normal}
\title{Sets normal distribution parameters for sampling}
\usage{
sf_normal(mean = 0, sd = 1)
}
\arguments{
\item{mean}{Mean of the normal distribution}

\item{sd}{Standard deviation of the normal distribution}
}
\value{
List with specified statistics and name of the model
}
\description{
Sets normal distribution parameters for sampling
}
\examples{
my_prior <- sf_normal()
}
