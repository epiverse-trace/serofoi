% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_serosurvey.R
\name{prob_seroprev_time_by_age}
\alias{prob_seroprev_time_by_age}
\title{Generate probabilities of seropositivity by age based on a time-varying FOI
model.}
\usage{
prob_seroprev_time_by_age(foi, seroreversion_rate)
}
\arguments{
\item{foi}{A dataframe containing the force of infection (FOI) values
for different years. It should have two columns: 'year' and 'foi'.}

\item{seroreversion_rate}{A non-negative numeric value representing the
rate of seroreversion.}
}
\value{
A dataframe with columns 'age' and 'seropositivity'.
}
\description{
This function calculates the probabilities of seropositivity by age based on
a time-varying FOI model.
It takes into account the FOI and the rate of seroreversion.
}
