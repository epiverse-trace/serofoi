% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_seromodel.R
\name{plot_rhats}
\alias{plot_rhats}
\title{Plot r-hats convergence criteria for the specified model}
\usage{
plot_rhats(seromodel, serosurvey, par_name = "foi_expanded", size_text = 11)
}
\arguments{
\item{seromodel}{stan_fit object obtained from sampling a model
with \link{fit_seromode}}

\item{serosurvey}{\describe{
\item{\code{survey_year}}{Year in which the survey took place
(only needed to plot time models)}
\item{\code{age_min}}{Floor value of the average between age_min and age_max}
\item{\code{age_max}}{The size of the sample}
\item{\code{n_sample}}{Number of samples for each age group}
\item{\code{n_seropositive}}{Number of positive samples for each age group}
}}

\item{par_name}{String specifying the parameter to be extracted
from \code{seromodel}}

\item{size_text}{Size of text for plotting (\code{base_size} in
\link[ggplot2:ggtheme]{ggplot2})}
}
\value{
ggplot object showing the r-hats of the model to be compared with the
convergence criteria (horizontal dashed line)
}
\description{
Plot r-hats convergence criteria for the specified model
}
\examples{
data(chagas2012)
seromodel <- fit_seromodel(
  serosurvey = chagas2012,
  model_type = "time",
  foi_index = data.frame(
    year = 1935:2011,
    foi_index = c(rep(1, 46), rep(2, 31))
  ),
  iter = 100,
  chains = 2
)
plot_rhats(seromodel, chagas2012)
}
