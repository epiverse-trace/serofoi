% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualisation.R
\name{plot_info_table}
\alias{plot_info_table}
\title{Generate plot summarizing a given table}
\usage{
plot_info_table(info_table, size_text)
}
\arguments{
\item{info_table}{Table with the information to be summarised}

\item{size_text}{Text size of the graph returned by the function}
}
\value{
ggplot object summarizing the information in \code{info_table}
}
\description{
Generate plot summarizing a given table
}
\examples{
serodata <- prepare_serodata(chagas2012)
seromodel_object <- fit_seromodel(
  serodata = serodata,
  foi_model = "constant",
  iter = 1000
)
seromodel_summary <- extract_seromodel_summary(
  seromodel_object = seromodel_object,
  serodata = serodata
)
info_table <- t(seromodel_summary)
plot_info_table(info_table, size_text = 15)
}
