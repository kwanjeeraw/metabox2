% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_stat_summary}
\alias{get_stat_summary}
\title{Get statistical summary of METBObj}
\usage{
get_stat_summary(METBObj)
}
\arguments{
\item{METBObj}{METBObj object contains list of data.}
}
\value{
data frame of statistical summary.
}
\description{
Get statistical summary of METBObj.
}
\details{
The function calculates median, mean, SD, coefficient of variation (cv) and normality.
}
\examples{
#sugar_dt = set_input_obj(sugar, 1,2,5)
#out=get_stat_summary(sugar_dt)
}
