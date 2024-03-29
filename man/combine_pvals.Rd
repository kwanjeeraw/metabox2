% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{combine_pvals}
\alias{combine_pvals}
\title{Combine p-values}
\usage{
combine_pvals(x)
}
\arguments{
\item{x}{numeric vectors of p-values}
}
\value{
numeric number of combined p-values
}
\description{
The function can combine p-values using Fisher's method.
}
\examples{
#result <- combine_pvals(c(0.01,0.005,0.1))
}
\references{
Fisher R. (1932) Statistical methods for research workers. Oliver and Boyd, Edinburgh.
}
\seealso{
\code{\link{pchisq}}
}
\author{
Kwanjeera W \email{kwanich@ucdavis.edu}
}
