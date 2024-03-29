% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiv_analyze.R
\name{multiv_analyze}
\alias{multiv_analyze}
\title{Multivariate analysis}
\usage{
multiv_analyze(METBObj, method="pca", scale="center")
}
\arguments{
\item{METBObj}{METBObj object contains list of data.}

\item{method}{name of multivariate analysis method. Choose one from the list: pca, pls, opls. Default is pca.}

\item{scale}{text indicating a scaling method for parameter \code{scaleC} of \code{ropls::opls()}: no centering nor scaling ('none'),
mean-centering only ('center') [default], mean-centering and pareto scaling ('pareto'), or mean-centering and auto scaling ('standard').}
}
\value{
a list of the following components:

model_summary = a data frame with the model overview.

score_val = a numerical matrix of x scores.

loading_val = a numerical matrix of x loadings.

oscore_val = a numerical matrix of orthogonal scores.

oloading_val = a numerical matrix of orthogonal loadings.

vip_val = a numerical vector of VIP.

ovip_val = a numerical vector of variable importance for orthogonal modeling.

details = a list of analysis details: itestMethod, scale.
}
\description{
Perform multivariate analysis.
}
\examples{
#sugar_dt = set_input_obj(sugar, 1,2,5)
#out=multiv_analyze(sugar_dt)
}
\references{
Thevenot EA, et al. (2015) ropls.
}
\seealso{
\code{\link[ropls:opls]{ropls::opls()}}
}
\author{
Kwanjeera W \email{kwanjeera.wan@mahidol.ac.th}
}
