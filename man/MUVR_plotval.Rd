% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MUVR_util.R
\name{MUVR_plotval}
\alias{MUVR_plotval}
\title{Validation plot}
\usage{
MUVR_plotval(MVObj,plot_title="")
}
\arguments{
\item{MVObj}{MVObj object contains list of data.}

\item{plot_title}{text indicating plot title.}
}
\value{
ggplot object
}
\description{
Plot validation metric from MUVR object.
}
\details{
Overriding \code{\link[MUVR:plotVAL]{MUVR::plotVAL()}}.
}
\examples{
#res_muvr = run_muvr(METBObj, method="PLS")
#MUVR_plotval(res_muvr)
}
\references{
Lin Shi, et al. MUVR (2019) \url{https://doi.org/10.1093/bioinformatics/bty710}.
}
\seealso{
\code{\link[MUVR:plotVAL]{MUVR::plotVAL()}}
}
\author{
Kwanjeera W \email{kwanjeera.wan@mahidol.ac.th}
}
