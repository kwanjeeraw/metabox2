% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_loess.R
\name{run_loess}
\alias{run_loess}
\title{LOESS normalization}
\usage{
run_loess(p, f, e, e_matrix)
}
\arguments{
\item{p}{is phenotypes. p is a data.frame containing following columns:
sampleType must include 'qc', 'sample'; It can include NA, which will not be normalized by LOESS, e.g. blank samples; It can include other type of samples like 'validate1', 'validate2', e.g. validate QCs.
time is the injection order of each sample. Must be numeric, and unique.
batch is batch information of samples. Each batch must contain at least 2 'qc's, although 5 qcs are recommended.}

\item{f}{is phenotypes.}

\item{e}{is the data matrix. Each row is a compound and each column is a sample.}

\item{e_matrix}{is the \code{e} data matrix. Each row is a compound and each column is a sample.}
}
\value{
matrix with samples in columns and compounds in rows.
}
\description{
Normalize samples by QC-based methods, LOESS
}
\examples{
#loess_dat = run_loess(p, f, e, e_matrix)
}
\references{
Sili Fan, et al. SERRF (2019) \url{https://doi.org/10.1021/acs.analchem.8b05592}.
}
\author{
Sili Fan \email{}
}
