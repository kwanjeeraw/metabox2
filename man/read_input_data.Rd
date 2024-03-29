% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_input_data.R
\name{read_input_data}
\alias{read_input_data}
\title{Read CSV file}
\usage{
read_input_data(filename)
}
\arguments{
\item{filename}{the file path to the CSV file.
The input file contains data, which samples are in rows and variables/features are in columns separated by comma (,).
Information such as batch, replication, sample type and injection order, each must be in a separate column.}
}
\value{
data frame; row: variables/features, col: samples
}
\description{
Read CSV file
}
\examples{
#out = read_input_data('filename')
#fat = set_input_obj(read_input_data('data/repeated_samples_classification.csv'),2,3,4)
#lung = set_input_obj(read_input_data('data/independent_samples_classification.csv'),,5,6)
#freelive = set_input_obj(read_input_data('data/regression_example.csv'),1,2,3)
}
\author{
Kwanjeera W \email{kwanjeera.wan@mahidol.ac.th}
}
