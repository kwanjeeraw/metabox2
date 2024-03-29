% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/overrep_analyze.R
\name{overrep_analyze}
\alias{overrep_analyze}
\title{Perform overrepresentation analysis}
\usage{
overrep_analyze(txtinput, nodetype="compound", settype="pathway",
organism="hsa", size=3, universe)
}
\arguments{
\item{txtinput}{a character vector of variable/feature IDs e.g. c('C12078', 'C02273'). See details below.}

\item{nodetype}{a string specifying a node type. It can be one of compound (default), protein, gene.}

\item{settype}{a string specifying a set type. It can be one of pathway (default), chemicalclass.}

\item{organism}{a string specifying organism code from KEGG database, the parameter will not be used for \code{settype = "chemicalclass"}. Choose one from hsa (default), tdc.}

\item{size}{a number specifying the minimum number of members in each annotation term to be used in the analysis. Default = 3.}

\item{universe}{a character vector of variable/feature IDs that represent the universe. Default's to all unique IDs in a set collection.}
}
\value{
a list of the following components:

enrichment = enrichment results.

network = network format for \code{\link[piano:networkPlot]{piano::networkPlot()}}.

details = a list of analysis details: univsize, bgsize, testMethod, pAdjusted, inputsize, numsets, minsize.
}
\description{
perform overrepresentation analysis to retrieve overrepresented annotation terms of the input entities. The function wraps around the main functions of \pkg{\link{piano}}.
}
\details{
For pathway analysis, Metabox uses KEGG ID (e.g.C12078) for compounds, UniProt entry (e.g.P0C9J6) for proteins, Ensembl (e.g.ENSG00000139618) for genes.
For chemical class analysis, HMDB ID (e.g.HMDB0000001) is used for compounds.
}
\examples{
#out=overrep_analyze(fnanal_data$compound_data$kegg, organism = "tdc", size = 5) #pathway ORA
#out=overrep_analyze(fnanal_data$combined_data$id[1:9], nodetype = "protein") #pathway ORA
#out=overrep_analyze(fnanal_data$compound_data$hmdb, settype="chemicalclass") #chemical class ORA
}
\references{
Fisher R. (1932) Statistical methods for research workers. Oliver and Boyd, Edinburgh.

Väremo L., Nielsen J., and Nookaew I. (2013) Enriching the gene set analysis of genome-wide data by incorporating directionality of gene expression and combining statistical hypotheses and methods. Nucleic Acids Research, 41(8), pp. 4378-4391.
}
\seealso{
\code{\link[piano:loadGSC]{piano::loadGSC()}}, \code{\link[piano:runGSAhyper]{piano::runGSAhyper()}}, \code{\link[piano:GSAsummaryTable]{piano::GSAsummaryTable()}}
}
\author{
Kwanjeera W \email{kwanjeera.wan@mahidol.ac.th}
}
