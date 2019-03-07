\name{GSE25055}
\alias{GSE25055}
\docType{data}
\title{
Example TNBC microarray data
}
\description{
This is a TNBC microarray dataset from GSE25055 contained in a
\code{SummarizedExperiment} object. It includes gene expression profiles and
clinical information which can be accessed by the \code{assays} and
\code{colData} functions, respectively.
We obtained gene expression profiles of entire breast cancer
samples from the \code{curatedBreastData} package and extracted
TNBC samples based on the expression profiles and
immunohistochemistry results.
}
\source{
https://bioconductor.org/packages/release/data/experiment/html/curatedBreastData.html
}
\references{
Hatzis, C. et al. (2011). A genomic predictor of response and survival following
taxane-anthracycline chemotherapy for invasive breast cancer \emph{JAMA}, 305, 1873-81.
}
\examples{
data(GSE25055)

#Access gene expression profiles
head(assays(GSE25055)[[1]])

#Access clinical information
head(colData(GSE25055))
}
\keyword{datasets}
