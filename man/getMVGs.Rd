% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{getMVGs}
\alias{getMVGs}
\title{Get the most variable genes (MVGs)}
\usage{
getMVGs(profile = NULL, coefVar = 0.1, no = 2000)
}
\arguments{
\item{profile}{code{matrix} The expression profile}

\item{coefVar}{\code{numeric} A threshold value for the coefficient variation, default is 0.1}

\item{no}{\code{numeric} Number of genes to return, default is 2000}
}
\value{
A vector of genes
}
\description{
Get the most variable genes (MVGs)
}
\examples{
nmfObj <- getMVGs(profile = expr, coefVar = 0.1, no = 2000)

}
