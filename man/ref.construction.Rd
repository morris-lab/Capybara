% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HighResolutionReferenceConstruction.R
\name{ref.construction}
\alias{ref.construction}
\title{Reference Construction}
\usage{
ref.construction(sc, sc.aux, criteria)
}
\arguments{
\item{sc}{The single-cell resolution dataset}

\item{sc.aux}{The auxiliary data frame that annotate the single-cell resolutiond dataset}

\item{criteria}{The column name to use for construction of the reference}
}
\description{
This function constructs reference from single-cell resolution reference data to be used for quadratic programming calculation
}
\examples{
ref.construction(single.ref.mtx, single.aux.df, "cell.type")
}
