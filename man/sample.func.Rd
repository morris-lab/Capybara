% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EmpiricalPValueCalculation.R
\name{sample.func}
\alias{sample.func}
\title{Empirical p-value of Identity Score Calculation}
\usage{
sample.func(dens.x, curr.val, prob, n = 1000)
}
\arguments{
\item{dens.x}{The background dataset that we compare the identity score of our sample with}

\item{curr.val}{The identity score of the sample dataset to be evaluated its significance}

\item{n}{The number of times of resampling. The default is n=1000}
}
\description{
This function resamples from a sample dataset and returns an empirical p-value of identity score
}
\examples{

}
\keyword{empirical}
\keyword{p-value}
\keyword{resampling,}
