% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FastMultiPCF.R
\name{calcScoresSegs}
\alias{calcScoresSegs}
\title{Takes a vector of values and marks positions at which the respective elements
are larger than a given quantile}
\usage{
calcScoresSegs(dist, sum, gamma)
}
\arguments{
\item{dist}{Vector with distances between breakpoints}

\item{sum}{matrix with summed intensities between each breakpoint}

\item{gamma}{threshold}
}
\value{
the modified vector of potential breakpoints
}
\description{
Takes a vector of values and marks positions at which the respective elements
are larger than a given quantile
}
