% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cumFreq.R
\name{cumFreq}
\alias{cumFreq}
\title{Generates the multi sample segmentation aberration frequency plot}
\usage{
cumFreq(
  mSetsAnno,
  seg_mpcf,
  target_ratios,
  colour.amplification,
  colour.loss,
  thresh,
  array_type
)
}
\arguments{
\item{mSetsAnno}{A list of the RGSet of the target data, the control data and the annotation data}

\item{seg_mpcf}{a dataframe of the segmentation results of the multi sample segmentation.}

\item{target_ratios}{a dataframe of the intesity of all bins for each sample.}

\item{colour.amplification}{Colour for amplification}

\item{colour.loss}{Colour for loss}

\item{thresh}{a float specifying the threshold when a segment is called as an aberration}
}
\value{
returns the frequency plot of the samples that were segmented togethher using the multi sample segmentation algorithm
}
\description{
Generates the multi sample segmentation aberration frequency plot
}
