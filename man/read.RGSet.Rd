% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cumul.CNV.R
\name{read.RGSet}
\alias{read.RGSet}
\title{Reads the specified methylation arrays into an RGSet. Note that any combination of 450k and EPIC arrays will be coerced into an RGSet of 450k type.}
\usage{
read.RGSet(dataFiles, ArrayType)
}
\arguments{
\item{dataFiles}{Dataframe with a column batch name as requested by minfi for reading in experiments.}

\item{ArrayTye}{A string (either "450k", "EPIC", "combined" or "EPIC2")}
}
\value{
A list of the RGSet of the target data, the control data and the annotation data
}
\description{
Reads the specified methylation arrays into an RGSet. Note that any combination of 450k and EPIC arrays will be coerced into an RGSet of 450k type.
}
