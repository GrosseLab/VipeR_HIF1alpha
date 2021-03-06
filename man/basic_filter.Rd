% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/general_functions.R
\name{basic_filter}
\alias{basic_filter}
\title{Basic row filter}
\usage{
basic_filter(row, min_reads = 5, min_prop = 0.47)
}
\arguments{
\item{row}{this is a vector of numerics that will be passedin}

\item{min_reads}{the minimum mean number of reads}

\item{min_prop}{the minimum proportion of reads to pass this filter}
}
\value{
a logical of length 1
}
\description{
A basic filter to be used.
}
\note{
https://github.com/pachterlab/sleuth/blob/048f0551a31c4aee6e59b75c86cab46ae1b3ca3a/R/sleuth.R#L28
}
