% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/edegR_function.R
\name{summarize_replicates}
\alias{summarize_replicates}
\title{Summarize replicates in an expression matrix}
\usage{
summarize_replicates(exp_mat, groups, method = mean, changeColNames = T)
}
\arguments{
\item{exp_mat}{is a matrix or a data.frame, in which expression of genes are described by rows}

\item{groups}{vector describing which columns should be summarized}

\item{method}{defines how the expression values shall be summarized. One of mean (default), median, geometric mean, min, max, sum or any other method that is implemented in R}

\item{changeColNames}{boolean for adjusting colnames}
}
\value{
a \code{vector} or \code{matrix} of normalized expression values per row
}
\description{
This method summarizes the replicated expression values in an expression matrix
}
\author{
Alexander Gabel,Claus Weinholdt
}
