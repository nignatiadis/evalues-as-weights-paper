% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{weighted_storey_pi0}
\alias{weighted_storey_pi0}
\title{Storey's pi0 estimator for weighted multiple testinng}
\usage{
weighted_storey_pi0(pvalues, weights = 1, tau = 0.5, m = length(pvalues))
}
\arguments{
\item{pvalues}{Numeric vector of unadjusted p-values.}

\item{weights}{Numeric vector of multiple testing weights}

\item{tau}{Numeric (default = 0.5), the level at which tau-censoring is applied.}

\item{m}{Total number of tests (default: `length(pvalues)`)}
}
\value{
Estimated null proportion
}
\description{
Storey's pi0 estimator for weighted multiple testinng
}
