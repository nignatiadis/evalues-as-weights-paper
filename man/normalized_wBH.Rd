% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{normalized_wBH}
\alias{normalized_wBH}
\title{Normalized wBH}
\usage{
normalized_wBH(Ps, Es, alpha, Storey = FALSE)
}
\arguments{
\item{Ps}{Numeric vector of unadjusted p-values.}

\item{Es}{Numeric vector of e-values}

\item{alpha}{Significance level at which to apply method}

\item{Storey}{Bool (default: FALSE): is the procedure pi0 adaptive or not?}
}
\value{
Binary vector of rejected/non-rejected hypotheses.
}
\description{
Normalized wBH
}
