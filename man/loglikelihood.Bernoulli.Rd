% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/commons.R
\name{loglikelihood.Bernoulli}
\alias{loglikelihood.Bernoulli}
\title{Dirichlet-Bernoulli loglikelihood}
\usage{
\method{loglikelihood}{Bernoulli}(modelBernoulli, V)
}
\arguments{
\item{V}{A test matrix V. Same as the original V but with NA
in the training values}
}
\description{
Compute likelihood of some matrix V  p(V* | V)
}
\details{
It uses the samples from WH and p(V* | V) = E[WH]
}
