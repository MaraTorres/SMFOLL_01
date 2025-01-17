% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/halfnormal.R
\name{HalfNormal}
\alias{HalfNormal}
\alias{dhn}
\alias{phn}
\alias{qhn}
\alias{rhn}
\title{The Half-Normal Distribution.}
\usage{
dhn(param, x, log = FALSE)

phn(param, q, lower.tail = TRUE, log = FALSE)

qhn(param, p, lower.tail = TRUE, log = FALSE)

rhn(param, n)
}
\arguments{
\item{param}{Vector with shape parameter first and scale parameter second.}

\item{x, q}{Vector of quantiles.}

\item{log}{Logical; if TRUE, probabilities p are given as log(p).}

\item{lower.tail}{Logical; if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].}

\item{p}{Vector of probabilities.}

\item{n}{Number of observations.}
}
\value{
dhn gives the density, phn gives the distribution function, qhn quantile function, rhn random generation function.
}
\description{
Provides density, distribution function, quantile function and random generation for the half-normal distribution.
}
\examples{
phn(c(.4,2.1,3),5)

dhn(c(1,1),1)
}
