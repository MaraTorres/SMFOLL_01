% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/th.R
\name{test}
\alias{test}
\title{Hypothesis test.}
\usage{
test(x, dist1, result1, dist2 = NULL, result2 = NULL, status = F)
}
\arguments{
\item{x}{Vector of quantiles.}

\item{dist1}{First distribution.}

\item{result1}{First result.}

\item{dist2}{Second distribution.}

\item{result2}{Second result.}

\item{status}{Data set with or without censure.}
}
\value{
Statistics, IC and Hypothesis test.
}
\description{
Provides the comparison between two fits
}
\examples{
data(atuaria)
result1=oll.fit(dist = 'bxii',starts = c(1,1,1),data = atuaria$x)
test(atuaria$x,'bxii',result1)
result2=oll.fit(dist = 'ollbxii',starts = c(result1$par,1),data = atuaria$x)
test(atuaria$x,'bxii',result1,'ollbxii',result2)
}
