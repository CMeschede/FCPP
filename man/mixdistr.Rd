% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mixdistr_functions.R
\name{mixdistr}
\alias{mixdistr}
\alias{rmixdistr}
\alias{pmixdistr}
\alias{qmixdistr}
\alias{dmixdistr}
\title{Distribution functions and random number generation}
\usage{
rmixdistr(n, tail, ei, scale = NULL)

pmixdistr(q, tail, ei, scale = NULL)

qmixdistr(p, tail, ei, scale = NULL)

dmixdistr(x, tail, ei, scale = NULL)
}
\arguments{
\item{n}{number of observations.}

\item{tail}{tail parameter}

\item{ei}{extremal index / the weighting}

\item{scale}{scale parameter, default NULL sets scale = ei ^ (-1 / tail)}

\item{p}{vector of probabilities.}

\item{x, q}{vector of quantiles.}
}
\value{
\code{rmixdistr} generates random variables,
\code{pmixdistr} returns the distribution function,
\code{qmixdistr} returns the quantile function, and
\code{dmixdistr} returns the density.
}
\description{
Probability density, cumulative distribution function, quantile function and
random variate generation for a mixture distribuion with the Dirac measure at
point zero and the Mittag-Leffler distribution as components.
}
\details{
Now the details are coming... coming soon.
}
\examples{
rmixdistr(n = 10, tail = 0.8, ei = 0.5)
pmixdistr(q = 10, tail = 0.8, ei = 0.5)
qmixdistr(p = 0:10 / 10, tail = 0.8, ei = 0.5)
dmixdistr(x = c(0, 1, 10), tail = 0.8, ei = 0.5)
}
\seealso{
\code{\link[MittagLeffleR]{MittagLeffleR}} for the Mittag-Leffler functions
}
