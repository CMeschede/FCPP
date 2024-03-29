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
rmixdistr(n, tail, ei, scale = 1)

pmixdistr(q, tail, ei, scale = 1, lower.tail = TRUE, log.p = FALSE)

qmixdistr(p, tail, ei, scale = 1, lower.tail = TRUE, log.p = FALSE)

dmixdistr(x, tail, ei, scale = 1, log.p = FALSE)
}
\arguments{
\item{n}{number of observations.}

\item{tail}{tail parameter \eqn{\beta}}

\item{ei}{extremal index / weighting \eqn{\theta}}

\item{scale}{scale parameter \eqn{\sigma} (default \code{scale = 1}) sets
\eqn{\sigma_* = \theta ^ {-1 / \beta}  \sigma},
where \eqn{\sigma_*} is the scale parameter of the Mittag-Leffler distribution}

\item{lower.tail}{logical; if TRUE, probabilities are \eqn{P(X \le x)},
otherwise, \eqn{P(X>x)}.}

\item{log.p}{logical; if \code{TRUE}, probabilities \eqn{p} are given as log(p)}

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
random variate generation for a mixture distribution with the Dirac measure at
point zero and the Mittag-Leffler distribution as components.
}
\details{
The mixed distribution is a mixture of a Dirac measure
at point zero and a Mittag-Leffler distribution with tail parameter
\eqn{\beta} and a scale parameter \eqn{\sigma_*} (short notation \eqn{\text{ML}(\beta,\sigma_*)}),
where \eqn{\sigma_*=\theta^{-1/\beta} \sigma}.

The c.d.f. of the mixture distribution is given by
\deqn{F_{\beta,\theta,\sigma}(x)=(1-\theta) 1_{[0,\infty)}(x)
+ \theta F^*_{\beta,\sigma_*}(x),}
where \eqn{F^*} is the c.d.f. of the Mittag-Leffler distribution.
For clarification, \eqn{\sigma} corresponds to the input \code{scale} and
\eqn{\sigma_*=\theta^{-1/\beta} \sigma}.
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
