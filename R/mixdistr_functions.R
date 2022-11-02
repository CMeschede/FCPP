#' Distribution functions and random number generation
#'
#' Probability density, cumulative distribution function, quantile function and
#' random variate generation for a mixture distribuion with the Dirac measure at
#' point zero and the Mittag-Leffler distribution as components.
#'
#' @seealso \code{\link[MittagLeffleR]{MittagLeffleR}} for the Mittag-Leffler functions
#' @aliases mixdistr rmixdistr pmixdistr qmixdistr dmixdistr
#' @param n number of observations.
#' @param p vector of probabilities.
#' @param x,q vector of quantiles.
#' @param tail tail parameter \eqn{\beta}
#' @param ei extremal index / the weighting \eqn{\theta}
#' @param scale scale parameter \eqn{\sigma} (default \code{scale = 1}) sets \eqn{sigma^\star = ei ^ (-1 / tail) * \sigma},
#'  where sigma^* is the scale parameter of the Mittag-Leffler distribution
#' @param log.p logical; if \code{TRUE}, probabilitied \eqn{p} are given as log(p)
#' @param lower.tail logical; if TRUE, probabilities are \eqn{P(X \le x)},
#'  otherwise, \eqn{P(X>x)}.
#'
#' @details The mixed distribution is a mixture of da Dirac measure
#' in 0 and a Mittag-Leffler distribution with tail parameter
#' \eqn{\beta} and a scale parameter \eqn{\sigma}. A Mittag-Leffler \eqn{ML(\beta,\sigma)}
#' distributed random variable \eqn{Z_{\beta}} with parameter
#' \eqn{\beta\in(0,1]} and \eqn{\sigma=1} can be defined on the positive real numbers
#' by its Laplace transform
#' \deqn{\mathcal{L}_{Z_{\beta}}(s)=\frac{1}{1+s^\beta}_.}
#' For an arbitrary \eqn{\sigma>0} write \eqn{\sigma Z_\beta}.
#' The mixed distribution has the c.d.f
#' \deqn{F_{\beta,\theta,\sigma}(x)=(1-\theta)\cdot 1_{[0,\infty)}(x)
#' +\theta\cdot F^\star_{\beta,\sigma^\star}(x),}
#' where \eqn{F^\star_{\beta,\sigma^\star}} is the c.d.f
#' of the \eqn{ML(\beta,\theta^{-1/\beta}\cdot\sigma)} distribution,
#' where \eqn{\sigma^\star} corresponds to the input \code{scale}.
#'
#'
#' @return \code{rmixdistr} generates random variables,
#' \code{pmixdistr} returns the distribution function,
#' \code{qmixdistr} returns the quantile function, and
#' \code{dmixdistr} returns the density.
#'
#' @examples
#' rmixdistr(n = 10, tail = 0.8, ei = 0.5)
#' pmixdistr(q = 10, tail = 0.8, ei = 0.5)
#' qmixdistr(p = 0:10 / 10, tail = 0.8, ei = 0.5)
#' dmixdistr(x = c(0, 1, 10), tail = 0.8, ei = 0.5)
#' @name mixdistr
#' @export


rmixdistr <- function(n, tail, ei, scale = 1) {
  stopifnot(as.integer(n) == n, tail > 0, tail <= 1, ei > 0, ei <= 1, scale > 0)
    scale1 <- ei ^ (-1 / tail) * scale

  if (tail == 1) {
    r <- ifelse(stats::runif(n) > ei, 0, stats::rexp(n, rate = 1 / scale1))
  } else{
    r <- ifelse(stats::runif(n) > ei, 0, MittagLeffleR::rml(n , tail = tail,
                                                     scale = scale1 ))
    }
  return(r)
}

#' @rdname mixdistr
#' @export
pmixdistr <-  function(q, tail, ei, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(tail > 0, tail <= 1, ei > 0, ei <= 1, scale > 0)
  scale1 <- ei ^ (-1 / tail) * scale
  if (tail == 1) {
    p <- 1 - ei + ei * stats::pexp(q, rate = 1 / scale1)
    p[q < 0] <- 0
  } else {
    p <- numeric(0)
    if (length(q[q >= 0]) > 0) {
      p[q >= 0] <- 1 - ei + ei * MittagLeffleR::pml(q[q >= 0], tail = tail,
                                                    scale = scale1)
    }
    p[q < 0] <- 0
  }
  if(!lower.tail){
    p <- 1 - p
  }
  if(log.p){
    p <- log(p)
  }
  return(p)
}

#' @rdname mixdistr
#' @export
qmixdistr <-  function(p, tail, ei, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(p >= 0, p <= 1, tail > 0, tail <= 1, ei > 0, ei <= 1, scale > 0)
  scale1 <- ei ^ (-1 / tail) * scale
  if(log.p){
    p <- log(p)
  }
  if(!lower.tail){
    p <- 1 - p
  }
  if (tail == 1) {
    p[p < (1 - ei)] <- 1 - ei
    q <- stats::qexp((p - (1 - ei)) / ei, rate = 1 / scale)
  } else {
    q <- numeric(length(p))
    if ( any(p > (1 - ei)) ) {
      q[p > (1 - ei) & p < 1] <- MittagLeffleR::qml((p[p > (1 - ei) & p < 1] -
                                                       (1 - ei)) / ei ,
                                                    tail = tail, scale = scale)}
    q[p == 1] <- Inf
  }
  return(q)
}

#' @rdname mixdistr
#' @export
dmixdistr <- function(x, tail, ei, scale = 1, log.p = FALSE) {
  stopifnot(tail > 0, tail <= 1, ei > 0, ei <= 1, scale > 0)
  scale1 <- ei^{-1/tail} * scale
  if(tail == 1) {
    d <- ei * stats::dexp(x , rate = 1 / scale1)
    d[x == 0] <- 1 - ei
  }
 else {
    d <- numeric(length(x))
    if(any(x > 0)) {
      d[x > 0] <- ei * MittagLeffleR::dml(x[x > 0], tail = tail, scale = scale1)
    }
    d[x == 0] <- 1 - ei
  }
  if(log.p){
    d <- log(d)
  }
  return(d)
}
