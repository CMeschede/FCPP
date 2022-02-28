#' Distribution functions and random number generation
#'
#' Probability density, cumulative distribution function, quantile function and
#' random variate generation for a mixture distribuion with the Dirac measure at
#' point zero and the Mittag-Leffler distribution as components.
#'
#' Now the details are coming... coming soon.
#' @seealso \code{\link[MittagLeffleR]{MittagLeffleR}} for the Mittag-Leffler functions
#' @aliases mixdistr rmixdistr pmixdistr qmixdistr dmixdistr
#' @param n number of observations.
#' @param p vector of probabilities.
#' @param x,q vector of quantiles.
#' @param tail tail parameter
#' @param ei extremal index / the weighting
#' @param scale scale parameter, default NULL sets scale = ei ^ (-1 / tail)
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
rmixdistr <- function(n, tail, ei, scale = NULL) {
  if (is.null(scale)) {
    scale <- ei ^ (-1 / tail)
  }
  if (tail == 1) {
    r <- ifelse(runif(n) > ei, 0, rexp(n, rate = 1 / scale))
  } else{
    r <- ifelse(runif(n) > ei, 0, MittagLeffleR::rml(n , tail = tail ,
                                                     scale = scale )
    )
  }
  return(r)
}

#' @rdname mixdistr
#' @export
pmixdistr <-  function(q, tail, ei, scale = NULL) {
  if (is.null(scale)) {
    scale <- ei ^ (-1 / tail)
  }
  if (tail == 1) {
    p <- 1 - ei + ei * pexp(q, rate = 1 / scale)
    p[q < 0] <- 0
  } else {
    p <- numeric(0)
    if (length(q[q >= 0]) > 0) {
      p[q >= 0] <- 1 - ei + ei * MittagLeffleR::pml(q[q >= 0], tail = tail,
                                                    scale = scale)
    }
    p[q < 0] <- 0
  }
  return(p)
}

#' @rdname mixdistr
#' @export
qmixdistr <-  function(p, tail, ei, scale = NULL) {
  if (is.null(scale)) {
    scale <- ei ^ (-1 / tail)
  }
  if (tail == 1) {
    p[p < (1 - ei)] <- 1 - ei
    q <- qexp( (p - (1 - ei)) / ei, rate = 1/scale)
  } else {
    q <- numeric(length(p))
    if ( any(p > (1 - ei)) ) {
      q[p > (1 - ei) & p < 1] <- MittagLeffleR::qml((p[p > (1 - ei) & p < 1] -
                                                       (1 - ei)) / ei ,
                                                    tail = tail, scale = scale)
    }
    q[p == 1] <- Inf
  }
  return(q)
}

#' @rdname mixdistr
#' @export
dmixdistr <- function(x, tail, ei, scale = NULL) {
  if (is.null(scale)) {
    scale <- ei^{-1/tail}
  }
  if(tail == 1){
    d <- ei * dexp(x , rate = 1/scale)
    d[x==0] <- 1-ei
  } else {
    d <- numeric(length(x))
    if (any(x > 0)) {
      d[x > 0] <- ei * MittagLeffleR::dml(x[x > 0], tail = tail, scale = scale)
    }
    d[x == 0] <- 1-ei
  }
  return(d)
}
