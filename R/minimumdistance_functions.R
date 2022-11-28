#' Minimum Distance Functions
#'
#' Functions that calculate the Cramer-von Mises distance (with modifications)
#' between the empirical distribution function (e.c.d.f.) of a positive valued sample \code{WW}
#' and a mixture distribution with the dirac measure in zero and Mittag-Leffler distribution
#' as parts.
#' The sample \code{WW} may be return times of events of interests.
#'
#'@aliases distance distance_cm distance_cm1 and distance_cm2
#'@param WW a vector or tibble (data.frame) with one column containing the
#'positive valued sample
#'@param tail tail parameter \eqn{\beta} of the Mittag-Leffler distribution
#'@param ei extremal index / weighting \eqn{\theta}
#'@param scale scale parameter \eqn{\sigma}
#'
#'
#'@details
#' The cumulative distribution function (c.d.f.) of the mixture distribution is given as
#'\deqn{F_{\beta,\theta,\sigma}(x)=(1-\theta)\cdot\text{I}_{[0,\infty)}
#'(x)+\theta \cdot F^*_{\beta,\theta,\sigma}(x),}
#'where \eqn{F^*_{\beta,\theta,\sigma}} is the c.d.f. of the Mittag-Leffler distribution
#'with tail-parameter \eqn{\beta <= 1} and scale parameter
#'\eqn{\sigma* = \theta^{-1/\beta}\cdot \sigma > 0}
#' (short notation \eqn{ML(\beta,\sigma*)}).
#'The Cramer-von Mises distance between the e.c.d.f. \eqn{F_{n}} of \eqn{WW_1,\dots,WW_n}
#'and the mixture distribution is then defined as
#'\deqn{\Delta^{[CM]}(F_{n},F_{\beta,\theta,\sigma})
#'=\int\limits_0^\infty(F_{n}(x)-F_{\beta,\theta,\sigma}(x))^2
#'dF_{\beta,\theta,\sigma}(x)
#'}
#'\deqn{
#'=\frac{1}{n}\sum\limits_{i=1}^n
#' \bigg(
#' \frac{i-\frac{1}{2}}{n}
#' -F_{\beta,\theta,\sigma}(t_{(i+1)})
#' \bigg)^2 +
#' \frac{1}{12n^2}+\frac{2}{3}(1-\theta)^3
#'}
#' where \eqn{n=}\code{length(WW)}.
#' The function \code{distance_cm_mod1} and \code{distance_cm_mod2}
#' calculate two modified version of the Cramer-von Mises distance
#' \deqn{
#' \Delta^{[CMmod1]}(F_{n},F_{\beta,\theta,\sigma})
#' =\int\limits_0^\infty(F_{n}(x)-F_{\beta,\theta,\sigma}(x))^2
#' dF^*_{\beta,\theta,\sigma}(x)
#' }
#' \deqn{
#'=\frac{1}{\theta}\frac{1}{n}\sum\limits_{i=1}^n
#' \bigg(
#' \frac{i-\frac{1}{2}}{n}
#' -F_{\beta,\theta,\sigma}(t_{(i+1)})
#' \bigg)^2 +
#' \frac{1}{12n^2}\frac{1}{\theta}+\frac{1}{3}\frac{(1-\theta)^3}{\theta}
#'}
#' and
#' \deqn{
#' \Delta^{[CMmod2]}(F_{n},F_{\beta,\theta,\sigma})
#' =\int\limits_0^\infty(\max\lbrace F_{n}(x),1-\theta\rbrace
#' -F_{\beta,\theta,\sigma}(x))^2
#' dF^*_{\beta,\theta,\sigma}(x)
#' }
#' \deqn{
#'=\frac{1}{\theta^3}\frac{1}{n}\sum\limits_{i=
#'\lceil n(1-\theta)\rceil+1}^n
#' \bigg(
#' \frac{i-\frac{1}{2}}{n}
#' -F_{\beta,\theta,\sigma}(t_{(i+1)})
#' \bigg)^2 +
#' \frac{n-\lceil n(1-\theta)\rceil}{12n^3\theta^3}
#'}
#'\deqn{
#'-\frac{(n(1-\theta)^3)-\lceil n(1-\theta)\rceil^3}
#'{3n^3\theta^3}
#'+
#'\frac{(n(1-\theta))^2-\lceil n(1-\theta)\rceil^2}
#'{n^2\theta^3}
#'F_{\beta,\theta,\sigma}(t_{(\lceil n(1-\theta)\rceil +1)})
#'}
#'\deqn{
#'-\frac{n(1-\theta)-\lceil n(1-\theta)\rceil}{n\theta^3}
#'F_{\beta,\theta,\sigma}(t_{(\lceil n(1-\theta)\rceil +1)})^2
#'}
#' respectively.
#'
#' @name distance

#' @export
distance_cm <- function(WW, tail, ei, scale = 1) {
  # WW = IETs
  # tail = beta; tail parameter
  # ei = theta; extremal index
  # ptail = P(JJ > u) = k / n (number of exceedances / number of observations)
  if(is.data.frame(WW) & length(WW) == 1) {
    WW <- dplyr::pull(WW)
  }
  if(any(WW < 0)){
    stop("WW (waiting times) should be positive values")
  }
  kstar <- length(WW) # number of IETs
  k <- kstar + 1 # number of exceedances
  WW_sort <- sort(WW, decreasing = F)
  pn <- stats::ppoints(n = kstar, a = 0.5) # (1:kstar - 0.5)/kstar
  pmisch <- pmixdistr(WW_sort, tail = tail, ei = ei, scale = scale) # F_{beta,theta, sigma}(t_(i))
  dist <- mean((pn - pmisch) ^ 2)
  s <- dist + 1 / (12 * (kstar ^ 2)) + (2 / 3) * (1 - ei) ^ 3
  return(s)
}

#' @rdname distance
#' @export
distance_cm_mod1 <- function(WW, tail, ei, scale = 1) {
  if(is.data.frame(WW) & length(WW) == 1) {
    WW <- dplyr::pull(WW)
  }
  if(any(WW < 0)){
    stop("WW (waiting times) should be positive values")
  }
  kstar <- length(WW) # number of IETs
  k <- kstar + 1 # number of exceedances
  WW_sort <- sort(WW, decreasing = F) # IETs
  pn <- stats::ppoints(n = kstar, a = 0.5) # (1:kstar - 0.5)/kstar
  pmisch <- pmixdistr(WW_sort, tail = tail, ei = ei, scale = scale) # F_{beta,theta, sigma}(t_(i))
  dist <- mean((pn - pmisch) ^ 2)
  s <- (dist + 1 / (12 * (kstar ^ 2)) - ((1 - ei) ^ 3) / 3) / ei
  return(s)
}

#' @rdname distance
#' @export
distance_cm_mod2 <- function(WW, tail, ei, scale = 1) {
  if(is.data.frame(WW) & length(WW) == 1) {
    WW <- dplyr::pull(WW)
  }
  if(any(WW < 0)){
    stop("WW (waiting times) should be positive values")
  }
  kstar <- length(WW) # number of IETs
  k <- kstar + 1 # number of exceedances
  WW_sort <- sort(WW, decreasing = F) # IETs
  l <- kstar * (1 - ei)
  m <- ceiling(l) # lceil kstar*(1 - theta) rceil
  if(m == kstar) {
    pml_maxWW <- MittagLeffleR::pml(max(WW_sort),
                                    tail = tail, scale = ei ^ (-1  /tail) * scale)
    s <- 1/3 - pml_maxWW + pml_maxWW ^ 2
    return(s)
  }
  if(m < kstar) {
    if(m == 0) {
      cdf_WW_m <- 0
    } else {
      cdf_WW_m <- pmixdistr(WW_sort[m], tail = tail, ei = ei, scale = scale)
    }
    WW_trunc <- WW_sort[(m + 1):kstar]
    pn_trunc <- (stats::ppoints(n = kstar, a = 0.5))[(m + 1):kstar] # (1:kstar - 0.5)/kstar
    pmisch_trunc <- pmixdistr(WW_trunc, tail = tail, ei = ei, scale = scale) # F_{beta,theta}(t_(i))
    dist <- sum((pn_trunc - pmisch_trunc) ^ 2) / kstar
    s <- (dist
          + (kstar - m) / (12 * (kstar ^ 3))
          - (l ^ 3 - m ^ 3) / (3 * (kstar ^ 3))
          + (l ^ 2 - m ^ 2) / (kstar ^ 2) * cdf_WW_m
          - (l - m) / kstar * (cdf_WW_m ^ 2)) / (ei ^ 3)
    return(s)
  }
}

# not in the article:
distance_cm_mod3 <- function(WW, tail, ei, scale = 1) {
  if(is.data.frame(WW) & length(WW) == 1) {
    WW <- dplyr::pull(WW)
  }
  if(any(WW < 0)){
    stop("WW (waiting times) should be positive values")
  }
  kstar <- length(WW) # number of IETs
  k <- kstar + 1 # number of exceedances
  WW_sort <- sort(WW, decreasing = F) # scaled IETs
  pn <- stats::ppoints(n = kstar, a = 0.5) # (1:kstar - 0.5)/kstar
  pmisch <- pmixdistr(WW_sort, tail = tail, ei = ei, scale = scale) # F_{beta,theta}(t_(i))
  dist <- mean((pn - pmisch) ^ 2)
  s <- (dist + 1 / (12 * (kstar ^ 2)) - ((1 - ei) ^ 3) / 3) /
    ((0.5 - 1 / 3 * ei) * ei ^ 2)
  s <- (1 / ((0.5 - 1 / 3 * ei) * ei)) *
    sum((((sort(c(1 - ei, pmisch)) - pn) ^ 2) * weights))
  return(s)
}

## minimum distance functions approximations

#' @export
distance_cm_approx <- function(WW, tail, ei, scale = 1) {
  if(is.data.frame(WW) & length(WW) == 1) {
    WW <- dplyr::pull(WW)
  }
  if(any(WW < 0)){
    stop("WW (waiting times) should be positive values")
  }
  k <- length(WW) + 1
  WW_sort <- sort(WW, decreasing  = F) # IETs
  pn <- (0:(k - 1)) / k
  pmls <- MittagLeffleR::pml(WW_sort, tail = tail, scale = ei ^ (-1 / tail) * scale)
  pmisch <- (1 - ei) + ei * pmls
  weights <- diff(c(0, pmls, 1))
  s <- (1 - ei) ^ 3 + ei * sum((((c(1 - ei, pmisch) - pn) ^ 2) * weights))
  return(s)
}

#' @export
distance_cm_mod1_approx <- function(WW, tail, ei, scale = 1) {
  if(is.data.frame(WW) & length(WW) == 1) {
    WW <- dplyr::pull(WW)
  }
  if(any(WW < 0)){
    stop("WW (waiting times) should be positive values")
  }
  k <- length(WW) + 1
  WW_sort <- sort(WW, decreasing = F) # IETs
  pn <- (0:(k - 1)) / k
  pmls <- MittagLeffleR::pml(WW_sort, tail = tail, scale = ei ^ (-1 / tail) * scale)
  pmisch <- (1 - ei) + ei * pmls
  weights <- diff(c(0, pmls, 1))
  s <- sum((((c(1 - ei, pmisch) - pn) ^ 2) * weights))
  return(s)
}

#' @export
distance_cm_mod2_approx <- function(WW, tail, ei, scale = 1) {
  if(is.data.frame(WW) & length(WW) == 1) {
    WW <- dplyr::pull(WW)
  }
  if(any(WW < 0)) {
    stop("WW (waiting times) should be positive values")
  }
  k <- length(WW) + 1
  WW_sort <- sort(WW, decreasing = F) # sorted IETs
  pn <- (0:(k - 1)) / k
  pmls <- MittagLeffleR::pml(WW_sort, tail = tail, scale = ei ^ (-1  /tail) * scale)
  weights <- diff(c(0, pmls, 1))
  s <- sum((c(0,pmls) -
              pmax(0, pmin(1, (1 / ei * (pn - 1 + ei)))))  ^ 2  * weights)
  return(s)
}

# not in the article:
distance_cm_mod3_approx <- function(WW, tail, ei, scale = 1) {
  if(is.data.frame(WW) & length(WW) == 1) {
    WW <- dplyr::pull(WW)
  }
  if(any(WW < 0)){
    stop("WW (waiting times) should be positive values")
  }
  k <- length(WW) + 1
  WW_sort <- sort(WW, decreasing = F) # scaled IETs
  pn <- (0:(k - 1)) / k
  pmls <- MittagLeffleR::pml(WW_sort, tail = tail, scale = ei ^ (-1 / tail) * scale)
  pmisch <- (1 - ei) + ei * pmls
  weights <- diff(c(0, pmls, 1))
  s <- (1 / ((0.5 - 1 / 3 * ei) * ei)) *
    sum((((c(1 - ei, pmisch) - pn) ^ 2) * weights))
  return(s)
}
