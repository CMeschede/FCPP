#' Minimum Distance Functions
#'
#' Functions that calculate the Cramér-von-Mises (or modifications of it) distance
#' between the empirical distribution funtion
#' of given waiting times and the Cumulative Distribution function
#' of a mixture distribution of a dirac measure in zero and a
#' Mittag-Leffler distribution.
#'
#'@aliases distance distance_cm distance_cm1 and distance_cm2
#'@param WW A data.frame or a tibble containing inter exceedance times
#'@param tail = \eqn{\beta}; tail parameter
#'@param ei = \eqn{\theta}; extremal index
#'@param ptail = P(JJ > u) = k / n (number of exceedances / number of observations)
#'
#'
#'@details
#'For the cumulative distribution function
#' \eqn{F^*_{\beta,\theta}} of a \eqn{ML\big(\beta,\theta^{-1/\beta}\big)}-distribution,
#' the cumulative distribution function of the mixture distribution is given as
#'\deqn{F_{\beta,\theta}(x)=(1-\theta)\cdot\text{I}_{[0,\infty)}
#'(x)+\theta \cdot F^*_{\beta,\theta}(x).}
#'The Cramér von Mises distance between the empirical distribution function
#'\eqn{F_{k_*}} of \code{WW}
#'and \eqn{F^*_{\beta,\theta}} is defined as
#'\deqn{\Delta^{[CM]}(F_{k_*},F_{\beta,\theta})
#'=\int\limits_0^\infty(F_{k_*}(x)-F_{\beta,\theta}(x))^2
#'dF_{\beta,\theta}(x)
#'}
#' where \eqn{k_*=}\code{length(WW)}.
#' The function \code{distance_cm_mod1} and \code{distance_cm_mod2}
#' calculating the modified version of
#' the Cramér von Mises Distance functions
#' \deqn{
#' \Delta^{[CMmod1]}(F_{k_*},F_{\beta,\theta})
#' =\int\limits_0^\infty(F_{k_*}(x)-F_{\beta,\theta}(x))^2
#' dF^*_{\beta,\theta}(x)
#' }
#' respectively
#' \deqn{
#' \Delta^{[CMmod2]}(F_{k_*},F_{\beta,\theta})
#' =\int\limits_0^\infty(\max\lbrace F_{k_*}(x),1-\theta\rbrace
#' -F_{\beta,\theta}(x))^2
#' dF^*_{\beta,\theta}(x)
#' }.
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
