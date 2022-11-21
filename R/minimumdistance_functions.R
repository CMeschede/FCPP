## minimum distance functions

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
