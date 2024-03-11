#'Optimization functions for known and unknown scale parameter
#'
#'An optimization function, that estimates the unkown
#'parameter of
#'the mixture distribution consisting of a
#'dirac measure in zero and a Mittag Leffler distribution.
#'The function is using a grid of multiple initialization values.
#'
#'
#' @param WW     a vector or tibble (data.frame) with one column containing the
#'positive valued sample
#' @param ptail  probability of exceeding P(JJ > u) = \eqn{p(u)}
#' @param distance_fct  Cramér-von Mises distance function calculating the distance between the asymptotic and empirical cdf;
#' can be chosen from [distance_cm], [distance_cm_mod1] and [distance_cm_mod2]
#' @param scale0 scaling parameter \eqn{\rho}
#' @param t vector of starting values in (\code{a_tail},1] for the tail parameter \eqn{\beta}
#' @param e vector of starting values in (\code{a_ei},1] for the extremal index \eqn{\theta}
#' @param s vector of starting values in \eqn{(0,\infty)} for the scaling parameter \eqn{\sigma = p(u)^{-1/\beta} \cdot \rho}
#' @param a_tail lowest possible estimation for the tail parameter
#' @param a_ei lowest possible estimation for the extremal index
#' @param method algorithm that is used by \code{stats::optim} to minimize the
#' distance function
#' @param ... Additional parameters passed on to [optim]
#'
#'
#' @details
#' These two optimization functions calculate the minimum distance estimations of
#' \eqn{\beta, \theta, \sigma}, the parameters of the mixture distribution with
#' cumulative probability function (c.d.f.).
#' \deqn{F_{\beta,\theta,\sigma}(x)=(1-\theta)\cdot\text{I}_{[0,\infty)}
#' (x)+\theta \cdot F^*_{\beta,\theta,\sigma}(x),}
#' where \eqn{F^*_{\beta,\theta,\sigma}} is the c.d.f. of the Mittag-Leffler distribution
#' with tail parameter \eqn{\beta \le 1} and scale parameter
#' \eqn{\sigma_* = \theta^{-1/\beta}\cdot \sigma > 0}
#' (short notation \eqn{\text{ML}(\beta,\sigma_*)}).
#'
#' Use \code{optim_bts}, if you do not have any prior information.
#' Use \code{optim_bt}, if you know the excess probability \eqn{p(u) =} \code{ptail}
#' and the partial scale parameter \eqn{\rho =} \code{scale0}. Then, \eqn{\sigma = p(u)^{-1/\beta} \cdot \rho}
#' only depends on the unkown tail parameter \eqn{\beta} and we can reduce the search from three to two parameters.
#'
#' All combinations of \code{e}, \code{t} (and \code{s}) are used as starting values
#' to compute the minimum of \code{distance_fct}.
#' The parameters corresponding to the lowest value found, are returned
#' as estimation for \code{ei} and \code{tail}.
#'
#' @return
#' A tibble with 7 or 5 columns and one row
#'  - \code{t}, \code{e} (and \code{s}) are the starting points for the optimization algorithm that has led to the result.
#'  - \code{beta}, \code{theta} (and \code{sigma}) are the estimations for the parameter \eqn{\beta, \theta} (and \eqn{\sigma}).
#'  - \code{value} is the distance between the sample and the estimated mixture distribution regarding the chosen minimum distance function \code{distance_fct}.
#'

#'
#' @examples
#' n <- 10000
#' k <- 100
#' p <- k / n
#' beta  <- 0.9
#' theta <- 0.7
#' rho <- 2
#' sigma <- p^{-1/beta} * rho
#' # true parameter values:
#' c(beta, theta, sigma)
#'
#' # Data generation:
#' dat <- data_generation(n = n, stability = beta, ei = theta, scale0 = rho,
#'           wait_dist = "ML")
#' # Thinning the data:
#' dat_thinned <- thin(dat, k = k)
#' # waiting times:
#' WW <- interarrivaltime(dat_thinned)
#'
#' # Optimization:
#' optim_bts(WW = WW, distance_fct = distance_cm_mod1)
#' optim_bt(WW = WW, distance_fct = distance_cm_mod1, ptail = p, scale0 = rho)
#
#' @name optim
#' @export
optim_bts <- function(WW, distance_fct,
                      t = c(0.25,0.55,0.85), e = c(0.25,0.55,0.85), s = NULL,
                      a_tail = 0.1, a_ei = 0.1, method = "L-BFGS-B", ...) {
  if(is.data.frame(WW) & length(WW) == 1) {
    WW <- dplyr::pull(WW)
  }
  # Starting points of beta (b) and theta (t):
  start <- tidyr::crossing(
    t = t,
    e = e
  )
  # + starting points of sigma (s)
  # logarithmized to expand the searching space from [0, inf) to (-inf, -inf)
  if(is.numeric(s)) {
    start <- start |> tidyr::crossing(
      s = log(s) # starting point if specified by the user
    )
  } else if(is.null(s)) {
    # logMoment estimation for the Mittag-Leffler distribution
    initial <- MittagLeffleR::logMomentEstimator(WW)
    # initial[1] = tail parameter (beta); initial[2] = scale parameter (sigma)
    start <- start |> tidyr::crossing(
      s = log(initial[2]) # starting point if not specified by the user
    )
  } else {
    stop("s should be a numeric vector")
  }
  est <- dplyr::mutate(start,
                       opt = purrr::pmap(list(t, e, s), function(t., e., s.) {
                         tryCatch(
                           stats::optim( par = c(t., e., s.), fn = function(x) {
                             do.call(
                               distance_fct, list(WW = WW,
                                                  tail = x[1], ei = x[2], scale = exp(x[3]))
                             )
                           },
                           lower = c(a_tail, a_ei, -Inf), upper = c(1, 1, Inf),
                           method = method, ...),
                           error = function(x) NA
                         )
                       }),
                       beta = purrr::map_dbl(opt, ~tryCatch(.x$par[1],
                                                            error = function(x) NA )),
                       theta = purrr::map_dbl(opt, ~tryCatch(.x$par[2],
                                                             error = function(x) NA )),
                       sigma =  purrr::map_dbl(opt, ~tryCatch(exp(.x$par[3]), # rescaling
                                                              error = function(x) NA )),
                       value = purrr::map_dbl(opt, ~tryCatch(.x$value,
                                                             error = function(x) NA ))
  )
  est <- dplyr::select(est, -opt)
  est2 <- est[which.min(est$value), ]
  return(est2)
}

#' @rdname optim
#' @export
optim_bt <- function(WW, distance_fct, ptail, scale0,
                                       t = c(0.25,0.55,0.85), e = c(0.25,0.55,0.85),
                                       a_tail = 0.1, a_ei = 0.1, method = "L-BFGS-B") {
  if(is.data.frame(WW) & length(WW) == 1) {
    WW <- dplyr::pull(WW)
  }
  start <- tidyr::crossing(
    t = t,
    e = e
  )
  est <- dplyr::mutate(start,
                       opt = purrr::pmap(list(t, e), function(t., e.) {
                         tryCatch(
                           stats::optim( par = c(t., e.), fn = function(x) {
                             do.call(
                               distance_fct, list(WW = WW, tail = x[1], ei = x[2],
                                                  scale = ptail ^ {-1 / x[1]} * scale0)
                             )
                           },
                           lower = c(a_tail, a_ei), upper = c(1, 1),
                           method = method),
                           error = function(x) NA
                         )
                       }),
                       beta = purrr::map_dbl(opt, ~tryCatch(.x$par[1],
                                                            error = function(x) NA )),
                       theta = purrr::map_dbl(opt, ~tryCatch(.x$par[2],
                                                             error = function(x) NA )),
                       value = purrr::map_dbl(opt, ~tryCatch(.x$value,
                                                             error = function(x) NA ))
  )
  est <- dplyr::select(est, -opt)
  est2 <- est[which.min(est$value), ]
  return(est2)
}
