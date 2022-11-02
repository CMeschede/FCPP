# WW - vector of inter-exceedance times
# ptail - P(JJ > u) = k/n (number of exceedances over number of observations)
# distance_fct - r-function calculating the distance between the asymptotic and empirical cdf
#' @export
optim_multistart_bt <- function(WW, distance_fct, ptail, rho,
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
                                                  scale = ptail ^ {-1 / x[1]} * rho)
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


#' @export
optim_multistart_btr <- function(WW, distance_fct, ptail = NULL, type = 1,
                                            t = c(0.25,0.55,0.85), e = c(0.25,0.55,0.85),
                                            r = NULL, s = NULL,
                                            a_tail = 0.1, a_ei = 0.1, method = "L-BFGS-B") {
  if(is.data.frame(WW) & length(WW) == 1) {
    WW <- dplyr::pull(WW)
  }
  initial <- MittagLeffleR::logMomentEstimator(WW) # for the starting point of rho
  if(type == 1) {
  if(is.null(ptail)) stop("for type = 1 you need to give the tail probability ptail")
  start <- tidyr::crossing(
    t = t,
    e = e
  )
  if(is.numeric(r)) {
    start <- start |> tidyr::crossing(
      r = r
    )
  } else if(is.numeric(s)) {
    start <- start |> tidyr::crossing(
      r = s * ptail ^ {1 / t}
    )
  } else {
    start <- start |> tidyr::crossing(
      r = initial[2] * ptail ^ {1 / t}
    )
  }


  # para_orig: parameter (tail, ei, scale)
  # para: a transformed parameter so that the optimation is unconstrained
  #para_orig <- function(para) {c()}
  est <- dplyr::mutate(start,
                       opt = purrr::pmap(list(t, e, r), function(t., e., r.) {
                         tryCatch(
                           stats::optim( par = c(t., e., r.), fn = function(x) {
                             do.call(
                               distance_fct, list(WW = WW, tail = x[1], ei = x[2],
                                                  scale = ptail ^ {-1 / x[1]} * x[3])
                             )
                           },
                           lower = c(a_tail, a_ei, 0), upper = c(1, 1, Inf),
                           method = method),
                           error = function(x) NA
                         )
                       }),
                       beta = purrr::map_dbl(opt, ~tryCatch(.x$par[1],
                                                            error = function(x) NA )),
                       theta = purrr::map_dbl(opt, ~tryCatch(.x$par[2],
                                                             error = function(x) NA )),
                       rho = purrr::map_dbl(opt, ~tryCatch(.x$par[3],
                                                             error = function(x) NA )),
                       value = purrr::map_dbl(opt, ~tryCatch(.x$value,
                                                             error = function(x) NA ))
  )
  est <- dplyr::select(est, -opt)
  est2 <- est[which.min(est$value), ]
  return(est2)
  }
  if(type == 2) {
    if(is.numeric(s)) {
      start <- start |> tidyr::crossing(
        s = s
      )
    } else if(is.numeric(r) & is.numeric(ptail) & length(ptail) == 1) {
      start <- start |> tidyr::crossing(
        s = r * ptail ^ {-1 / t}
      )
    } else {
      start <- start |> tidyr::crossing(
        s = initial[2]
      )
    }

    # para_orig: parameter (tail, ei, scale)
    # para: a transformed parameter so that the optimation is unconstrained
    #para_orig <- function(para) {c()}
    est <- dplyr::mutate(start,
                         opt = purrr::pmap(list(t, e, s), function(t., e., s.) {
                           tryCatch(
                             stats::optim( par = c(t., e., s.), fn = function(x) {
                               do.call(
                                 distance_fct, list(WW = WW, tail = x[1], ei = x[2], scale = x[3])
                               )
                             },
                             lower = c(a_tail, a_ei, 0), upper = c(1, 1, Inf),
                             method = method),
                             error = function(x) NA
                           )
                         }),
                         beta = purrr::map_dbl(opt, ~tryCatch(.x$par[1],
                                                              error = function(x) NA )),
                         theta = purrr::map_dbl(opt, ~tryCatch(.x$par[2],
                                                               error = function(x) NA )),
                         rho = purrr::map_dbl(opt, ~tryCatch(.x$par[3] * ptail ^ {1 / .x$par[1]},
                                                             error = function(x) NA )),
                         sigma =  purrr::map_dbl(opt, ~tryCatch(.x$par[3],
                                                                error = function(x) NA )),
                         value = purrr::map_dbl(opt, ~tryCatch(.x$value,
                                                               error = function(x) NA ))
    )
    est <- dplyr::select(est, -opt)
    est2 <- est[which.min(est$value), ]
    return(est2)
  }
}

# WW - vector of inter-exceedance times
# ptail - P(JJ > u) = k/n (number of exceedances over number of observations)
# distance_fct - r-function calculating the distance between the asymptotic and empirical cdf
# start in (0,1)^2
#' @export
optim_singlestart <- function(WW, ptail, distance_fct, start = NULL,
                              a_tail = 0.1, a_ei = 0.1, method = "L-BFGS-B", ...) {
  if(is.null(start)) {
    start <- c((1 + a_tail) / 2, (1 + a_ei) / 2)
  }
  est <- tryCatch(
    stats::optim(par = start, fn = function(x) {
      do.call(
        distance_fct, list(WW = WW, tail = x[1], ei = x[2], ptail = ptail, ...)
        )
      },
      lower = c(a_tail, a_ei), upper = c(1, 1), method = method)$par,
    error = function(x) NA
    )
  return(est)
  }
