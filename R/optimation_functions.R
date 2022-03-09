# WW - vector of inter-exceedance times
# n - number of observations
# distance_fct - r-function calculating the distance between the asymptotic and empirical cdf
optim_multistart <- function(WW, n, distance_fct,
                                 t = c(0.25,0.55,0.85), e = c(0.25,0.55,0.85),
                                 a_tail = 0.1, a_ei = 0.1, method = "L-BFGS-B", ...) {
  start <- tidyr::crossing(
    t = t,
    e = e
  )
  est <- dplyr::mutate(start,
                       opt = purrr::pmap(list(t, e), function(t., e.) {
                         tryCatch(
                           optim( par = c(t., e.), fn = function(x) {
                             do.call(
                               distance_fct, list(
                                 WW = WW, tail = x[1], ei = x[2], n = n, ...
                                 )
                               )
                             },
                             lower = c(a_tail, a_ei), upper = c(1, 1),
                             method = method),
                           error = function(x) NA
                           )
                         }
                       ),
                       beta = purrr::map_dbl(est , ~tryCatch(.x$par[1],
                                                             error = function(x) NA )),
                       theta = purrr::map_dbl(est , ~tryCatch(.x$par[2],
                                                              error = function(x) NA )),
                       value = purrr::map_dbl(est , ~tryCatch(.x$value,
                                                              error = function(x) NA ))
  )
  est <- dplyr::select(est, -opt)
  est2 <- est[which.min(est$value), ]
  return(est2)
}

# WW - vector of inter-exceedance times
# n - number of observations
# distance_fct - r-function calculating the distance between the asymptotic and empirical cdf
# start in (0,1)^2
optim_singlestart <- function(WW, n, distance_fct, start = c(0.55,0.55),
                              a_tail = 0.1, a_ei = 0.1, method = "L-BFGS-B", ...) {
  est <- tryCatch(
    optim( par = c(0.55,0.55), fn = function(x) {
      do.call(
        distance_fct, list(WW = WW, tail = x[1], ei = x[2], n = n, ...)
        )
      },
      lower = c(a_tail, a_ei), upper = c(1, 1), method = method)$par,
    error = function(x) NA
    )
  return(est)
  }
