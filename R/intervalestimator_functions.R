#' Modified Intervals estimator of Ferro and Segers (2003)
#'
#' The interval estimator of Ferro and Segers for estimating the extremal index \eqn{\theta}
#' using waiting times in continuous time.
#'
#' @param WW  a vector or tibble (data.frame) with one column containing the
#' positive valued sample
#'
#' @details
#' The original Interval estimator is for events that occurr in equidistant time intervals and thus
#' interexceedance times (the waiting time between two extremes) are integeres
#' and always one or larger.
#'
#' In this function we allow inter-exceedance times that are continuous in time and
#' can be lower than one.
#'
#' Original interval estimator:
#' \deqn{
#'  \hat{\theta}_n = \begin{cases}
#'  \max \left\{ 1, \, \dfrac{2 \cdot \left(\sum_{i = 1}^k W_i - 1 \right)^2}{k \cdot \sum_{i=1}^k (W_i - 1)(W_i -2)}\right\}, & \text{if  } \max\{W_i\} > 2\\
#'  \max \left\{ 1, \, \dfrac{2 \cdot \left(\sum_{i = 1}^k W_i \right)^2}{k \cdot \sum_{i=1}^k W_i^2}\right\}, & \text{if  } \max\{W_i\} \le 2 \end{cases}
#' }
#'
#' Modified interval estimator:
#'
#' \deqn{
#'  \tilde{\theta}_n = \begin{cases}
#'  \max \left\{ 1, \, \dfrac{2 \cdot \left(\sum_{i = 1}^k \max\{0, \, W_i - 1\} \right)^2}{k \cdot \sum_{i=1}^k \max\{0, \, W_i - 1\}\max\{0, \, W_i - 2\}}\right\}, & \text{if  } \max\{W_i\} > 2\\
#'  \max \left\{ 1, \, \dfrac{2 \cdot \left(\sum_{i = 1}^k W_i \right)^2}{k \cdot \sum_{i=1}^k W_i^2}\right\}, & \text{if  } \max\{W_i\} \le 2 \end{cases}
#' }
#'
#' @export
#'
#' @examples
#' # Data generation:
#' dat <- data_generation(n = 10000, stability = 1, ei = 0.7, wait_dist = "ML")
#' # Thinning the data:
#' dat_thinned <- thin(dat, k = 500)
#' # waiting times:
#' WW <- interarrivaltime(dat_thinned)
#'
#' # Optimization:
#' interval_estimator(WW = WW)
#'
interval_estimator <- function(WW){
  if(is.data.frame(WW) & length(WW) == 1) {
    WW <- dplyr::pull(WW)
  }
  kstar <- length(WW)
  if(any(WW > 2)) {
    WW1 <- WW2 <- WW
    WW1[WW <= 1] <- 1
    WW2[WW <= 2] <- 2
    est <- log(2) + 2 * log(sum(WW1 - 1)) - log(kstar) - log(sum((WW1 - 1) * (WW2 - 2)))
  }
  else{
    est <- log(2) + 2 * log(sum(WW)) - log(kstar) - log(sum(WW ^ 2))
  }
  est <- min(exp(est), 1)
  return(est)
}
