#' Data Generating Function
#'
#' Function that generates a marked point process \eqn{(JJ[i], WW[i]), i = 1,...,n}.
#' The JJ[i] is the i-th mark and WW[i] is the waiting time between the (i-1)-th and
#' i-th event. The marks (JJ[i]) form a MAR-process which is a stationary time series
#' with extremal index \eqn{\theta \in (0,1]}.
#' The generated process fulfills all assumptions so that the inter-exceedance times (IETs)
#' are asymptotically mixture distributed with the dirac measure at point zero and the
#' Mittag-Leffler distribution as parts.
#' The IETs are the times between two consequtive extremes which are identified
#' by the peak-over-threshold (POT) method.
#'
#' @param n number of observations
#' @param ei extremal index \eqn{\theta \in (0,1]}
#' @param stability stability parameter \eqn{\alpha > 0}
#' used for the waiting time distribution (see 'Details')
#' @param scale scale parameter \eqn{\sigma > 0}. Default \code{scale = 1}.
#' @param wait_dist distribution of the waiting times WW (see 'Details')
#'
#' @details
#' The event magnitudes/ marks (JJ) form a Max-Autoregression(MAR)-process.
#' This process is stationary with extremal index \eqn{\theta}.
#'
#' The waiting times (WWs) are i.i.d. and stochastically independent to the
#' event magnitudes (JJ).
#' The waiting time distribution can be chosen as "\code{stable}", "\code{ML}" or
#' "\code{pareto}".
#'
#' If \code{wait_dist = "stable"}, the \code{stability} parameter \eqn{\alpha} has to be in
#' the interval (0,1].
#' In case of \eqn{\alpha < 1}, the WWs are stable distributed with stability parameter \eqn{\alpha}.
#' Then, the distribution is heavy-tailed in the sense that mean and variance do
#' not exist.
#' In this case \code{WW} are generated using the package \code{stabledist}.
#' In the special case of \code{stability = 1}, the waiting times are equidistant equal \eqn{\sigma},
#' since the stable distribution converges in distribution to the degenerative dirac measure in \eqn{\sigma}.
#'
#' If \code{wait_dist = "ML"}, the \code{stability} parameter \eqn{\alpha} has to be in
#' the interval (0,1].
#' The WWs are Mittag-Leffler distributed with stability/tail parameter \eqn{\alpha}.
#' For \eqn{\alpha < 1} the distribution is heavy-tailed in the sense that mean
#' and variance do not exist.
#' In case of \eqn{\alpha = 1}, the Mittag-Leffler distribution corresponds to
#' an exponential distribution with rate \eqn{1/\sigma}.
#' For the data generation of the Mittag-Leffler distribution
#' the package \code{MittagLeffleR} is used.
#'
#' If \code{wait_dist = "pareto"}, the \code{stability} parameter \eqn{\alpha}
#' has to be in \eqn{(0, \infty)/{1}}.
#' If \eqn{\alpha < 1}, a Pareto distribution with shape parameter \eqn{\alpha}
#' and scale parameter \eqn{(1/(\Gamma(1-\alpha))^{1/\alpha} * \sigma} is used.
#' Then the distribution is heavy-tailed in the sense that mean and variance do not exist.
#' If \eqn{\alpha > 1}, the scale parameter becomes \eqn{(\alpha-1)/\alpha * \sigma}.
#' For \eqn{\alpha > 1}, the expected value exists and equals \eqn{\sigma}. For \eqn{\alpha > 2}
#' the variance exists additionally.
#' For the generation of the Pareto distribution the Package \code{ReIns} is used.
#'
#' For all three waiting time distributions with \code{stability} parameter \eqn{\alpha < 1}
#' it holds that they are in the domain of attraction of the stable distribution
#' which is a positively skewed sum-stable distribution.
#'
#' The generated process (JJ,WW) fulfills all assumptions so that the inter-exceedance times (IETs)
#' are asymptotically mixture distributed with the dirac measure at point zero and the
#' Mittag-Leffler distribution as parts.
#' The IETs are the times between two consequtive extremes which are identified
#' by the peak-over-threshold (POT) method.
#'
#' @return
#' A tibble with two columns
#'
#' \code{JJ} a stationary MAR-process with extremal index \code{ei}
#'
#' \code{WW} independent waiting times with stability parameter \code{stability}
#' and scale parameter \code{scale}
#'
#' @export
#'
#' @examples
#' dat <- data_generation(20, 0.5, 0.5)
#' dat
#' dat2 <- data_generation(15, 0.2, 0.85, wait_dist = "pareto")
#' dat2


data_generation <- function(n, ei = 1, stability = 1, scale = 1,
                            wait_dist = "stable") {
  ## input control:
  # 'n' number of observations
  if(!isTRUE(all.equal(n , round(n))) || length(n) != 1)
    stop("sample size has to be an integer")
  # 'ei'/ true parameter theta in (0,1]
  if(ei <= 0 || ei > 1 || length(ei) != 1)
    stop("Extremal index ei should be a single value in (0,1].")
  # 'wait_dist' distribution of waiting times
  if(! (wait_dist %in% c("stable", "ML", "pareto")) ||
     length(wait_dist) != 1 )
    stop("wait_dist should be one of following characters: stable, ML, pareto")
  # 'stability'/ true stability parameter alpha > 0
  if(stability <= 0 || length(stability) != 1)
    stop("Stability parameter should be a single positiv value.")
  # stability larger 1 only for Pareto
  if(!(wait_dist == "pareto") & stability > 1)
    stop("A stability parameter larger than 1 is only implemented within the use
         of Pareto distributed waiting times")
  if(wait_dist == "pareto" & stability == 1)
    stop("The special case of pareto waiting times with stability parameter equals
         one is unfortunately not implemented, yet.
         Please, choose a stability parameter smaller or larger than one")

  if(stability <= 1) {
    tail <- stability
  } else {
    tail <- 1
  }

  # generating event values (magnitudes):
    J <- rep(0, n)
    eps <- extRemes::revd(n, scale = 1, shape = 1, loc = 1)
    J[1] <- eps[1]
    if(n >= 2) {
      for (i in 2:n){
        J[i] <- max((1 - ei) * J[i - 1], ei * eps[i])
        }
      }
  # generating waiting times:
  # tail == 1:
if(wait_dist == "stable") { # stable

  if(stability < 1){
    # stability < 1
    sigma <- (cos(pi * tail / 2)) ^ (1 / tail)
    W <- stabledist::rstable(n, alpha = tail, beta = 1 , gamma = sigma,
                             delta = 0, pm = 1)
  } else if(stability == 1){
    W <- rep(1, n)
  }

  }
  if(wait_dist == "pareto") { # pareto

    if(stability < 1){
      tail <- stability
      sigma <- (1 / gamma(1 - tail)) ^ (1 / tail) # stability < 1
    }
    if(stability > 1){
      sigma <- (stability - 1) / stability # stability > 1
    }
    W <- ReIns::rpareto(n, shape = stability, scale = sigma)
  }


    if(wait_dist == "ML") { # Mittag-Leffler distribution
     # stability = tail < 1
       sigma <- 1

     if(stability < 1){
       tail <- stability
       W <- MittagLeffleR::rml(n, tail = tail, scale = sigma )
     }
     if(stability == 1){
       W <- stats::rexp(n, rate = 1)
     }
  }
  W <- scale * W
  return(tibble::tibble(JJ = J, WW = W))
}
