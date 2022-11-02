#' Data Generating Function
#'
#' Function that generates a marked point process (JJ[i], WW[i]) \eqn{i = 1,...,n}.
#' The JJ[i] is the i-th mark and WW[i] is the waiting time between the i-1- and
#' i-th event.
#'
#' @param n number of observations
#' @param ei true parameter theta in (0,1] (extremal index)
#' @param stability The stability parameter of the distribution of the waiting times.
#' The stability parameter is in \eqn{\mathbb{R}^+\setminus\lbrace 1\rbrace} if \code{wait_dist}
#' is "pareto" and in (0,1] else.
#' @param wait_dist distribution of waiting times. Should be "stable",
#' "ML" (Mittag-Leffler Distribution) or "pareto".
#'
#' @details
#' If the distribution of the waiting times is chosen as "stable" and the stability parameter
#' is smaller than 1, then the waiting times are stable distributed.
#' In this case \code{WW} are generated using the
#' package \code{stabledist}. In the special case of \code{stability}=1, the stable distribution
#' simplifies to a dirac measure in 1.
#'
#' For the case of \code{wait_dist="ML"} a Mittag-Leffler distribution where the
#' tail parameter and the stability parameter are the same. If the stability parameter
#' is equal to 1, then the Mittag-Leffler distribution corresponds to an exponential
#' distribution with
#' rate 1. For the data generation of the Mittag-Leffler distribution
#' the package \code{MittagLeffleR} is used.
#'
#' For \code{wait_dist="pareto"} and a stability parameter
#' \eqn{\alpha} smaller than 1, a Pareto distribution
#'  with shape parameter \eqn{\alpha} and scale parameter
#'  \eqn{\big(1/\Gamma(1-\alpha)\big)^{1/\alpha}} is used
#'  to generate the waiting times. For \eqn{\alpha} larger than 1 the
#'  scale parameter becomes \eqn{(\alpha-1)/\alpha}.For the
#'  generation of the waiting times the Package \code{ReIns} is used.
#'
#'
#' @return
#' A tibble with
#'
#' \code{JJ} a stationary MAR-process with extremal index \eqn{\theta} (ei).
#'
#' \code{WW} independent waiting times between the events.
#'
#' @export
#'
#' @examples
#' dat <- data_generation(20, 0.5, 0.5)
#' dat
#' dat2 <- data_generation(15, 0.2, 0.85, wait_dist = "pareto")
#' dat2


data_generation <- function(n, ei = 1, stability = 1,
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


    if(wait_dist == "ML") { # Mittag-Leffler
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
  return(tibble::tibble(JJ = J, WW = W))
}
