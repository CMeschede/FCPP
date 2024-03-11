#' Data Generating Function
#'
#' A tibble / data.frame with two columns \code{JJ} and \code{WW} where \code{JJ} are
#' the marks / magnitudes and \code{WW} are the waiting times between the (i-1)-th and
#' i-th event.
#' The generated process fulfills all assumptions so that the inter-exceedance times (IETs)
#' regarding a threshold \eqn{u} are asymptotically mixture distributed with the
#' dirac measure at point zero and the Mittag-Leffler distribution as parts.
#' The IETs are the times between two consequtive extremes which are identified
#' by the peak-over-threshold (POT) method.
#'
#' @param n number of observations
#' @param ei extremal index \eqn{\theta \in (0,1]}
#' @param stability stability parameter \eqn{\alpha > 0}
#' used for the waiting time distribution (see 'Details')
#' @param scale0 scale parameter \eqn{\rho > 0} of the waiting time distribution.
#' Default \code{scale0 = 1}
#' @param wait_dist distribution of the waiting times \code{WW}.
#' The waiting time distribution can be chosen as "\code{stable}", "\code{ML}",
#' "\code{pareto}" or "\code{shifted_pareto}". (see 'Details')
#' @param u threshold  (default NULL); if \eqn{u} is numeric, it holds \code{JJ[1] > u}
#' and the length of the data.frame is \code{n + 1}
#'
#' @return
#' A tibble with two columns
#'
#' -  \code{JJ} a stationary MAR-process with extremal index \code{ei}
#'
#' -  \code{WW} independent waiting times with stability parameter \code{stability}
#' and scale parameter \code{scale0}
#'
#' @details
#' The marks \code{JJ} form a max-autoregressive process (MAR) which is a
#' stationary time series with extremal index \eqn{\theta \in (0,1]}.
#' This process is stationary with extremal index \eqn{\theta}.
#'
#' The waiting times \code{WW} are i.i.d. and stochastically independent to the
#' marks \code{JJ}.
#' The waiting time distribution can be chosen as "\code{stable}", "\code{ML}",
#' "\code{pareto}" or "\code{shifted_pareto}".
#'
#' If \code{wait_dist = "stable"}, the \code{stability} parameter \eqn{\alpha}
#' has to be in the interval (0,1].
#' In case of \eqn{\alpha < 1}, the WWs are stable distributed with stability parameter \eqn{\alpha}.
#' Then, the distribution is heavy-tailed in the sense that mean and variance do
#' not exist.
#' In this case \code{WW} are generated using the package \code{stabledist}.
#' In the special case of \code{stability = 1}, the waiting times are equidistant equal \eqn{\rho},
#' since the stable distribution converges in distribution to the degenerative dirac measure in \eqn{\rho}.
#'
#' If \code{wait_dist = "ML"}, the \code{stability} parameter \eqn{\alpha} has to be in
#' the interval (0,1].
#' The WWs are Mittag-Leffler distributed with stability/tail parameter \eqn{\alpha}.
#' For \eqn{\alpha < 1} the distribution is heavy-tailed in the sense that mean
#' and variance do not exist.
#' In case of \eqn{\alpha = 1}, the Mittag-Leffler distribution corresponds to
#' an exponential distribution with rate \eqn{1/\rho}.
#' For the data generation of the Mittag-Leffler distribution
#' the package \code{MittagLeffleR} is used.
#'
#' If \code{wait_dist = "pareto"}, the \code{stability} parameter \eqn{\alpha}
#' has to be in \eqn{(0, \infty)/ \{1 \}}.
#' If \eqn{\alpha < 1}, a Pareto distribution with shape parameter \eqn{\alpha}
#' and scale parameter \eqn{(1/(\Gamma(1-\alpha))^{1/\alpha} \cdot \rho} is used.
#' Then the distribution is heavy-tailed in the sense that mean and variance do not exist.
#' If \eqn{\alpha > 1}, the scale parameter becomes \eqn{(\alpha-1)/\alpha \cdot \rho}.
#' For \eqn{\alpha > 1}, the expected value exists and equals \eqn{\rho}. For \eqn{\alpha > 2}
#' the variance exists additionally.
#' For the generation of the Pareto distribution the Package \code{ReIns} is used.
#'
#' If \code{wait_dist = "shifted_pareto"}, the \code{stability} parameter \eqn{\alpha}
#' has to be in \eqn{(0, 1]}.
#' If \eqn{\alpha < 1}, the waiting times are as described as for \code{wait_dist = "pareto"},
#' but are additively shifted by one. That shift leads to the special case of \code{stability = 1},
#' where the waiting times are equidistant equal \eqn{\rho},
#' since the shifted Pareto distribution converges in distribution to the
#' degenerative dirac measure in \eqn{\rho}
#' while the normal/unshifted Pareto distribution converges to the degenerative
#' dirac measure in zero as \eqn{\alpha \to 1}.
#'
#' For all four waiting time distributions with \code{stability} parameter \eqn{\alpha < 1}
#' it holds that they are in the domain of attraction of the stable distribution
#' which is a positively skewed sum-stable distribution.
#'
#' The generated process \code{(JJ,WW)} is a marked point process and fulfills all assumptions so that the inter-exceedance times (IETs)
#' are asymptotically mixture distributed with the dirac measure at point zero and the
#' Mittag-Leffler distribution as parts.
#'
#'
#' @export
#'
#' @examples
#' dat <- data_generation(n = 1000, stability = 0.9, ei = 0.8)
#' dat
#' dat2 <- data_generation(n = 200, stability = 1, ei = 0.7, scale0 = 2, wait_dist = "ML")
#' dat2


data_generation <- function(n, stability = 1, ei = 1, scale0 = 1,
                            wait_dist = "stable", u = NULL) {
  ## input control:
  # 'n' number of observations
  if(!isTRUE(all.equal(n , round(n))) || length(n) != 1)
    stop("The sample size has to be an integer.")
  # 'ei'/ true parameter theta in (0,1]
  if(ei <= 0 || ei > 1 || length(ei) != 1)
    stop("Extremal index ei should be a single value in (0,1].")
  # 'wait_dist' distribution of waiting times
  if(! (wait_dist %in% c("stable", "ML", "pareto", "shifted_pareto")) ||
     length(wait_dist) != 1 )
    stop("wait_dist should be one of following characters: stable, ML, pareto, shifted_pareto.")
  # 'stability'/ true stability parameter alpha > 0
  if(stability <= 0 || length(stability) != 1)
    stop("Stability parameter should be a single positiv value.")
  # stability larger 1 only for Pareto
  if(!(wait_dist == "pareto") & stability > 1)
    stop("A stability parameter larger than 1 is only implemented within the use
         of Pareto distributed waiting times (without a shift).")
  if(wait_dist == "pareto" & stability == 1)
    stop("The special case of pareto waiting times with stability parameter equals
         one is unfortunately not implemented, yet.
         Please, choose a stability parameter smaller or larger than one.")
  if(stability <= 1) {
    tail <- stability
  } else {
    tail <- 1
  }

  # generating event values (magnitudes):
    eps <- extRemes::revd(n, scale = 1, shape = 1, loc = 1)
    # if a threshold u is given the first given magnitude exceeds the threshold
    # and the sample size becomes n + 1
    if(is.numeric(u)) {
      p <-  extRemes::pevd(u, scale = 1, shape = 1, loc = 1)
      r <-  stats::runif(1, p, 1)
      J0 <- extRemes::qevd(r, scale = 1, shape = 1, loc = 1)
      eps <- c(J0, eps)
    }
    m <- length(eps)
    J <- rep(0, m)
    J[1] <- eps[1]
    if(m >= 2) {
      for (i in 2:m){
        J[i] <- max((1 - ei) * J[i - 1], ei * eps[i])
        }
      }
  # generating waiting times:
  # tail == 1:
if(wait_dist == "stable") { # stable

  if(stability < 1){
    # stability < 1
    sigma <- (cos(pi * tail / 2)) ^ (1 / tail)
    W <- stabledist::rstable(m, alpha = tail, beta = 1 , gamma = sigma,
                             delta = 0, pm = 1)
  } else if(stability == 1){
    W <- rep(1, m)
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
    W <- ReIns::rpareto(m, shape = stability, scale = sigma)
  }

  if(wait_dist == "shifted_pareto") { # pareto

    if(stability < 1){
      tail <- stability
      sigma <- (1 / gamma(1 - tail)) ^ (1 / tail) # stability < 1
      W <- ReIns::rpareto(m, shape = stability, scale = sigma) + 1
    }
  }

    if(wait_dist == "ML") { # Mittag-Leffler distribution
     # stability = tail < 1
       sigma <- 1

     if(stability < 1){
       tail <- stability
       W <- MittagLeffleR::rml(m, tail = tail, scale = sigma )
     }
     if(stability == 1){
       W <- stats::rexp(m, rate = 1)
     }
  }
  W <- scale0 * W
  return(tibble::tibble(JJ = J, WW = W))
}
