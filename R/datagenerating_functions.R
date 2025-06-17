#' Data Generating Function
#'
#' A tibble / data.frame with two columns \code{JJ} and \code{WW} where \code{JJ} are
#' the marks / magnitudes and \code{WW} are the waiting times between the (i-1)-th and
#' i-th event.
#' The generated process fulfills all assumptions so that the interexceedance times (IETs)
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
#' @param mag_dist process of the magnitudes \code{JJ} (default \code{"MAR"}).
#' The process can be chosen as \code{"MAR"} or "\code{MM}". (see 'Details')
#' @param wait_dist distribution of the waiting times \code{WW} (default \code{"stable"}).
#' The waiting time distribution can be chosen as "\code{stable}", "\code{ML}",
#' or "\code{pareto}". (see 'Details')
#' @param u threshold  (default NULL); if \eqn{u} is numeric, it holds \code{JJ[1] > u}
#' and the length of the data.frame is \code{n + 1}
#' @param alpha is a auxility vector for the moving maxima process (default NULL).
#' It is only needed if \code{mag_dist = "MM"}.
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
#' The marks \code{JJ} either form a max-autoregressive process (MAR) which is a
#' stationary time series with extremal index \eqn{\theta \in (0,1]}.
#' This process is stationary with extremal index \eqn{\theta}.
#'
#' Alternatively, the marks \code{JJ} can form a moving maximum process (MM) which
#' is also stationary time series with extremal index \eqn{\theta \in (0,1]}.
#' A auxility parameter \code{alpha} must be set such that \code{sum(alpha)} = 1
#' and \code{max(alpha)} = \eqn{\theta}. If \code{alpha = NULL}, the auxility parameter is chosen for you
#' (see \code{FCPP:::MMprocess})
#'
#' The waiting times \code{WW} are i.i.d. and stochastically independent to the
#' marks \code{JJ}.
#' The waiting time distribution can be chosen as "\code{stable}", "\code{ML}",
#' or "\code{pareto}".
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
#' For all three waiting time distributions with \code{stability} parameter \eqn{\alpha < 1}
#' it holds that they are in the domain of attraction of the stable distribution
#' which is a positively skewed sum-stable distribution.
#'
#' The generated process \code{(JJ,WW)} is a marked point process and fulfills all assumptions so that the interexceedance times (IETs)
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

data_generation <- function(n, stability = 1, ei = 1, scale0 = 1, mag_dist = "MAR",
                            wait_dist = "stable", u = NULL, alpha = NULL) {
  ## input control:
  # 'n' number of observations
  if(!isTRUE(all.equal(n , round(n))) || length(n) != 1)
    stop("The sample size has to be an integer.")
  # 'ei'/ true parameter theta in (0,1]
  if(ei <= 0 || ei > 1 || length(ei) != 1)
    stop("Extremal index ei should be a single value in (0,1].")
  # 'wait_dist' distribution of waiting times
  if(! (wait_dist %in% c("stable", "ML", "pareto")) ||
     length(wait_dist) != 1 )
    stop("wait_dist should be one of following characters: stable, ML, pareto.")
  # 'stability'/ true stability parameter alpha > 0
  if(stability <= 0 || length(stability) != 1)
    stop("Stability parameter should be a single positiv value.")
  # stability larger 1 only for Pareto
  if(!(wait_dist == "pareto") & stability > 1)
    stop("A stability parameter larger than 1 is only implemented within the use
         of Pareto distributed waiting times.")
  if(wait_dist == "pareto" & stability == 1)
    stop("The special case of pareto waiting times with stability parameter equals
         one is unfortunately not implemented, yet.
         Please, choose a stability parameter smaller or larger than one.")
  if(mag_dist == "MAR" & is.numeric(alpha))
    warning("For the MAR-process no alpha is needed")
  if(mag_dist == "MM" & is.numeric(alpha)) {
    if(sum(alpha) != 1) {
      alpha <- NULL
      warning("The vector alpha has to sum to one. Therefore the standard choice is used.")
    }
  } else {u <- NULL}

  if(stability <= 1) {
    tail <- stability
  } else {
    tail <- 1
  }
  if(mag_dist == "MAR") {
    J <- MARprocess(n, ei, u)
  }
  if(mag_dist == "MM") {
    J <- MMprocess(n, ei, u)
  }
  # generating waiting times:
  # tail == 1:
  m <- length(J)
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

MMprocess <- function(n , ei, u = NULL) {
  if(ei == 1) {alpha <- 1}
  if(ei < 1 & ei >= 0.8) {alpha <- c(ei, 1-ei)}
  if(ei < 0.8 & ei >= 0.6){alpha <- c(ei, 0.2, 1-ei-0.2)}
  if(ei < 0.6 & ei >= 0.5){alpha <- c(ei, 0.2, 0.2, 1-ei-0.4)}
  if(ei < 0.5 & ei >= 0.4 ){alpha <- c(ei, 0.2, 0.2, 0.1, 1-ei-0.5)}
  if(ei < 0.4 & ei >= 0.3 ){alpha <- c(ei, 0.2, 0.2, 0.1, 0.1, 1-ei-0.6)}
  if(ei < 0.3 & ei >= 0.2 ){alpha <- c(ei, 0.2, 0.2, 0.1, 0.1, 0.1, 1-ei-0.7)}
  if(ei < 0.2 & ei >= 0.1 ){alpha <- c(ei, rep(0.1, 5), rep(0.05, 7), 1-ei-0.85)}
  k <- length(alpha)
  if(is.numeric(u)) {
    p <-  extRemes::pevd(u / ei, scale = 1, shape = 1, loc = 1)
    r <-  stats::runif(1, p, 1)
    J0 <- extRemes::qevd(r, scale = 1, shape = 1, loc = 1)
    eps <- extRemes::revd(n + k + 1, scale = 1, shape = 1, loc = 1)
    eps[k] <- J0
  }
  m <- length(eps)
  J <- rep(0, m)
  for(i in seq_along(J)) {
    J[i] <- max(alpha * eps[(i+k-1):i])
  }
  return(J)
}

MARprocess <- function(n, ei, u) {
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
  return(J)
}
