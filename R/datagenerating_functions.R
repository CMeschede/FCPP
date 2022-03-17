#' Data Generating Function
#'
#'
#' @param n number of observations
#' @param tail true tail parameter beta in (0,1]
#' @param ei true parameter theta in (0,1] (extremal index)
#' @param wait_dist distribution of waiting times
#' @param mag_dist distribution of event magnitudes
#'
#' @return
#' @export
#'
#'
data_generation <- function(n, tail = 1, ei = 1,
                            wait_dist = "stable", mag_dist = "MAR") {
  ## input control:
  # 'n' number of observations
  if(!isTRUE(all.equal(n , round(n))) || length(n) != 1)
    stop("sample size has to be an integer")
  # 'tail'/ true parameter beta in (0,1]
  if(tail <= 0 || tail > 1 || length(tail) != 1)
    stop("tail parameter tail should be a single value in (0,1].")
  # 'ei'/ true parameter theta in (0,1]
  if(ei <= 0 || ei > 1 || length(ei) != 1)
    stop("Extremal index ei should be a single value in (0,1].")
  # 'wait_dist' distribution of waiting times
  if(! (wait_dist %in% c("stable", "ML", "pareto", "exponential", "dirac")) ||
     length(wait_dist) != 1 )
    stop("wait_dist should be one of following characters: stable, ML, pareto,
         exponential , dirac.")
  # 'wait_dist' distribution of waiting times
  if((wait_dist %in% c("exponential", "dirac") & tail < 1) ||
     (wait_dist %in% c("stable", "ML", "pareto") & tail == 1))
    stop("if tail = 1, wait_dist should be one of the following characters:
         exponential, dirac. if tail < 1, wait_dist schould be one of the
         following characters: stable , ML , pareto.")
  # 'mag_dist' distribution of events/magnitudes
  if(! (mag_dist %in% c("MAR","MM")) || length(mag_dist) != 1 )
    stop("mag_dist should be one of following characters: MAR , MM.")

  # generating event values (magnitudes):
  if(mag_dist == "MAR") {
    J <- rep(0, n)
    eps <- extRemes::revd(n, scale = 1, shape = 1, loc = 1)
    J[1] <- eps[1]
    if(n >= 2) {
      for (i in 2:n){
        J[i] <- max((1 - ei) * J[i - 1], ei * eps[i])
        }
      }
    } else if(mag_dist == "MM") {
    a <- c(ei, (1 - ei) / 4, (1 - ei) / 4, (1 - ei) / 4)
    m <- length(a)
    J <- rep(0, n)
    eps <- extRemes::revd(n + m, scale = 1, shape = 1, loc = 1)
    for(i in 1:n){
      J[i] <- max(a[1] * eps[i + m - 1], a[2] * eps[i + m - 2],
                  a[3] * eps[i + m - 3], a[4] * eps[i + m - 4])
      }
    }
  # generating waiting times:
  # tail == 1:
  if(wait_dist == "exponential") { # exponential waiting times
    sigma <- 1
    W <- rexp(n, rate = 1 / sigma)
  } else if(wait_dist == "dirac") { # equidistant waiting times
    W <- rep(1, n)
  } else if(wait_dist == "stable") { # stable
    sigma <- (cos(pi * tail / 2)) ^ (1 / tail)
    W <- stabledist::rstable(n, alpha = tail, beta = 1 , gamma = sigma, delta = 0, pm = 1)
  } else if(wait_dist == "pareto") { # pareto
    sigma <- (1 / gamma(1 - tail)) ^ (1 / tail)
    W <- ReIns::rpareto(n, shape = tail, scale = sigma)
  } else if(wait_dist == "ML") { # Mittag-Leffler
    sigma <- 1
    n <- as.numeric(n)
    W <- MittagLeffleR::rml(n, tail = tail, scale = sigma )
  }
  return(tibble::tibble(JJ = J, WW = W))
}
