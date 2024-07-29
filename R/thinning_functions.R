#' Thinning function
#'
#' This function returns all events that
#' exceed a specified threshold \code{u} from a tibble, a data.frame or a
#' two-column matrix \code{data}.
#' Alternatively, the function can also return the \code{k} largest observations.
#'
#' @param data data.frame, tibble or matrix with two columns.
#' The first column contains the event magnitudes (class "numeric") and
#' in the second column contains teh corresponding waiting times (class "numeric")
#' where the i-th waiting time is the time between the (i-1)-th
#' and i-th event
#' @param k number of exceedances
#' @param u threshold
#'
#' @return
#' A tibble with two columns:
#' -  the exceedances \code{newJJ} and
#' -  the corresponding waiting times \code{newWW}
#' @export
#'
#' @examples
#' data_generation(n = 10000, stability = 0.9, ei = 0.7, wait_dist = "ML") |>
#'   thin(k = 100)
#' data_generation(n = 10000, stability = 0.9, ei = 0.7, wait_dist = "ML") |>
#'   thin(u = 100)
#'



thin <- function(data, k = NULL, u = NULL) {
  # input control
  if(length(dim(data)) != 2) {
    stop("'data' should be a matrix, tibble or data.frame with two columns")
  }

  data <- as.matrix(data)
  JJ <- data[, 1]
  WW <- data[, 2]

  if(anyNA(WW)) {
    JJ <- JJ[!is.na(WW)]
    WW <- WW[!is.na(WW)]
    warning("The second column contains missing values (NA). These entries have been deleted.")
  }
  if(any(WW <= 0)) {
    stop("The second column of 'data' should contain the waiting times being
         greater than zero.")
  }
  if(anyNA(JJ)) {
    ind <- which(is.na(JJ))
    for(i in seq_along(ind)) {
      WW[ind[i] + 1] <- WW[ind[i] + 1] + WW[ind[i]]
    }
    WW <- WW[-ind]
    JJ <- JJ[-ind]
    warning("The first column contains missing values (NA). These entries have been deleted
            and the corresponding waiting times has been added to the next observation.")
  }
  n <- length(JJ)

  if(is.null(k) & is.null(u)) {
    stop("Input of arguments k and u are missing. You have to set exactly one
         of both.")
  } else if(!is.null(k) & !is.null(u)) {
    k_1 <- sum(JJ > u)
    if(k_1 != k) {
      warning("Number of exceedances k = ", k, " and given threshold u =", u,
              "do not match. For computation, only the input of k = ", k,
              " will be used.")
    }
  } else if(is.null(k)) {
    k <- sum(JJ > u)
    if( k < 2) {
      stop("The threshold u = ", u, " is too high. There are less
                      than two magnitudes that exceed u")
    }
  }
  if(k > n) {
    stop("Can't threshold to ", k, " exceedances if I only have ", n, " observations.")
  }
  k0 <- k
  while(sort(JJ, decreasing = T)[k] == sort(JJ, decreasing = T)[k+1]) {
    k <- k + 1
  }
  if(k0 != k) {
    warning("The k-th and (k + 1)-st greatest magnitudes have the same value.
            All magnitudes equal or greater than the k-th larges magnitude are kept.
            Thus, the number of exceedances does not match k" )
  }
  idxJ <- sort(order(JJ, decreasing = TRUE)[1:k]) # Index of the k highest events
  b <- rep(0, times = n)
  b[idxJ + 1] <- 1
  b <- b[1:n]
  # b a dichotom vector of length n with 1 located one entry after an
  # exceedance occurs, otherwise 0.
  a <- 1 + cumsum(b == 1)
  newJJ <- JJ[idxJ]
  firstJJ <- newJJ[1]
  newWW <- stats::aggregate(WW, list(a), sum)$x[1:k]
  out <- tibble::tibble(newJJ = newJJ, newWW = newWW)
  return(out)
}

#' Arrivaltimes
#'
#' This function transforms a tibble, a data.frame or a
#' two-column matrix of event magnitudes and waiting times to a vector that
#' contains the cumulative waiting times till the i-th event.
#'
#' @param data data.frame, tibble or matrix with two columns.
#' The first column contains the event magnitudes (class "numeric") and
#' in the second column contains the corresponding waiting times (class "numeric")
#' where the i-th waiting time is the time between the (i-1)-th
#' and i-th event
#'
#' @return A vector that contains the arrivaltimes.
#' @export
#'
#' @examples
#' dat <- data_generation(n = 10000, stability = 0.9, ei = 0.7, wait_dist = "ML") |>
#'   thin(k = 100)
#' arrivaltime(dat)
#'

arrivaltime <- function(data) {
  # input control
  if(length(dim(data)) != 2) {
    stop("'data' should be a matrix, tibble or data.frame with two columns")
  }
  data <- as.matrix(data)
  WW <- data[, 2]
  TT <- cumsum(WW)
  return(TT)
}

#' magnitudes
#'
#' This function extracts the event magnitudes of a tibble, a data.frame or a
#' two-column matrix of event magnitudes and waiting times.
#'
#' @param data data.frame, tibble or matrix with two columns.
#' The first column contains the event magnitudes (class "numeric") and
#' in the second column contains the corresponding waiting times (class "numeric")
#' where the i-th waiting time is the time between the (i-1)-th
#' and i-th event
#'
#' @return A vector that contains the magnitudes.
#' @export
#'
#' @examples
#' dat <- data_generation(n = 10000, stability = 0.9, ei = 0.7, wait_dist = "ML") |>
#'   thin(k = 100)
#' magnitudes(dat)
#'

magnitudes <- function(data) {
  # input control
  if(length(dim(data)) != 2) {
    stop("'data' should be a matrix, tibble or data.frame with two columns")
  }
  data <- as.matrix(data)
  JJ <- data[, 1]
  return(JJ)
}

#' Interarrivaltime
#'
#' This function extracts the waiting times / interarrivaltimes of a tibble, a data.frame or a
#' two-column matrix of event magnitudes and waiting times.
#'
#' @param data data.frame, tibble or matrix with two columns.
#' The first column contains the event magnitudes (class "numeric") and
#' in the second column contains the corresponding waiting times (class "numeric")
#' where the i-th waiting time is the time between the (i-1)-th
#' and i-th event
#' @param skip_first logical; if TRUE, the first entry will be deleted.
#'
#' @return A vector that contains the interarrivaltimes.
#' @export
#'
#' @examples
#' dat <- data_generation(n = 10000, stability = 0.9, ei = 0.7, wait_dist = "ML") |>
#'   thin(k = 100)
#' interarrivaltime(dat)
#'

interarrivaltime <- function(data, skip_first = T) {
  # input control
  if(length(dim(data)) != 2) {
    stop("'data' should be a matrix, tibble or data.frame with two columns")
  }
  data <- as.matrix(data)
  WW <- data[, 2]
  if(skip_first == TRUE) {
    WW <- WW[-1] # first time excluded
  }
  return(WW)
}


