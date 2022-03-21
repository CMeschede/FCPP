#' Thinning function
#'
#' @param data data.frame, tibble or matrix with two columns.
#' In the first column are the event magnitudes (JJ) stored and
#' in the second column are the corresponding waiting times (WW)
#' where the i-th waiting time (WW) is the time between the (i-1)-th
#' and i-th event
#' @param k number of exceedances
#' @param u threshold
#'
#' @return
#' @export
#'


thin <- function(data, k = NULL, u = NULL) {
  # input control
  if(length(dim(data)) != 2) {
    stop("'data' should be a matrix, tibble or data.frame with two columns")
  }
  if(dim(data)[2] != 2) {
    stop("'data' should be a matrix, tibble or data.frame with two columns ")
  }
  data <- as.matrix(data)
  JJ <- data[, 1]
  WW <- data[, 2]
  n <- length(JJ)
  if(any(WW <= 0)) {
    stop("the second column of 'data' should contain the waiting times being
         greater than zero")
  }
  if(is.null(k) & is.null(u)) {
    stop("Input of arguments k and u are missing. You have to set at least one
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
  idxJ <- sort(order(JJ, decreasing = TRUE)[1:k]) #Index of the k highest events
  b <- rep(0, times = n)
  b[idxJ + 1] <- 1
  b <- b[1:n]
  # b a dichotom vector of length n with 1 located one entry after an
  # exceedance occurs, otherwise 0.
  a <- 1 + cumsum(b == 1)
  newJJ <- JJ[idxJ]
  firstJJ <- newJJ[1]
  newWW <- aggregate(WW, list(a), sum)$x[1:k]
  out <- tibble::tibble(newJJ = newJJ, newWW = newWW)
  return(out)
}

#' Arrivaltime
#'
#' @param data dataframe with two columns
#'
#' @return
#' @export
#'

arrivaltime <- function(data) {
  # input control
  if(length(dim(data)) != 2) {
    stop("'data' should be a matrix, tibble or data.frame with two columns")
  }
  if(dim(data)[2] != 2) {
    stop("'data' should be a matrix, tibble or data.frame with two columns ")
  }
  data <- as.matrix(data)
  JJ <- data[, 1]
  WW <- data[, 2]
  TT <- cumsum(WW)
  return(TT)
}

#' magnitudes
#'
#' @param data dataframe with two columns
#'
#' @return
#' @export
#'

magnitudes <- function(data) {
  # input control
  if(length(dim(data)) != 2) {
    stop("'data' should be a matrix, tibble or data.frame with two columns")
  }
  if(dim(data)[2] != 2) {
    stop("'data' should be a matrix, tibble or data.frame with two columns ")
  }
  data <- as.matrix(data)
  JJ <- data[, 1]
  WW <- data[, 2]
  return(JJ)
}

#' Interarrivaltime
#'
#' @param data dataframe with two columns
#'
#' @return
#' @export
#'

interarrivaltime <- function(data) {
  # input control
  if(length(dim(data)) != 2) {
    stop("'data' should be a matrix, tibble or data.frame with two columns")
  }
  if(dim(data)[2] != 2) {
    stop("'data' should be a matrix, tibble or data.frame with two columns")
  }
  data <- as.matrix(data)
  JJ <- data[, 1]
  WW <- data[, 2]
  k <- length(WW)
  WW <- WW[-1]
  return(WW)
}
