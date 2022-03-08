thin <- function(data, k = NULL, u = NULL, type = c("mod1", "mod2")) {
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
  if(any(WW <= 0)){
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
    if( k < 2) { stop("The threshold u = ", u, " is too high. There are less
                      than two magnitudes that exceed u") }
  }
  if (k > n) {
    stop("Can't threshold to ", k, " exceedances if I only have ", n,
         " observations.")
  }
  idxJ <- sort(order(JJ, decreasing = TRUE)[1:k]) #Index of the k highest events
  b <- rep(0, times = n)
  if(type == "mod1") {
    b[idxJ + 1] <- 1
    b <- b[1:n]
  } else if(type == "mod2") {
    b[idxJ] <- 1
    #'b' a dichotom vector of length n with 1 where exceedances are located
    # otherwise 0
  }
  a <- 1 + cumsum(b == 1)
  newJJ <- JJ[idxJ]
  firstJJ <- newJJ[1]
  # newJJ <- newJJ[-1]
  if(type = "mod1") {
    newWW <- aggregate(WW, list(a), sum)$x[1:k]
    # newWW <- aggregate(WW, list(a), sum)$x[2:k]
  } else if(type = "mod2") {
    newWW <- aggregate(WW, list(a), sum)$x
    if(length(newWW) == k) {
      newWW <- newWW[1:k]
      # newWW <- newWW[1:(k-1)]
    } else if(length(newWW == (k + 1))) {
      newWW <- newWW[2:(k + 1)]
      # newWW <- newWW[2:k]
    }
  }
  out <- tibble(newJJ = newJJ, newWW = newWW)
  return(out)
}

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

interarrivaltime <- function(data, type) {
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
  if(type == "mod1") {
    WW <- WW[-1]
  } else if(type == "mod2") {
    WW <- WW[-k]
  }
  return(WW)
}
