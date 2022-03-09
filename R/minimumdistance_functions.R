## minimum distance functions

distance_cm <- function(WW, tail, ei, n) {
  s <- (1 - ei) + ei * distance_cm_mod1(WW, tail, ei, n)
  return(s)
}

distance_cm_mod1 <- function(WW, tail, ei, n) {
  if(is.WW.frame(WW) & length(WW) == 1) {
    WW <- dplyr::pull(WW)
  }
  if(any(WW < 0)){
    stop("WW (waiting times) should be positive values")
  }
  k <- length(WW) + 1
  WW_trans <- WW * (k / n) ^ (1 / tail)
  pn <- (0:(k - 1)) / k
  pmls <- MittagLeffleR::pml(WW_trans, tail = tail, scale = ei ^ (-1 / tail))
  pmisch <- (1 - ei) + ei * pmls
  weights <- diff(sort(c(0, pmls, 1), decreasing = F))
  s <- sum( (((sort(c(1 - ei, pmisch)) - pn) ^ 2) * weights) )
  return(s)
}

distance_cm_mod2 <- function(WW, tail, ei, n) {
  if(is.WW.frame(WW) & length(WW) == 1) {
    WW <- dplyr::pull(WW)
  }
  if(any(WW < 0)) {
    stop("WW (waiting times) should be positive values")
  }
  k <- length(WW) + 1
  WW_trans <- WW * (k / n) ^ (1 / tail)
  pn <- (0:(k - 1)) / k
  pmls <- MittagLeffleR::pml(WW_trans, tail = tail, scale = ei ^ (-1  /tail))
  weights <- diff(sort(c(0, pmls, 1), decreasing = F))
  s <- sum( ( sort(c(0,pmls)) -
                pmax(0, pmin(1, (1 / ei * (pn - 1 + ei)))) )  ^ 2  * weights )
  return(s)
}

distance_cm_mod3 <- function(WW, tail, ei, n) {
  s <- (1 / ((0.5 - 1 / 3 * ei) * ei)) * distance_cm_mod1(WW, tail, ei, n)
  return(s)
}
