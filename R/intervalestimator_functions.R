# modified Intervals estimator of Ferro and Segers (2003)

interval_estimator <- function(WW){
  if(is.data.frame(WW) & length(WW) == 1) {
    WW <- dplyr::pull(WW)
  }
  kstar <- length(WW)
  if(any(WW > 2)) {
    WW1 <- WW2 <- WW
    WW1[WW <= 1] <- 1
    WW2[WW <= 2] <- 2
    est <- (2 * sum(WW1 - 1) ^ 2) / (kstar * sum((WW1 - 1) * (WW2 - 2)))
  }
  else{
    est <- (2 * sum(WW) ^ 2) / (kstar * sum(WW ^ 2))
  }
  est <- min(est, 1)
  return(est)
}
