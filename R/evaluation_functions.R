## bias und rmse

rmse <- function(x, x.wahr, na.rm = T, ...) {
  if( length(x.wahr) != 1 & length(x.wahr) != length(x))
    warning("x.wahr should be of length 1 or of the same length as x")
  return(sqrt(mean((x - x.wahr) ^ 2, na.rm = na.rm, ...)))
}

bias <- function(x, x.wahr, na.rm = T, ...) {
  if( length(x.wahr) != 1 & length(x.wahr) != length(x))
    warning("x.wahr should be of length 1 or of the same length as x")
  return(mean(x - x.wahr, na.rm = na.rm, ...))
}
