rmixdistr <- function(n, tail, ei, scale = NULL) {
  if (is.null(scale)) {
    scale <- ei ^ (-1 / tail)
  }
  if (tail == 1) {
    r <- ifelse(runif(n) > ei, 0, rexp(n, rate = 1 / scale))
  } else{
    r <- ifelse(runif(n) > ei, 0, MittagLeffleR::rml(n , tail = tail ,
                                                     scale = scale )
    )
  }
  return(r)
}

pmixdistr <-  function(q, tail, ei, scale = NULL) {
  if (is.null(scale)) {
    scale <- ei ^ (-1 / tail)
  }
  if (tail == 1) {
    p <- 1 - ei + ei * pexp(q, rate = 1 / scale)
    p[q < 0] <- 0
  } else {
    p <- numeric(0)
    if (length(q[q >= 0]) > 0) {
      p[q >= 0] <- 1 - ei + ei * MittagLeffleR::pml(q[q >= 0], tail = tail,
                                                    scale = scale)
    }
    p[q < 0] <- 0
  }
  return(p)
}

qmixdistr <-  function(p, tail, ei, scale = NULL) {
  if (is.null(scale)) {
    scale <- ei ^ (-1 / tail)
  }
  if (tail == 1) {
    p[p < (1 - ei)] <- 1 - ei
    q <- qexp( (p - (1 - ei)) / ei, rate = 1/scale)
  } else {
    q <- numeric(length(p))
    if ( any(p > (1 - ei)) ) {
      q[p > (1 - ei) & p < 1] <- MittagLeffleR::qml((p[p > (1 - ei) & p < 1] -
                                                       (1 - ei)) / ei ,
                                                    tail = tail, scale = scale)
    }
    q[p == 1] <- Inf
  }
  return(q)
}

dmixdistr <- function(x, tail, ei, scale = NULL) {
  if (is.null(scale)) {
    scale <- ei^{-1/tail}
  }
  if(tail == 1){
    d <- ei * dexp(x , rate = 1/scale)
    d[x==0] <- 1-ei
  } else {
    d <- numeric(length(x))
    if (any(x > 0)) {
      d[x > 0] <- ei * MittagLeffleR::dml(x[x > 0], tail = tail, scale = scale)
    }
    d[x == 0] <- 1-ei
  }
  return(d)
}
