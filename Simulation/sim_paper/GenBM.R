#' Simulate Brownian motions.
#'
#' @param n         Sample size (n = 100 by default).
#' @param p         Number of regularly spaced grid points within [0,1] (p = 100
#'                  by default).
#' @param type.miss Type of missigness mechanism either 'mcar' (default) or 
#'                  'censored' (see paper).
#' @param a         Parameter specifying lower bound of censoring interval.
#' @param b         Parameter specifying upper bound of censoring interval.
#'
#' @return A list with X (true data matrix without NA's) and X.obs (with NA's).
#' @export
#'
#' @examples
GenBM <- function(n = 100, p = 100, type.miss = 'mcar', a = -1, b = 1) {
  # Generate observations
  X <- matrix(0, nrow = n, ncol = p)
  
  grid <- seq(0, 1, length.out = p)
  dt <- 1 / (p - 1)
  
  # Simulate Brownian motion
  for (i in 1:n) {
    dB <- sqrt(dt) * rnorm(p - 1)
    X[i, 2:p] <- cumsum(dB)
  }
  
  # Miss observations
  X.obs <- X
  
  if (type.miss != "mcar" & type.miss != "censored" & type.miss != "interval") {
    stop("Unknown missingness mechansim")
  }
  
  if(type.miss == 'mcar') {
    for(i in 1:n) {
      if(sample(c(TRUE, FALSE), size = 1)) {
        u1 <- sample(grid, size = 1)
        u2 <- sample(grid, size = 1)
        obs <- which(grid >= min(u1,u2) & grid < max(u1,u2))
        X.obs[i, -obs] <- NA
      }
    }
  }
  
  if(type.miss == 'censored') {
    obs <- which(X > a & X < b)
    X.obs[-obs] <- NA
  }
  
  if(type.miss == 'interval') {
    for(i in 1:n) {
      miss <- which(X[i,] <= a | X[i,] >= b)
      if(length(miss) > 0) {
        X.obs[i,min(miss):p] <- NA
      }
    }
  }
  
  return(list(X = X, X.obs = X.obs))
}
