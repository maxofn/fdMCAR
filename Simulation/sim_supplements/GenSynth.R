GenSynth <- function(n = 100, X, prob = 1, type.miss = 'mcar', a = -5, b = 15) {
  # Generate observations
  
  p <- dim(X)[2]
  grid <- 1:p
  sampleidx <- sample(dim(X)[1], size = n, replace = TRUE)
  mu <- colMeans(X, na.rm = TRUE)
  X <- X[sampleidx,]
  
  # Miss observations
  X.obs <- X
  
  if(type.miss != 'mcar' & type.miss != 'censored') {
    stop("Unknown missingness mechansim")
  }
  
  if(type.miss == 'mcar') {
    for(i in 1:n) {
      u2 <- sample(grid[floor(p/4):p], size = 1)
      obs <- which(grid <= u2)
      X.obs[i, -obs] <- NA
    }
  }
  
  if(type.miss == 'censored') {
    X_sweep <- sweep(X, 2, mu)
    for(i in 1:n) {
      if(sample(c(TRUE, FALSE), 1, prob = c(prob, 1 - prob))) {
        if(X_sweep[i,1] < 10) {
          u2 <- sample(grid[floor(p/4):p], size = 1)
          obs <- which(grid <= u2)
          X.obs[i, -obs] <- NA
        }
        
        #miss <- which(X_sweep[i,] <= a | X_sweep[i,] >= b)
        #X.obs[i,][miss] <- NA
      }
    }
  }
  
  return(list(X = X, X.obs = X.obs))
}

