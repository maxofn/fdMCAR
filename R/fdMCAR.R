# Functions ---------------------------------------------------------------

#' MCAR tests for functional data
#' 
#' @param X.obs       Incompletely observed data matrix with NA's;
#' @param A           Set containing the indices of curves belonging to the
#'                    first group. If A = NULL (default), it is chosen via
#'                    a clustering procedure;
#' @param method      Method for approximating the distribution of the test
#'                    statistic. Either method = "asymptotic" (default) for
#'                    asymptotic distribution or "bootstrap" for a bootstrap
#'                    procedure;
#' @param nSim        Number of simulated examples for approximating the
#'                    asymptotic distribution (if method = "asymptotic");
#' @param nBoot       Number of bootstrap repetitions (if method = "bootstrap");
#' @param return.SCB  Return difference in means and required quantile for 
#'                    constructing simultaneous confidence bands.
#'
#' @returns A data frame with the value of test statistics and p-values.
#' @export
#'
#' @examples
#' set.seed(1)
#' data <- GenBM(n = 100, p = 100, type.miss = 'mcar')
#' fdMCAR.test(data$X.obs)

fdMCAR.test <- function(X.obs,
                        A = NULL,
                        method = "asymptotic",
                        nSim = 1000,
                        nBoot = 1000, 
                        return.SCB = FALSE) {
  
  # Parameters
  
  nMC <- 500 # number of points for approximating the double integral
  
  n <- dim(X.obs)[1] # sample size
  p <- dim(X.obs)[2] # number of grid points
  
  # sample_A ... Indices of curves in group A
  # sample_B ... Indices of curves in group B
  
  centers <- NULL
  if(is.null(A)) {
    kMed <- kMedians(!is.na(X.obs)) 
    cluster <- kMed$cluster
    centers <- kMed$centers
    sample_A <- which(cluster == 1)
  } else {
    sample_A <- A
  }
  
  sample_B <- (1:n)[-sample_A] # put the rest into the second group
  
  nA <- length(sample_A) # number of curves in group A
  nB <- length(sample_B) # number of curves in group B
  
  # Choose domain for the comparison
  I <- which(colSums(!is.na(X.obs[sample_A,,drop = FALSE])) >= nA * 0.25 &
             colSums(!is.na(X.obs[sample_B,,drop = FALSE])) >= nB * 0.25)

  # Check if A or I are degenerate
  if(length(I) <= 1 | length(sample_A) <= 1 | length(sample_A) >= n-1) {
    warning("A and/or domain I are degenerate. Try a different partition.")
    df.test <- data.frame("T.value" = c(NA, NA, NA),
                          "p.value" = c(NA, NA, NA))
    list <- list(df.test = df.test, nA = nA, nB = nB, mean.diff = NA,
                 lI = length(I))
    return(list)
  }
  
  # Estimate mean functions
  
  mu_A <- colMeans(X.obs[sample_A,], na.rm = TRUE)
  mu_B <- colMeans(X.obs[sample_B,], na.rm = TRUE)
  
  diff = mu_A - mu_B
  diff[-I] <- NA
  
  T_L2 <- n * sum((mu_A[I] - mu_B[I])^2) / p
  T_mu  <- sqrt(n) * max(abs(mu_A[I] - mu_B[I]))
  
  # Estimate distribution functions
  z.sample <- stats::rnorm(nMC,
                           mean = mean(colMeans(X.obs, na.rm = TRUE), na.rm = TRUE), 
                           sd = sqrt(max(apply(X.obs, 2, var, na.rm = TRUE), na.rm = TRUE)))
  
  z.sample  <- sort(z.sample) # Gaussian probability measure
  p.sample  <- sample(I, nMC, replace = TRUE)
  mc.points <- cbind(p.sample, z.sample)
  
  p.hat_A <- 1/n * colSums(!is.na(X.obs[sample_A,]))
  p.hat_B <- 1/n * colSums(!is.na(X.obs[sample_B,]))
  
  Xind <- t(apply(X.obs[, mc.points[,1]], 1, 
                  function(x) {(x <= mc.points[,2])}))
  F.hat_A <- colMeans(Xind[sample_A,], na.rm = TRUE)
  F.hat_B <- colMeans(Xind[sample_B,], na.rm = TRUE)
  
  T_F <- n * sum((F.hat_A - F.hat_B)^2) / nMC
  
  if(method == "bootstrap") {
    
    ############################################################################
    # We approximate the distributions via bootstrap.
    # This implements Algorithms 5 and 6 of the supplementary material.
    ############################################################################
        
    X_sweep <- X.obs
    X_sweep[sample_A,] <- sweep(X.obs[sample_A,], 2, mu_A)
    X_sweep[sample_B,] <- sweep(X.obs[sample_B,], 2, mu_B)
    
    Xind[sample_A,] <- sweep(Xind[sample_A,], 2, F.hat_A)
    Xind[sample_B,] <- sweep(Xind[sample_B,], 2, F.hat_B)
    
    T_L2.b <- numeric(nBoot)
    T_mu.b  <- numeric(nBoot)
    T_F.b  <- numeric(nBoot)
    
    for(b in 1:nBoot) {
      # Draw resamples separately
      idx.b_A <- sample(sample_A, size = nA, replace = TRUE)
      idx.b_B <- sample(sample_B, size = nB, replace = TRUE)
      
      # Estimate mean functions
      
      mu.b_A <- colMeans(X_sweep[idx.b_A,], na.rm = TRUE)
      mu.b_B <- colMeans(X_sweep[idx.b_B,], na.rm = TRUE)
      
      T_L2.b[b] <- n * sum((mu.b_A[I] - mu.b_B[I])^2) / p
      T_mu.b[b]  <- sqrt(n) * max(abs(mu.b_A[I] - mu.b_B[I]))
      
      # Estimate distribution function

      F.hat.b_A <- colMeans(Xind[idx.b_A,], na.rm = TRUE)
      F.hat.b_B <- colMeans(Xind[idx.b_B,], na.rm = TRUE)

      T_F.b[b] <- n * sum((F.hat.b_A - F.hat.b_B)^2) / nMC
    }
    
    p_L2 <- length(which(T_L2.b > T_L2))/nBoot
    p_mu  <- length(which(T_mu.b > T_mu))/nBoot
    p_F  <- length(which(T_F.b > T_F))/nBoot
    
  } else {
    
    ############################################################################
    # We consider the asymptotic distribution.
    # This implements Algorithm 2 and 3 of the supplementary material.
    ############################################################################
    
    # (1) Mean comparison wrt L2 norm
    
    XO <- X.obs[,I]
    XO[sample_A,] <- sweep(XO[sample_A,], 2, mu_A[I])
    XO[sample_B,] <- sweep(XO[sample_B,], 2, mu_B[I])
    XO[is.na(X.obs[,I])] <- 0
    
    covL2 <- (t(XO[sample_A,]) %*% XO[sample_A,] / (n * p.hat_A[I] %*% t(p.hat_A[I])) +
              t(XO[sample_B,]) %*% XO[sample_B,] / (n * p.hat_B[I] %*% t(p.hat_B[I]))) / p
    
    fit.svd <- svd(covL2, nu = 0)
    
    fit.ev   <- fit.svd$d
    n.ev <- max(min(which(cumsum(fit.ev)/sum(fit.ev) > 0.99)), 2) # FEV: explain 99% of var
    
    xi <- matrix(stats::rnorm(nSim * n.ev), nrow = nSim, ncol = n.ev)
    chisq.mat <- xi^2
    T_L2.sim <- chisq.mat %*% fit.ev[1:n.ev]
    p_L2 <- length(which(T_L2.sim > T_L2))/nSim
    
    # (2) Mean comparison wrt sup norm
    
    phi <- fit.svd$v[, 1:n.ev]
    Zs <- xi %*% diag(sqrt(fit.ev[1:n.ev] * p)) %*% t(phi)
    Zsmax <- apply(Zs, 1, function(v){max(abs(v))})

    p_mu <- length(which(Zsmax > T_mu))/nSim
    
    # (3) Comparison of distributions
    
    Xind <- t(apply(X.obs[, mc.points[,1]], 1,
                    function(x) {(x <= mc.points[,2])}))
    Xind[sample_A,] <- sweep(Xind[sample_A,], 2, F.hat_A)
    Xind[sample_B,] <- sweep(Xind[sample_B,], 2, F.hat_B)
    Xind[is.na(Xind)] <- 0
    
    covmu <- (t(Xind[sample_A,]) %*% Xind[sample_A,] / (n * p.hat_A[mc.points[,1]] %*% t(p.hat_A[mc.points[,1]])) +
              t(Xind[sample_B,]) %*% Xind[sample_B,] / (n * p.hat_B[mc.points[,1]] %*% t(p.hat_B[mc.points[,1]]))) / nMC
    fit.svd.mu <- svd(covmu, nu = 0, nv = 0)
    
    fit.ev.mu <- fit.svd.mu$d
    n.ev.mu <- max(min(which(cumsum(fit.ev.mu)/sum(fit.ev.mu) > 0.99)),2) # FEV: explain 99% of var
    
    chisq.mat <- matrix(stats::rchisq(nSim * n.ev.mu, df = 1), nrow = nSim,
                        ncol = n.ev.mu)
    
    T_F.sim <- chisq.mat %*% fit.ev.mu[1:n.ev.mu]
    p_F <- length(which(T_F.sim > T_F))/nSim
  }
  
  # Return objects
  
  df.test <- data.frame("T.value" = c(T_L2, T_mu, T_F),
                        "p.value" = c(p_L2, p_mu, p_F))
  rownames(df.test) <- c("T_L2", "T_mu", "T_F")
  
  mean.diff <- list()
  
  if(return.SCB) {
    
    ############################################################################
    # We compute an asymptotic simultaneous confidence band.
    # This implements Algorithm 4 and 7 of the supplementary material.
    ############################################################################
    
    if(method == "bootstrap") {
      quant <- quantile(T_mu.b, 0.95, na.rm = TRUE)
    }
    else {
      quant <- quantile(Zsmax, 0.95, na.rm = TRUE)
    }
    mean.diff <- list(q = quant, diff = diff)
  }
  
  list <- list(df.test = df.test, nA = nA, nB = nB, mean.diff = mean.diff,
               lI = length(I), sample_A = sample_A, centers = centers)
  return(list)
}

#' Clustering observation sets based on kMedians
#'
#' @param obs.mat Data matrix consisting of TRUE (observed) and FALSE (missing).
#'
#' @returns       A vector containing values 1 (group A) and 2 (group B).
#' @export
#'
#' @examples
#' set.seed(1)
#' data <- (n = 100, p = 100, type.miss = 'mcar')
#' kMedians(!is.na(data$X.obs))
#' 
kMedians <- function(obs.mat) {
  n <- dim(obs.mat)[1]
  p <- dim(obs.mat)[2]
  
  # Use 5 starts
  best_totwss <- Inf
  best_result <- NULL
  
  for(start in 1:5) {
    converged <- FALSE
    
    # Initialize clusters via kmeans++
    k1 <- sample(n, 1)
    centers <- matrix(NA, nrow = 2, ncol = p)
    centers[1,] <- obs.mat[k1,]
    dists <- as.matrix(dist(rbind(centers[1,], obs.mat), method = "manhattan"))[1,-1]
    
    k2 <- sample(1:n, 1, prob = dists^2/(sum(dists^2)))
    centers[2,] <- obs.mat[k2,]
    
    while(converged == FALSE) {
      cluster <- rep(NA, n)
      
      # Assign sets to nearest center
      dists <- as.matrix(dist(rbind(centers, obs.mat), method = "manhattan"))[1:2, (2+1):(2+n)]
      cluster <- max.col(-t(dists))
      
      # Update centers as Vorob'ev medians of each cluster
      new.centers <- do.call(rbind, lapply(1:2, function(j) {
        rows <- obs.mat[cluster == j, , drop = FALSE]
        round(colMeans(rows), 0)
      }))
      
      # Check for convergence
      if(identical(new.centers, centers)) {
        converged <- TRUE
      } else {
        centers <- new.centers
      }
    }
    
    totwss <- sum(sapply(1:2, function(j) {
      rows <- obs.mat[cluster == j, , drop = FALSE]
      if (nrow(rows) == 0) return(0)
      sum(rowSums(abs(rows - matrix(centers[j, ], nrow(rows), p, byrow = TRUE))))
    }))
    
    if (totwss < best_totwss) {
      best_totwss <- totwss
      best_result <- list(cluster = cluster,
                          centers = centers,
                          tot.withinss = totwss)
    }
  }
  
  cluster <- best_result$cluster
  centers <- best_result$centers
  
  # Ensure that largest sets are in group A
  max_size_A <- max(rowSums(obs.mat[which(cluster == 1),,drop = FALSE]))
  max_size_B <- max(rowSums(obs.mat[which(cluster == 2),,drop = FALSE]))
  
  if(max_size_B > max_size_A) {
    cluster <- 3 - cluster
  }
  
  return(list("cluster" = cluster, "centers" = centers))
}

#' Plot barcodes of clusters
#'
#' @param obs.mat Data matrix consisting of TRUE (observed) and FALSE (missing).
#' @param nheight Height of y-axis in plot.
#' @param center  Left and right point of set-valued center (only 
#'                intervals are supported).
#' @param label   Label of set-valued center               

#' @return Plots the observations sets.
#' @export
#'
#' @examples
#' set.seed(1)
#' data <- (n = 100, p = 100, type.miss = 'mcar')
#' PlotSets(!is.na(data$X.obs), dim(data$X.obs)[1]) 
PlotSets <- function(obs.mat, nheight, center = NULL, label = "") {
  n <- dim(obs.mat)[1]
  p <- dim(obs.mat)[2]
  obs_NA <- obs.mat
  obs_NA[obs.mat == FALSE] <- NA
  
  size_ordered <- order(rowSums(obs_NA, na.rm = TRUE))
  
  par(mar = c(5.1, 1, 4.1, 1))
  plot(1:p, obs_NA[size_ordered[1],], ylim = c(0, nheight), type = "l",
   ylab = "", xlab = "", xaxt = "n", yaxt = "n")
  for(idx in 2:n) {
    lines(1:p, obs_NA[size_ordered[idx],] + (idx - 1))
  }
  
  if(!is.null(center)) {
    par(xpd = TRUE) 
    
    y_margin <- -(nheight * 0.1) 
    
    segments(x0 = center[1], y0 = y_margin, x1 = center[2], y1 = y_margin, lwd = 2)
    
    text(x = (center[2] - center[1]) / 2, y = y_margin, labels = label, pos = 1) 
  }
  
  return()
}
