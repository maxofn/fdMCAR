################################################################################
# This code companions the paper "Testing the Missing Completely at Random
# Assumption for Functional Data" by M. Ofner, S. Hoermann, D. Kraus, and
# D. Liebl.
################################################################################

library(foreach)
library(doParallel)
library(latex2exp)

source("../fdMCAR.R")
source("GenBM.R")


# Setup parallel computation
num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
cl <- makeCluster(num_cores)
registerDoParallel(cl)

cat("Running Sim Case 2: Censored\n")

set.seed(1)

# Simulation case 2 Censored: Rejection probability ----------------------

nr.sim <- 1000
# This takes some time! Decrease nr.sim to reduce computation time.

# a = -1
# b in [1, 2]

alpha <- 0.05
a <- -1
b.grid <- seq(1, 2, length.out = 11)
n.grid <- length(b.grid)

# Pre-allocate vector for results
p_mu.100 <- numeric(length(b.grid))
p_F.100 <- numeric(length(b.grid))

# Pre-allocate data frames for results
names(p_mu.100) <- b.grid
names(p_F.100) <- b.grid

for(b in seq_along(b.grid)) {
  
  cat(sprintf("Processing b = %d/%d...\n", b, length(b.grid)))
  
  sim_results <- foreach(i = 1:nr.sim, .combine = rbind) %dopar% {
    # Generate Data
    X.obs <- GenBM(n = 100, type.miss = 'censored', a = a, b = b.grid[b])$X.obs
    
    # Perform Test
    list  <- fdMCAR.test(X.obs, nSim = 10000, method = "asymptotic")
    
    p_mu_sim <- list$df.test["T_mu", "p.value"] < alpha
    p_F_sim <- list$df.test["T_F", "p.value"] < alpha
    
    # Return a named vector for this iteration
    c(p_mu_sim = p_mu_sim,
      p_F_sim = p_F_sim)
  }
  
  # Aggregate results
  p_mu.100[b] <- mean(sim_results[,"p_mu_sim"], na.rm = TRUE)
  p_F.100[b] <- mean(sim_results[,"p_F_sim"], na.rm = TRUE)
}

# Pre-allocate vector for results
p_mu.500 <- numeric(length(b.grid))
p_F.500 <- numeric(length(b.grid))

# Pre-allocate data frames for results
names(p_mu.500) <- b.grid
names(p_F.500) <- b.grid

for(b in seq_along(b.grid)) {
  
  cat(sprintf("Processing b = %d/%d...\n", b, length(b.grid)))
  
  sim_results <- foreach(i = 1:nr.sim, .combine = rbind) %dopar% {
    # Generate Data
    X.obs <- GenBM(n = 500, type.miss = 'censored', a = a, b = b.grid[b])$X.obs
    
    # Perform Test
    list  <- fdMCAR.test(X.obs, nSim = 10000, method = "asymptotic")
    
    p_mu_sim <- list$df.test["T_mu", "p.value"] < alpha
    p_F_sim <- list$df.test["T_F", "p.value"] < alpha
    
    # Return a named vector for this iteration
    c(p_mu_sim = p_mu_sim,
      p_F_sim = p_F_sim)
  }
  
  # Aggregate results
  p_mu.500[b] <- mean(sim_results[,"p_mu_sim"], na.rm = TRUE)
  p_F.500[b] <- mean(sim_results[,"p_F_sim"], na.rm = TRUE)
}

# Simulation Case 2 censored: Sensitivity clustering ---------------------

set.seed(1)

alpha <- 0.05
a <- -1
b <- 1
d.grid <- seq(0.5, 1, by = 0.1)

p <- 100

p_Fclust.100  <- numeric(length(d.grid))
p_Fdeter.100  <- numeric(length(d.grid))

# Pre-allocate data frames for results
names(p_Fclust.100) <- d.grid
names(p_Fclust.100) <- d.grid

for(d in seq_along(d.grid)) {
  
  cat(sprintf("Processing d = %d/%d...\n", d, length(d.grid)))
  
  sim_results <- foreach(i = 1:nr.sim, .combine = rbind) %dopar% {
    # Generate Data
    X.obs <- GenBM(n = 100, p = p, type.miss = 'censored', a = a, b = b)$X.obs
    
    list  <- fdMCAR.test(X.obs, nSim = 10000, method = "asymptotic")
    p_Fclust <- list$df.test["T_F", "p.value"] < alpha
    
    obs.mat <- !is.na(X.obs)
    A <- which(rowSums(obs.mat) < d.grid[d] * p)
    list  <- fdMCAR.test(X.obs, A = A, nSim = 10000, method = "asymptotic")
    p_Fdeter  <- list$df.test["T_F", "p.value"] < alpha
    
    # Return a named vector for this iteration
    c(p_Fclust = p_Fclust,
      p_Fdeter = p_Fdeter)
  }
  
  # Aggregate results
  p_Fdeter.100[d] <- mean(sim_results[,"p_Fdeter"], na.rm = TRUE)
  p_Fclust.100[d] <- mean(sim_results[,"p_Fclust"], na.rm = TRUE)
}

p_Fclust.250  <- numeric(length(d.grid))
p_Fdeter.250  <- numeric(length(d.grid))

# Pre-allocate data frames for results
names(p_Fclust.250) <- d.grid
names(p_Fclust.250) <- d.grid

for(d in seq_along(d.grid)) {
  
  cat(sprintf("Processing d = %d/%d...\n", d, length(d.grid)))
  
  sim_results <- foreach(i = 1:nr.sim, .combine = rbind) %dopar% {
    # Generate Data
    X.obs <- GenBM(n = 250, p = p, type.miss = 'censored', a = a, b = b)$X.obs
    
    list  <- fdMCAR.test(X.obs, nSim = 10000, method = "asymptotic")
    p_Fclust <- list$df.test["T_F", "p.value"] < alpha
    
    obs.mat <- !is.na(X.obs)
    A <- which(rowSums(obs.mat) < d.grid[d] * p)
    list  <- fdMCAR.test(X.obs, A = A, nSim = 10000, method = "asymptotic")
    p_Fdeter  <- list$df.test["T_F", "p.value"] < alpha
    
    # Return a named vector for this iteration
    c(p_Fclust = p_Fclust,
      p_Fdeter = p_Fdeter)
  }
  
  # Aggregate results
  p_Fdeter.250[d] <- mean(sim_results[,"p_Fdeter"], na.rm = TRUE)
  p_Fclust.250[d] <- mean(sim_results[,"p_Fclust"], na.rm = TRUE)
}

stopCluster(cl)

pdf(paste0("Figures/Censored_rej_prob.pdf"), height = 8, width = 8)
par(mfrow = c(2,2), ps = 14, mar = c(4, 4.1, 4.1, 2.1))
plot(b.grid, p_mu.100, ylim = range(0, 1), xlab = "b",
     ylab = "rejection probability", cex = 2,
     main = TeX(r'($n = 100$)'))
points(b.grid, p_F.100, pch = 3, cex = 2)
legend("bottomright", legend=c(TeX(r'($T_F$)'),
                               TeX(r'($T_{\mu}$)')), pch = c(3,1),
       horiz = TRUE)
plot(b.grid, p_mu.500, ylim = range(0, 1), xlab = "b",
     ylab = "rejection probability", cex = 2,
     main = TeX(r'($n = 500$)'))
points(b.grid, p_F.500, pch = 3, cex = 2)
legend("bottomright", legend=c(TeX(r'($T_F$)'),
                               TeX(r'($T_{\mu}$)')), pch = c(3,1),
       horiz = TRUE)
plot(d.grid, p_Fdeter.100, ylim = range(0, 1), xlab = TeX(r'($\delta$)'),
     ylab = "rejection probability", cex = 2,
     main = TeX(r'($n = 100$)'), pch = 3)
abline(h = mean(p_Fclust.100), lty = 1, lwd = 3)
legend("bottomright", legend = c("det","clust"), pch = c(3,NA), lty = c(NA,1),
       lwd = c(NA, 2), horiz = TRUE)
plot(d.grid, p_Fdeter.250, ylim = range(0, 1), xlab = TeX(r'($\delta$)'),
     ylab = "rejection probability", cex = 2,
     main = TeX(r'($n = 250$)'), pch = 3)
abline(h = mean(p_Fclust.250), lty = 1, lwd = 3)
legend("bottomright", legend = c("det","clust"), pch = c(3,NA), lty = c(NA,1),
       lwd = c(NA, 2), horiz = TRUE)
dev.off()

set.seed(1)
a <- -1
b <- 2

X.obs <- GenBM(n = 100, type.miss = 'censored', a = a, b = b)$X.obs
list <- fdMCAR.test(X.obs, A = which(rowSums(is.na(X.obs))==0), nSim = 1000, return.SCB = TRUE)

n <- 100
diff <- list$mean.diff$diff
q    <- list$mean.diff$q

grid <- seq(0, 1, length.out = 100)

pdf(paste0("Figures/Censored_example.pdf"), height = 5, width = 12.27)
par(mfrow = c(1,2), ps = 14)

sample_A <- which(rowSums(is.na(X.obs)) == 0)
sample_B <- which(rowSums(is.na(X.obs)) != 0)

distance <- (max(X.obs, na.rm = TRUE) - min(X.obs, na.rm = TRUE))/2
rg <- c(-distance, distance)

matplot(grid, t(X.obs[sample_A,]), lty = 1, type = "l",
        xlab = "time", ylab = "value", col = "gray", lwd = 2,
        yaxt = "n")
axis(2, at = c(-1, 0, 1, 2), labels = c(-1,0,1, 2))
matlines(grid, t(X.obs[sample_B,]), lty = 1, type = "l", col = "black",
         lwd = 2)
plot(grid, diff, type = "l", xlab = "time", ylab = "difference in means",
     ylim = rg, lwd = 4, yaxt = "n")
axis(2, at = c(-1, 0, 1), labels = c(-1,0,1))
lines(grid, diff +  q/sqrt(n), lty = 2, lwd = 2)
lines(grid, diff -  q/sqrt(n), lty = 2, lwd = 2)
abline(h = 0, lty = 1, lwd = 1)
dev.off()