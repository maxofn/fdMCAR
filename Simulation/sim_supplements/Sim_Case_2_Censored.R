################################################################################
# This code companions the paper "Testing the Missing Completely at Random
# Assumption for Functional Data" by M. Ofner, S. Hoermann, D. Kraus, and
# D. Liebl.
################################################################################

library(foreach)
library(doParallel)
library(latex2exp)

source("../fdMCAR.R")
source("GenSynth.R")

# Setup parallel computation
num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
cl <- makeCluster(num_cores)
registerDoParallel(cl)

cat("Running Sim Case 2: Censored\n")

set.seed(1)

# Simulation study: Case 2 Censoring --------------------------------------

nr.sim <- 1000
# This takes some time! Decrease nr.sim to reduce computation time.

X <- as.matrix(unname(read.csv("Data/heart_rate.csv", header = FALSE)))

alpha <- 0.05
prob.grid <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1)

# Pre-allocate vector for results
p_mu.100 <- numeric(length(prob.grid))
p_F.100 <- numeric(length(prob.grid))

# Pre-allocate data frames for results
names(p_mu.100) <- prob.grid
names(p_F.100) <- prob.grid

for(j in seq_along(prob.grid)) {
  
  cat(sprintf("Processing j = %d/%d...\n", j, length(prob.grid)))
  
  sim_results <- foreach(i = 1:nr.sim, .combine = rbind) %dopar% {
    # Generate Data
    X.obs <- GenSynth(n = 100, X, prob = prob.grid[j], type.miss = 'censored')$X.obs
    X.obs <- X.obs[which(rowSums(!is.na(X.obs))>0),]
    
    # Perform test
    list  <- fdMCAR.test(X.obs, nSim = 10000, method = "asymptotic")
    
    p_mu_sim <- list$df.test["T_mu", "p.value"] < alpha
    p_F_sim <- list$df.test["T_F", "p.value"] < alpha
    
    # Return a named vector for this iteration
    c(p_mu_sim = p_mu_sim,
      p_F_sim = p_F_sim)
  }
  
  # Aggregate results
  p_mu.100[j] <- mean(sim_results[,"p_mu_sim"], na.rm = TRUE)
  p_F.100[j] <- mean(sim_results[,"p_F_sim"], na.rm = TRUE)
}

# Pre-allocate vector for results
p_mu.250 <- numeric(length(prob.grid))
p_F.250 <- numeric(length(prob.grid))

# Pre-allocate data frames for results
names(p_mu.100) <- prob.grid
names(p_F.100) <- prob.grid

for(j in seq_along(prob.grid)) {
  
  cat(sprintf("Processing j = %d/%d...\n", j, length(prob.grid)))
  
  sim_results <- foreach(i = 1:nr.sim, .combine = rbind) %dopar% {
    # Generate Data
    X.obs <- GenSynth(n = 250, X, prob = prob.grid[j], type.miss = 'censored')$X.obs
    X.obs <- X.obs[which(rowSums(!is.na(X.obs))>0),]
    
    # Perform test
    list  <- fdMCAR.test(X.obs, nSim = 10000, method = "asymptotic")
    
    p_mu_sim <- list$df.test["T_mu", "p.value"] < alpha
    p_F_sim <- list$df.test["T_F", "p.value"] < alpha
    
    # Return a named vector for this iteration
    c(p_mu_sim = p_mu_sim,
      p_F_sim = p_F_sim)
  }
  
  # Aggregate results
  p_mu.250[j] <- mean(sim_results[,"p_mu_sim"], na.rm = TRUE)
  p_F.250[j] <- mean(sim_results[,"p_F_sim"], na.rm = TRUE)
}

stopCluster(cl)

pdf("./Figures/Censored_rej_prob_test.pdf", height = 5, width = 12.27)
par(mfrow = c(1,2), ps = 14)
plot(prob.grid, p_mu.100, ylim = range(0, 1), xlab = "p",
     ylab = "rejection probability", cex = 2,
     main = TeX(r'($n = 100$)'))
points(prob.grid, p_F.100, pch = 3, cex = 2)
legend("bottomright", legend=c(TeX(r'($T_{\mu}$)'),
                               TeX(r'($T_F$)')), pch = c(1,3), horiz = TRUE)
plot(prob.grid, p_mu.250, ylim = range(0, 1), xlab = "p",
     ylab = "rejection probability", cex = 2,
     main = TeX(r'($n = 250$)'))
points(prob.grid, p_F.250, pch = 3, cex = 2)
legend("bottomright", legend=c(TeX(r'($T_{\mu}$)'),
                               TeX(r'($T_F$)')), pch = c(1,3), horiz = TRUE)
dev.off()