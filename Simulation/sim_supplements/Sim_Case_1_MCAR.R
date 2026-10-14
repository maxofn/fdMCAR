################################################################################
# This code companions the paper "Testing the Missing Completely at Random
# Assumption for Functional Data" by M.Ofner, S. Hoermann, D. Kraus, and
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

cat("Running Sim Case 1: MCAR\n")

set.seed(1)

# Simulation study: Case 1 MCAR -------------------------------------------

X <- as.matrix(unname(read.csv("Data/heart_rate.csv", header = FALSE)))

nr.sim <- 1000
# This takes some time! Decrease nr.sim to reduce computation time.

n.vec <- c(100, 250, 500)
p <- 100
alpha <- 0.05

sim.result_asymp <- data.frame("T_mu"  = rep(NA, length(n.vec)),
                               "T_F"  = rep(NA, length(n.vec)))

sim.result_boot  <- data.frame("T_mu"  = rep(NA, length(n.vec)),
                               "T_F"  = rep(NA, length(n.vec)))

# Pre-allocate data frames for results
rownames(sim.result_asymp) <- n.vec
rownames(sim.result_boot)  <- n.vec

for(n.idx in seq_along(n.vec)) {
  n <- n.vec[n.idx]
  
  cat(sprintf("Processing n = %d [%d/%d]...\n", n, n.idx, length(n.vec)))
  
  sim_results <- foreach(i = 1:nr.sim, .combine = rbind) %dopar% {
    # Generate Data
    X.obs <- GenSynth(n, X, type.miss = 'mcar')$X.obs
    
    # Run Asymptotic Test
    res_asymp <- fdMCAR.test(X.obs, nSim = 10000, method = "asymptotic")
    p_mu_asymp <- res_asymp$df.test["T_mu", "p.value"]
    p_F_asymp <- res_asymp$df.test["T_F", "p.value"]
    
    # Run Bootstrap Test
    res_boot <- fdMCAR.test(X.obs, nBoot = 10000, method = "bootstrap")
    p_mu_boot <- res_boot$df.test["T_mu", "p.value"]
    p_F_boot <- res_boot$df.test["T_F", "p.value"]
    
    # Return a named vector for this iteration
    c(p_mu_asymp = p_mu_asymp, 
      p_F_asymp = p_F_asymp, 
      p_mu_boot  = p_mu_boot, 
      p_F_boot  = p_F_boot)
  }
  
  # Aggregate results
  sim.result_asymp$T_mu[n.idx] <- mean(sim_results[, "p_mu_asymp"] < alpha, na.rm = TRUE)
  sim.result_asymp$T_F[n.idx] <- mean(sim_results[, "p_F_asymp"] < alpha, na.rm = TRUE)
  
  sim.result_boot$T_mu[n.idx]  <- mean(sim_results[, "p_mu_boot"] < alpha,  na.rm = TRUE)
  sim.result_boot$T_F[n.idx]  <- mean(sim_results[, "p_F_boot"] < alpha,  na.rm = TRUE)
}

stopCluster(cl)

rownames(sim.result_asymp) <- n.vec
write.csv(sim.result_asymp, "Results/MCAR_rej_prob_asymptotic_synth.csv")

rownames(sim.result_boot) <- n.vec
write.csv(sim.result_boot, "Results/MCAR_rej_prob_bootstrap_synth.csv")
