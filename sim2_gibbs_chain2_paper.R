# Simulation study 1 - Gibbs

library(fda)
library(invgamma)
library(dqrng)
library(mvnfast)

source("gibbs_basis.R")
source("run_sims_gibbs.R")

# chain 1

# All data was simulated (including parameters)
m = 5
K = 10

coef <- c(-2, 0, 1.5, 1.5, 0,-1,-0.5,-1, 0, 0)
z_i <- c(1, 0, 1, 1, 0, 1, 1, 1, 0, 0)

set.seed(1268)
z_initial <- rbinom(10, 1, 0.7)

# chain 2

initial_values2 <- list(tau2 = 5, sigma2 = 5, theta = rep(4/5,m*K), beta = rep(1,m*K), z = rep(z_initial, m))

tp1 <- run_sims_gibbs(S = 100, seed = 2024, save_folder = "Simulation2_Gibbs_paper", scenario = 2, m = 5, times = 100, ordem = 4, K = 10, coef = coef, sigma = 0.5, mu_ki = 0.1, lambda1 = 0, lambda2 = 0, delta1 = 0, delta2 = 0, maxIter = 10000, initial_values = initial_values2)

write(tp1, file = paste0("Simulation2_Gibbs_paper/", "chain", 2, "/tempo", ".txt"))
