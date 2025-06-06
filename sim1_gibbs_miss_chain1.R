# Simulation study 1 - Gibbs - misspecified

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

z_initial <- rep(1, 10)

initial_values <- list(tau2 = 1, sigma2 = 1, theta = rep(1/5, m*K), beta = rep(-1, m*K), z = rep(z_initial, m))

tp <- run_sims_gibbs(S = 100, seed = 2024, save_folder = "Simulation1_Gibbs_miss_paper",
                     scenario = 1, m = 5, 
                     times = 100, ordem = 4, K = 10, coef = coef, sigma = 0.1, mu_ki = 0.5, 
                     lambda1 = 0, lambda2 = 0, delta1 = 0, delta2 = 0, 
                     maxIter = 10000, initial_values = initial_values, w = 6, 
                     corr = "miss", basis_type = "B-splines")

write(tp, file = paste0("Simulation1_Gibbs_miss_paper/", "chain", 1, "/tempo", ".txt"))
