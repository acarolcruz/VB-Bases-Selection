# Simulation study 1 - VB compared to Gibbs

library(fda)
library(invgamma)
library(dqrng)
library(mvnfast)

source("VB_bs_main.R")
source("run_sims.R")
source("expected_values.R")
source("elbo_formulas.R")
source("elbo.R")

# Running VB with sigma2 = 0.01 and independence
coef <- c(-2, 0, 1.5, 1.5, 0,-1,-0.5,-1, 0, 0)
z_i <- c(1, 0, 1, 1, 0, 1, 1, 1, 0, 0)

p_initial <-  rep(1, 10)

initial_values <- list(p = rep(p_initial, 5),  delta2 = 5, lambda2 = 10000) #sigma = 0.1

tp1 <- run_sims(S = 100, seed = 3000, save_folder = "Simulation1_VB_gibbs", scenario = 1, m = 5, times = 100, ordem = 4, K = 10, coef = coef, sigma = 0.1, mu_ki = 0.1, lambda_1 = 1e-6, lambda_2 = 1e-6, delta_1 = 100, delta_2 = 99*0.01, maxIter = 500, initial_values = initial_values, convergence_threshold = 0.001, w = NULL, corr = "FALSE", basis_type = "B-splines")

write(tp1, file = paste0("Simulation1_VB_gibbs/", "case", 1, "/tempo", ".txt"))

# Simulation Study 2 - VB compared to Gibbs ----

# Running VB with sigma2 = 0.25 and independence
coef <- c(-2, 0, 1.5, 1.5, 0,-1,-0.5,-1, 0, 0)
z_i <- c(1, 0, 1, 1, 0, 1, 1, 1, 0, 0)

p_initial <-  rep(1, 10)

initial_values2 <- list(p = rep(p_initial, 5),  delta2 = 93, lambda2 = 10000)

tp1 <- run_sims(S = 100, seed = 2024, save_folder = "Simulation1_VB_gibbs", scenario = 2, m = 5, times = 100, ordem = 4, K = 10, coef = coef, sigma = 0.5, mu_ki = 0.1, lambda_1 = 1e-6, lambda_2 = 1e-6, delta_1 = 100, delta_2 = 99*0.25, maxIter = 500, initial_values = initial_values2, convergence_threshold = 0.001, w = NULL, corr = "FALSE", basis_type = "B-splines")

write(tp1, file = paste0("Simulation1_VB_gibbs/", "case", 2, "/tempo", ".txt"))
