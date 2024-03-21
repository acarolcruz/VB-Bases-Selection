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

# tau2 = 100, sigma2 = 0.01
set.seed(1269)
betas <- round(rnorm(10,0,sqrt(100*(0.5)^2)), 1)

set.seed(1269)
z_i <- rbinom(10, 1, 0.7)

coef <- betas*z_i

set.seed(1268)
z_initial <- rbinom(10, 1, 0.7)

# chain 2

initial_values2 <- list(tau2 = 5, sigma2 = 5, theta = rep(4/5,m*K), beta = rep(1,m*K), z = rep(z_initial, m))

tp1 <- run_sims_gibbs(S = 50, seed = 3000, save_folder = "Simulation2_Gibbs", scenario = 2, m = 5, 
                      times = 100, ordem = 4, K = 10, coef = coef, sigma = 0.5, mu_ki = 0.7, 
                      lambda1 = 0, lambda2 = 0, delta1 = 0, delta2 = 0, 
                      maxIter = 10000, initial_values = initial_values2)

write(tp1, file = paste0("Simulation2_Gibbs/", "chain", 2, "/tempo", ".txt"))
