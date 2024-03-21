# Simulation study 1

library(fda)
library(invgamma)
library(dqrng)
library(mvnfast)

source("VB_bs_main.R")
source("run_sims.R")
source("expected_values.R")
source("elbo_formulas.R")
source("elbo.R")


# All data was simulated (including parameters)

# tau2 = 100, sigma2 = 0.01
set.seed(1269)
betas <- round(rnorm(10,0,sqrt(100*(0.5)^2)), 1)

set.seed(1269)
z_i <- rbinom(10, 1, 0.7)

coef <- betas*z_i

# Running VB with sigma2 = 0.01 (p all 1)


#z_initial <-  rbinom(10, 1, 0.7)
p_initial <-  rep(1, 10)


initial_values  <- list(p = rep(p_initial, 5),  delta2 = 90, lambda2 = 100000)

tp1 <- run_sims(S = 50, seed = 3000, save_folder = "Simulation2_VB", scenario = 1, m = 5,  times = 100, ordem = 4, K = 10, coef = coef, sigma = 0.5, mu_ki = 0.7, lambda_1 = 1000, lambda_2 = 999*100, delta_1 = 110, delta_2 = 109*0.25, maxIter = 100, initial_values = initial_values, convergence_threshold = 0.0001, w = NULL)

write(tp1, file = paste0("Simulation2_VB/", "case", 1, "/tempo", ".txt"))



# lambda_1 = 1000, lambda_2 = 999*100, delta_1 = 110, delta_2 = 109*0.01,