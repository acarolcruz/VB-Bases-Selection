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
betas <- round(rnorm(10,0,sqrt(200*(0.1)^2)), 1)

set.seed(1269)
z_i <- rbinom(10, 1, 0.7)

coef <- betas*z_i

# Running VB with sigma2 = 0.01 (10 random initializations)

n_case = 0
N_case = c()
for(case in 1:10){
  z_vector = rep(0, 10)
  while(sum(z_vector != z_i) > 3){
    n_case = n_case + 1
    set.seed(1268 + n_case)
    #z_initial <-  rbinom(10, 1, 0.7)
    p_initial <-  rbeta(10, 0.7, 0.3)
    
    z_vector = ifelse(p_initial > 0.5, 1, 0)
    #print(n_case)
  }  
  N_case = c(N_case, n_case)
  
  #initial_values  <- list(p = rep(p_initial, 5),  delta2 = 3, lambda2 = 200000)
  
  # tp1 <- run_sims(S = 50, seed = 3000, save_folder = "Simulation1p_info2", scenario = case, m = 5,  times = 100, ordem = 4, K = 10, coef = coef, sigma = 0.1, mu_ki = 0.7, lambda_1 = 1000, lambda_2 = 999*200, delta_1 = 110, delta_2 = 109*0.01, maxIter = 100, initial_values = initial_values, convergence_threshold = 0.0001)
  # 
  # write(tp1, file = paste0("Simulation1p_info2/", "case", case, "/tempo", ".txt"))
  
  initial_values  <- list(p = rep(p_initial, 5),  delta2 = 1, lambda2 = 2400)
  
  tp1 <- run_sims(S = 50, seed = 3000, save_folder = "Simulation1p_noinfo2", scenario = case, m = 5,  times = 100, ordem = 4, K = 10, coef = coef, sigma = 0.1, mu_ki = 0.7, lambda_1 = 10^(-10), lambda_2 = 10^(-10), delta_1 = 10^(-10), delta_2 = 10^(-10), maxIter = 100, initial_values = initial_values, convergence_threshold = 0.0001)
  
  write(tp1, file = paste0("Simulation1p_noinfo2/", "case", case, "/tempo", ".txt"))
  
}

# lambda_1 = 1000, lambda_2 = 999*100, delta_1 = 110, delta_2 = 109*0.01,