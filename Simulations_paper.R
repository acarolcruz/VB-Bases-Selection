library(fda)
library(invgamma)
library(dqrng)
library(mvnfast)

source("run_sims.R")

#Pedro's coef
coef <- c(-2, 0, 1.5, 1.5, 0,-1,-0.5,-1, 0, 0)
z_i <- c(1, 0, 1, 1, 0, 1, 1, 1, 0, 0)

p_initial <- rep(1, 10)

# sigma small
initial_values <- list(p = rep(p_initial, 5),  delta2 = 5, lambda2 = 10000, w = 10) #sigma = 0.1

tp1 <- run_sims(S = 100, seed = 2024, save_folder = "Simulation1_corr_w", scenario = 1, m = 5, times = 100, ordem = 4, K = 10, coef = coef, sigma = 0.1, mu_ki = 0.5, lambda_1 = 1e-6, lambda_2 = 1e-6, delta_1 = 100, delta_2 = 99*0.01, maxIter = 500, initial_values = initial_values, convergence_threshold = 0.001, w = 2, corr = TRUE, basis_type = "B-splines")

write(tp1, file = paste0("Simulation1_corr_w/", "case", 1, "/tempo", ".txt"))

# misspecified
initial_values <- list(p = rep(p_initial, 5),  delta2 = 5, lambda2 = 10000) #sigma = 0.1

tp1 <- run_sims(S = 100, seed = 2024, save_folder = "Simulation1_corr_w", scenario = 3, m = 5, times = 100, ordem = 4, K = 10, coef = coef, sigma = 0.1, mu_ki = 0.5, lambda_1 = 1e-6, lambda_2 = 1e-6, delta_1 = 100, delta_2 = 99*0.01, maxIter = 500, initial_values = initial_values, convergence_threshold = 0.001, w = 6, corr = "miss", basis_type = "B-splines")

write(tp1, file = paste0("Simulation1_corr_w/", "case", 3, "/tempo", ".txt"))


# sigma large ----
initial_values <- list(p = rep(p_initial, 5),  delta2 = 14, lambda2 = 10000, w = 20) #sigma = 0.2

tp1 <- run_sims(S = 100, seed = 2024, save_folder = "Simulation1_corr_w", scenario = 2, m = 5, times = 100, ordem = 4, K = 10, coef = coef, sigma = 0.2, mu_ki = 0.5, lambda_1 = 1e-6, lambda_2 = 1e-6, delta_1 = 100, delta_2 = 99*(0.2^2), maxIter = 100, initial_values = initial_values, convergence_threshold = 0.01, w = 6, corr = TRUE, basis_type = "B-splines")

write(tp1, file = paste0("Simulation1_corr_w/", "case", 2, "/tempo", ".txt"))

# misspecified
initial_values <- list(p = rep(p_initial, 5),  delta2 = 14, lambda2 = 10000, w = 10) #sigma = 0.2

tp1 <- run_sims(S = 100, seed = 2024, save_folder = "Simulation1_corr_w", scenario = 4, m = 5, times = 100, ordem = 4, K = 10, coef = coef, sigma = 0.2, mu_ki = 0.5, lambda_1 = 1e-6, lambda_2 = 1e-6, delta_1 = 100, delta_2 = 99*(0.2^2), maxIter = 100, initial_values = initial_values, convergence_threshold = 0.01, w = 6, corr = "miss", basis_type = "B-splines")

write(tp1, file = paste0("Simulation1_corr_w/", "case", 4, "/tempo", ".txt"))


# Fourier

coef # no need to specify the coef vector for the only fourier scenario 
z_i <- c(0,1,1,0,0,0,0,0,0,0) #no intercept
p_initial <- rep(1, 10)

initial_values <- list(p = rep(p_initial, 5),  delta2 = 14, lambda2 = 10000, w = 6) #sigma = 0.1

tp1 <- run_sims(S = 100, seed = 2024, save_folder = "Simulation3_corr_w", scenario = 1, m = 5, times = 100, ordem = 4, K = 10, coef = NULL, sigma = 0.1, mu_ki = 0.5, lambda_1 = 1e-6, lambda_2 = 1e-6, delta_1 = 100, delta_2 = 99*(0.1^2), maxIter = 100, initial_values = initial_values, convergence_threshold = 0.001, w = 6, corr = TRUE, basis_type = "Fourier")

write(tp1, file = paste0("Simulation3_corr_w/", "case", 1, "/tempo", ".txt"))

# Inverse (Fourier with Bsplines to fit the data)
coef = NULL # no need to specify the coef vector for the fourier scenario 
z_i <- c(0,1,1,0,0,0,0,0,0,0) # true z vector with no intercept
p_initial <- rep(1, 15)

initial_values <- list(p = rep(p_initial, 5),  delta2 = 3.9, lambda2 = 10000, w = 10) #sigma = 0.1

tp1 <- run_sims(S = 100, seed = 2024, save_folder = "Simulation4_corr_w", scenario = 1, m = 5, times = 100, ordem = 4, K = 15, coef = NULL, sigma = 0.1, mu_ki = 0.5, lambda_1 = 1e-6, lambda_2 = 1e-6, delta_1 = 100, delta_2 = 99*(0.1^2), maxIter = 100, initial_values = initial_values, convergence_threshold = 0.001, w = 6, corr = TRUE, basis_type = "inverse")

write(tp1, file = paste0("Simulation4_corr_w/", "case", 1, "/tempo", ".txt"))


