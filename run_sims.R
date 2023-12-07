# run a simulated scenario with multiples datasets

run_sims <- function(S = S, save_folder = save_folder, scenario = scenario, m = m, times = times, ordem = ordem, K = K, y = y, B = B, 
                     coef = coef, seed = seed, 
                     sigma = sigma, mu_ki = mu_ki, lambda_1 = lambda_1, lambda_2 = lambda_2, delta_1 = delta_1, delta_2 = delta_2, 
                     maxIter = maxIter, initial_values = initial_values){

  source("sim_data.R")
  source("run_VB.R")
  
  dir.create(save_folder, showWarnings = FALSE)
  dir.create(paste0(save_folder, "/case", scenario), showWarnings = FALSE)

  
  out <- sapply(1:S, run_VB, save_folder = save_folder, scenario = scenario, m = m, times = times, ordem = ordem, K = K, y = y, B = B, coef = coef, seed = seed, 
                sigma = sigma, mu_ki = mu_ki, lambda_1 = lambda_1, lambda_2 = lambda_2, delta_1 = delta_1, delta_2 = delta_2, 
                maxIter = maxIter, initial_values = initial_values)
  
  return(out)
}  
