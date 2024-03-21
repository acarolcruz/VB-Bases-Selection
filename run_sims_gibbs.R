# run gibbs a simulated scenario with multiples datasets

run_sims_gibbs <- function(S = S, save_folder = save_folder, scenario = scenario, m = m, times = times, ordem = ordem, K = K, y = y, B = B, coef = coef, seed = seed, sigma = sigma, mu_ki = mu_ki, lambda1 = lambda1, lambda2 = lambda2, delta1 = delta1, delta2 = delta2, maxIter = maxIter, initial_values = initial_values){
  
  source("sim_data.R")
  source("run_gibbs.R")
  
  
  dir.create(save_folder, showWarnings = FALSE)
  dir.create(paste0(save_folder, "/chain", scenario), showWarnings = FALSE)
  
  
  library(parallel)
  cls <- makeForkCluster(detectCores())
  outt <- parSapply(cls, 1:S, run_gibbs, save_folder = save_folder, scenario = scenario, m = m, times = times, ordem = ordem, K = K, 
                    coef = coef, seed = seed, sigma = sigma, mu_ki = mu_ki, lambda1 = lambda1, lambda2 = lambda2, delta1 = delta1, delta2 = delta2,
                    maxIter = maxIter, initial_values = initial_values)
  stopCluster(cls)
  
  return(outt)
}  
