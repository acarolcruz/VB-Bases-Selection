# Run VB algorithm for a simulation and save results in a folder

run_VB <- function(sim = sim, save_folder = save_folder, scenario = scenario, m = m, times = times, ordem = ordem, K = K, y = y, B = B, coef = coef, seed = seed, 
                         sigma = sigma, mu_ki = mu_ki, lambda_1 = lambda_1, lambda_2 = lambda_2, delta_1 = delta_1, delta_2 = delta_2, 
                         maxIter = maxIter, initial_values = initial_values){
  
  tt <- Sys.time()
  
  # Generate the data for one dataset 
  data <- sim_data(m = m, times = times, ordem = ordem, K = K, coef = coef, seed = seed + sim, sigma = sigma)
  
  # Run VB basis selection algorithm
  out <- vb_bs(y = data$y, B = data$B, m = m, mu_ki = mu_ki, lambda_1 = lambda_1, lambda_2 = lambda_2, delta_1 = delta_1, delta_2 = delta_2, maxIter = maxIter, K = K, initial_values = initial_values)
  
  tt <- difftime(Sys.time(), tt, units = 'mins')
  
  save(out, file = paste0(save_folder, "/case", scenario, "/Sim", sim, "VBresults.RData"))
  
  return(tt)
}