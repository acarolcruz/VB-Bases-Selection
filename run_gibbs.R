run_gibbs <- function(sim = sim, save_folder = save_folder, scenario = scenario, m = m, times = times, ordem = ordem, K = K, coef = coef, seed = seed, sigma = sigma, mu_ki = mu_ki, lambda1 = lambda1, lambda2 = lambda2, delta1 = delta1, delta2 = delta2, maxIter = maxIter,
initial_values = initial_values, w = w, corr = corr, basis_type = basis_type){
  
  tt <- Sys.time()
  
  # Generate the data for one dataset 
  if(corr == "miss"){
    source("sim_data_corr.R")
    source("aux_fun_corr.R")
    data <- sim_data_corr(m = m, times = times, ordem = ordem, K = K, coef = coef, seed = seed + sim, sigma = sigma, w = w, basis_type = basis_type)
  } else{
    source("sim_data.R")
    data <- sim_data(m = m, times = times, ordem = ordem, K = K, coef = coef, seed = seed + sim, sigma = sigma)}
  
  # Run Gibbs basis selection algorithm
  out <- gibbs_basis(y = data$y, B = data$B, m = m, mu_ki = mu_ki, lambda1 = lambda1, lambda2 = lambda2, delta1 = delta1, delta2 = delta2, maxIter = maxIter, K = K, initial_values = initial_values)
  
  tt <- difftime(Sys.time(), tt, units = 'mins')
  
  save(out, file = paste0(save_folder, "/chain", scenario, "/Sim", sim, "Gibbsresults.RData"))
  
  cat("Simulations saved in folder")
  cat("Processing time", tt)

  return(tt)
}