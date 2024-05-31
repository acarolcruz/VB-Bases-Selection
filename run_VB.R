# Run VB algorithm for a simulation and save results in a folder

run_VB <- function(sim = sim, save_folder = save_folder, scenario = scenario, m = m, times = times, ordem = ordem, K = K, y = y, B = B, coef = coef, seed = seed, 
                         sigma = sigma, mu_ki = mu_ki, lambda_1 = lambda_1, lambda_2 = lambda_2, delta_1 = delta_1, delta_2 = delta_2, 
                         maxIter = maxIter, initial_values = initial_values, convergence_threshold = convergence_threshold, w = w, corr = corr, basis_type = basis_type){
  
  tt <- Sys.time()
  
  # Generate the data for one dataset 
  
  if(!is.null(w) & corr == TRUE){
    source("sim_data_corr.R")
    source("elbo_formulas_corr.R")
    source("aux_fun_corr.R")
    source("VB_corr_vem.R")
    data <- sim_data_corr(m = m, times = times, ordem = ordem, K = K, coef = coef, seed = seed + sim, sigma = sigma, w = w, basis_type = basis_type)
    out <- vb_bs_corr(y = data$y, B = data$B, m = m, mu_ki = mu_ki, lambda_1 = lambda_1, lambda_2 = lambda_2, delta_1 = delta_1, delta_2 = delta_2, maxIter = maxIter, K = K, initial_values = initial_values, convergence_threshold = convergence_threshold, Xt = data$Xt)
    
  } else{
    if(!is.null(w) & corr == FALSE){
      source("sim_data_corr.R")
      source("elbo_formulas_corr.R")
      source("aux_fun_corr.R")
      source("vb_bs_corr_fixw.R")
      data <- sim_data_corr(m = m, times = times, ordem = ordem, K = K, coef = coef, seed = seed + sim, sigma = sigma, w = w, basis_type = basis_type)
      out <- vb_bs_corr_fixw(y = data$y, B = data$B, m = m, mu_ki = mu_ki, lambda_1 = lambda_1, lambda_2 = lambda_2, delta_1 = delta_1, delta_2 = delta_2, maxIter = maxIter, K = K, initial_values = initial_values, convergence_threshold = convergence_threshold, Xt = data$Xt, w = w)
      }
    if(!is.null(w) & corr == "miss"){
      source("sim_data_corr.R")
      source("expected_values.R")
      source("elbo_formulas.R")
      source('elbo.R')
      source('VB_bs_main.R')
      data <- sim_data_corr(m = m, times = times, ordem = ordem, K = K, coef = coef, seed = seed + sim, sigma = sigma, w = w, basis_type = basis_type)
      out <- vb_bs(y = data$y, B = data$B, m = m, mu_ki = mu_ki, lambda_1 = lambda_1, lambda_2 = lambda_2, delta_1 = delta_1, delta_2 = delta_2, maxIter = maxIter, K = K, initial_values = initial_values, convergence_threshold = convergence_threshold)
      
      
    }else{
        source("sim_data.R")
        source("expected_values.R")
        source("elbo_formulas.R")
        source('elbo.R')
        source('VB_bs_main.R')
        data <- sim_data(m = m, times = times, ordem = ordem, K = K, coef = coef, seed = seed + sim, sigma = sigma)
        out <- vb_bs(y = data$y, B = data$B, m = m, mu_ki = mu_ki, lambda_1 = lambda_1, lambda_2 = lambda_2, delta_1 = delta_1, delta_2 = delta_2, maxIter = maxIter, K = K, initial_values = initial_values, convergence_threshold = convergence_threshold)
        }
    }
  
  # Run VB basis selection algorithm
  #out <- vb_bs(y = data$y, B = data$B, m = m, mu_ki = mu_ki, lambda_1 = lambda_1, lambda_2 = lambda_2, delta_1 = delta_1, delta_2 = delta_2, maxIter = maxIter, K = K, initial_values = initial_values, convergence_threshold = convergence_threshold)
  
  tt <- difftime(Sys.time(), tt, units = 'mins')
  
  save(out, file = paste0(save_folder, "/case", scenario, "/Sim", sim, "VBresults.RData"))
  
  # for(param in names(data)){
  #   if(is.matrix(data[[param]])){
  #     write(data[[param]], file = paste0(save_folder, "/case", scenario, "/Sim", sim, param, ".txt"), ncolumns = ncol(data[[param]]))
  #   } else{
  #     write(data[[param]], file = paste0(save_folder, "/case", scenario, "/Sim", sim, param, ".txt"), ncolumns = length(data[[param]]))
  #   }
  # }
  
  return(tt)
}
