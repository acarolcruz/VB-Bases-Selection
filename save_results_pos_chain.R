
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


save_results_pos <- function(scenario = "Simulation1_Gibbs_paper", chain = 1, sim = 1, startpoint_sampling = 10, maxIter = 10000){
  
  load(paste0(paste0(scenario, "/chain", chain),"/Sim", sim, "Gibbsresults.RData"))
  
  index <- seq(from = (maxIter/2)+1, to = maxIter, by = startpoint_sampling)
  
  res1 <- list(pars_mat = out[[1]][index,], z_mat = out[[2]][index,])
  
  save(res1, file = paste0(scenario, "/chain", chain, "/Sim", sim, "IterGibbs_sample.RData"))
  
  #computing posterior means and modes for all parameters except alpha and T
  out2 <- list(pars = apply(res1$pars_mat, 2, mean),
               z_final = apply(res1$z_mat, 2, Mode))
  
  out_final <- list(betas = out2$pars[grep("^beta_",names(out2[[1]]))], thetas = out2$pars[grep("^theta_",names(out2[[1]]))], sigma2 = out2$pars[grep("^sigma2",names(out2[[1]]))], tau2 = out2$pars[grep("^tau2",names(out2[[1]]))], zs = out2$z_final)

  #for(g in names(out_final)){
  #    write(out_final[[g]], file = paste0(scenario, "/chain", chain, "/Sim", sim,  g, ".txt"), ncolumns = length(out_final[[g]]))
  #}
  
  cat("Results saved in folder \n")
}

for(sim in 1:100){save_results_pos(scenario = "Simulation2_Gibbs_miss_paper", chain = 2, sim = sim, startpoint_sampling = 50)}
