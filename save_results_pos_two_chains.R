for(sim in 1:100){
  
  # chain1 ----
  #res 1: sample parameters after burnin 50%
  load(paste0("Simulation2_Gibbs_paper 3/chain1/Sim", sim, "IterGibbs_sample.RData"))
  
  out_chain1 <- res1
  
  # chain2 ----
  
  #res 1: sample parameters after burnin 50% and thinning k = 25
  load(paste0("Simulation2_Gibbs_paper 3/chain2/Sim", sim, "IterGibbs_sample.RData"))
  
  out_chain2 <- res1
  
  # generating pos results per sim ----
  
  out <- list(pars = rbind(out_chain1$pars_mat, out_chain2$pars_mat), z_mat = rbind(out_chain1$z_mat, out_chain2$z_mat))
  
  # posterior means and modes
  out2 <- list(pars = apply(out$pars, 2, mean),
               z_mat = apply(out$z_mat, 2, Mode))
  
  out_final <- list(betas = out2$pars[grep("^beta_",names(out2[[1]]))], thetas = out2$pars[grep("^theta_",names(out2[[1]]))], sigma2 = out2$pars[grep("^sigma2",names(out2[[1]]))], tau2 = out2$pars[grep("^tau2",names(out2[[1]]))], Z = out2$z_mat)
  
 
  
  #CI <- lapply(1:out2$Nclusters, function(y){sapply(1:3, function(x){quantile(avalues[[y]][,x], c(0.025, 0.975))})})
  #alpha_sds <- lapply(1:out2$Nclusters, function(x){apply(avalues[[x]], 2, sd)})
  
  for(g in names(out_final)){
    write(out_final[[g]], file = paste0("Simulation2_Gibbs_paper 3", "/", "Sim", sim,  g, ".txt"), ncolumns = length(out_final[[g]]))
  }
  
  cat("Posterior results for simulation", sim,"were saved", "\n")
  
}

