
library(postpack)
library(coda)


for(sim in 1:50){
    
    cat("Convergence results simulation", sim, "\n")
    
    # chain1 ----
    #res 1: sample parameters after burnin 50%
    load(paste0("Simulation1_Gibbs/chain1/Sim", sim, "IterGibbs_sample.RData"))
    
    out_chain1 <- res1

    # chain2 ----
    
    #res 1: sample parameters after burnin 50% and thinning k = 25
    load(paste0("Simulation1_Gibbs/chain2/Sim", sim, "IterGibbs_sample.RData"))
    
    out_chain2 <- res1
    
    
    
    # convergence for beta
    mcmc <- list(out_chain1$pars_mat[,1:50], out_chain2$pars_mat[,1:50])
    
    results_mcmc <- post_convert(mcmc)
    
    print(gelman.diag(results_mcmc))
    
    # convergence for z - It does not work (is it because z is discrete?)
    mcmc <- list(out_chain1$z_mat, out_chain2$z_mat)

    try(results_mcmc <- post_convert(mcmc))
    try(print(gelman.diag(results_mcmc)))
    
    
    
}


