# Run this code after checking convergence and merging chain results

betas <- matrix(NA, nrow = 50, ncol = 50)
z <- matrix(NA, nrow = 50, ncol = 50)

for(sim in 1:50){
    betas[sim,] <- scan(paste0("Simulation2_Gibbs/Sim", sim, "betas.txt"))
    z[sim,] <- scan(paste0("Simulation2_Gibbs/Sim", sim, "Z.txt"))
}

tempo <- matrix(NA, nrow = 2, ncol = 50)
for(chain in 1:2){
  tempo[chain,] <- scan(paste0("Simulation2_Gibbs/chain", chain, "/tempo.txt"))
}

rowSums(tempo)
colSums(tempo)

betas_global <- matrix(NA, 50, 10)
z_global <- matrix(NA, 50, 10)
for(k in 1:10){
  betas_global[,k] <- apply(betas[,grep(paste0("^beta_", k, "_[0-9]+$"),
                                             colnames(res1$pars_mat))], 1, mean)
  z_global[, k] <- apply(z[,grep(paste0("^beta_", k, "_[0-9]+$"),
                                 colnames(res1$pars_mat))], 1, Mode)
} 


#bases = rep(1:10, 5)
#coef = rep(coef, 5)



# plot with same format as pedro
boxplot(t(t(betas_global)*t(z_global)), horizontal = TRUE, yaxt="n")

axis(2, 1:10, labels=c(expression(xi[1]),
                       expression(xi[2]),
                       expression(xi[3]),
                       expression(xi[4]),
                       expression(xi[5]),
                       expression(xi[6]),
                       expression(xi[7]),
                       expression(xi[8]), 
                       expression(xi[9]),
                       expression(xi[10])), cex.axis=1.2)
abline(v = c(-10.4,  -3.3,  -2.0,  -1.8,   0.0 , 2.4,   5.5,   8.4), lty = 2)
#c( -2.9, -0.9, -0.6, -0.5,  0.0, 0.7,  1.6,  2.4)
#c(-10.4,  -3.3,  -2.0,  -1.8,   0.0 , 2.4,   5.5,   8.4)


Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


# mean curve and credible interval 
beta_final <- apply(betas_global, 2, mean)
mode_z_gibbs <- apply(z_global, 2, Mode)





Xt <- seq(0, 1, length = 100)
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt), norder = 4, nbasis = 10)
B <- getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)


# table with estimates

data.frame(Coefficients = 1:10, True = coef, Estimate = as.vector(round(beta_final*mode_z_gibbs, 4)), SE = round(apply(t(t(betas_global)*t(z_global)), 2, sd),4))


# credible band for MCMC
gi_t_MCMC <- t(t(betas_global)*t(z_global))%*%t(B) # 50x100


ETI_MCMC <- apply(gi_t_MCMC, 2, function(x){quantile(x, c(0.025, 0.975))})

plot(x = seq(0, 1, length = 100), y = as.numeric(as.vector(beta_final*mode_z_gibbs)%*% t(B)), lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t))
lines(x = seq(0, 1, length = 100), y = coef%*%t(B), lwd = 1, lty = 4)
lines(x = seq(0, 1, length = 100), y = ETI_MCMC[1,], lwd = 1, lty = 2)
lines(x = seq(0, 1, length = 100), y = ETI_MCMC[2,], lwd = 1, lty = 2)
legend("topright", lwd = 1, lty = c(4,1,2), col = c("black", "red", "black"), legend= c("True curve", "Estimated curve", "Credible bands"))


# comparing VB and MCMC (in terms of mean curve) for this to work you need to run the code in save_res.R first
plot(x = seq(0, 1, length = 100), y = as.numeric(as.vector(beta_final*mode_z_gibbs)%*% t(B)), lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t))
lines(x = seq(0, 1, length = 100), y = as.numeric(as.vector(mu_beta_final*mode_z)%*% t(B)), lwd = 2, col= "blue", type = "l", ylab = expression(g[t]), xlab = expression(t))
legend("topright", lwd = 1, col = c("red", "blue"), legend= c("MCMC", "VB"))

