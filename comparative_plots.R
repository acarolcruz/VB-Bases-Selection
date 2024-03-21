# sigma small (Simulation 1)----

{
  pdf(file = "output_plots_sim_sigma_small_1.pdf", width = 12.5, height = 7.2)
  par(mfrow=c(1,2), mar = c(5, 2, 2, 0.1), cex = 1) 
  
  # VB
  boxplot(t(t(mu_betas)*t(z_star)), horizontal = TRUE, yaxt="n", main = "VB")
  
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
  abline(v = c( -2.9, -0.9, -0.6, -0.5,  0.0, 0.7,  1.6,  2.4), lty = 2)
  
  # Gibbs
  boxplot(t(t(betas_global)*t(z_global)), horizontal = TRUE, yaxt="n", main = "Gibbs")
  
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
  abline(v = c( -2.9, -0.9, -0.6, -0.5,  0.0, 0.7,  1.6,  2.4), lty = 2)
  
  # curve 
  # VB
  plot(x = seq(0, 1, length = 100), y = as.numeric(as.vector(mu_beta_final*mode_z)%*% t(B)), lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t), main = "VB")
  lines(x = seq(0, 1, length = 100), y = coef%*%t(B), lwd = 1, lty = 4)
  lines(x = seq(0, 1, length = 100), y = ETI[1,], lwd = 1, lty = 2)
  lines(x = seq(0, 1, length = 100), y = ETI[2,], lwd = 1, lty = 2)
  legend("topright", lwd = 1, lty = c(4, 1, 2), col = c("black", "red", "black"), legend= c("True curve", "Estimated curve", "Credible bands"))
  
  
  
  # Gibbs
  plot(x = seq(0, 1, length = 100), y = as.numeric(as.vector(beta_final*mode_z_gibbs)%*% t(B)), lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t), main = "Gibbs")
  lines(x = seq(0, 1, length = 100), y = coef%*%t(B), lwd = 1, lty = 4)
  lines(x = seq(0, 1, length = 100), y = ETI_MCMC[1,], lwd = 1, lty = 2)
  lines(x = seq(0, 1, length = 100), y = ETI_MCMC[2,], lwd = 1, lty = 2)
  legend("topright", lwd = 1, lty = c(4,1,2), col = c("black", "red", "black"), legend= c("True curve", "Estimated curve", "Credible bands"))
}  
dev.off()



{
  pdf(file = "output_plots_sim_sigma_small_2.pdf", width = 13, height = 7.2)
  par(mfrow=c(1,1), mar = c(5, 2, 2, 0.1), cex = 1) 
  
  # Comparing VB and Gibbs
  plot(x = seq(0, 1, length = 100), y = as.numeric(as.vector(beta_final*mode_z_gibbs)%*% t(B)), lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t))
  lines(x = seq(0, 1, length = 100), y = as.numeric(as.vector(mu_beta_final*mode_z)%*% t(B)), lwd = 2, col= "blue", ylab = expression(g[t]), xlab = expression(t), lty = 2)
  legend("topright", lwd = 1, lty = c(1,2), col = c("red", "blue"), legend= c("MCMC", "VB"))
}  
dev.off()



# sigma large (Simulation 2) ----
{
  pdf(file = "output_plots_sim_sigma_large_1.pdf", width = 13, height = 7.2)
  par(mfrow=c(1,2), mar = c(5, 2, 2, 0.1), cex = 1) 
  
  # VB
  boxplot(t(t(mu_betas)*t(z_star)), horizontal = TRUE, yaxt="n", main = "VB")
  
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
  
  
  
  # Gibbs
  boxplot(t(t(betas_global)*t(z_global)), horizontal = TRUE, yaxt="n", main = "Gibbs")
  
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
  
  
  # curve 
  # VB
  plot(x = seq(0, 1, length = 100), y = as.numeric(as.vector(mu_beta_final*mode_z)%*% t(B)), lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t), main = "VB")
  lines(x = seq(0, 1, length = 100), y = coef%*%t(B), lwd = 1, lty = 4)
  lines(x = seq(0, 1, length = 100), y = ETI[1,], lwd = 1, lty = 2)
  lines(x = seq(0, 1, length = 100), y = ETI[2,], lwd = 1, lty = 2)
  legend("topright", lwd = 1, lty = c(4, 1, 2), col = c("black", "red", "black"), legend= c("True curve", "Estimated curve", "Credible bands"))
  
  
  
  # Gibbs
  plot(x = seq(0, 1, length = 100), y = as.numeric(as.vector(beta_final*mode_z_gibbs)%*% t(B)), lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t))
  lines(x = seq(0, 1, length = 100), y = coef%*%t(B), lwd = 1, lty = 4)
  lines(x = seq(0, 1, length = 100), y = ETI_MCMC[1,], lwd = 1, lty = 2)
  lines(x = seq(0, 1, length = 100), y = ETI_MCMC[2,], lwd = 1, lty = 2)
  legend("topright", lwd = 1, lty = c(4,1,2), col = c("black", "red", "black"), legend= c("True curve", "Estimated curve", "Credible bands"))
}  
dev.off()


{
  pdf(file = "output_plots_sim_sigma_large_2.pdf", width = 13, height = 7.2)
  par(mfrow=c(1,1), mar = c(5, 2, 2, 0.1), cex = 1) 
  
  # Comparing VB and Gibbs
  plot(x = seq(0, 1, length = 100), y = as.numeric(as.vector(beta_final*mode_z_gibbs)%*% t(B)), lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t))
  lines(x = seq(0, 1, length = 100), y = as.numeric(as.vector(mu_beta_final*mode_z)%*% t(B)), lwd = 2, col= "blue", ylab = expression(g[t]), xlab = expression(t), lty = 2)
  legend("topright", lwd = 1, lty = c(1,2), col = c("red", "blue"), legend= c("MCMC", "VB"))
}  
dev.off()
