# Read all files and create plots
res <- list()
res[['mu']] <- list()
res[['p']] <- list()
res[['sigma2_beta']] <- list()
res[["conv_elbo"]] <- c()
res[["delta2"]] <- c()

for(sim in 1:50){
  
  load(paste0("Scenario1c2/case1","/Sim", sim, "VBresults.RData"))
  
  res[['mu']][[sim]] <- out$mu_beta
  res[['p']][[sim]] <- out$p
  res[['sigma2_beta']][[sim]] <- diag(out$Sigma_beta) # change how the this information is being stored (use array)
  res[['conv_elbo']][sim] <- out$conv_elbo
  res[['delta2']][sim] <- out$delta2
}

tempo <- scan(paste0("Scenario_p1/case1/tempo.txt"))
mean(tempo) # less than a second to run each dataset 

mu_beta_values <- do.call("rbind", res[['mu']])
p_values <- do.call("rbind", res[['p']])
sigma2_beta_values <- do.call("rbind", res[['sigma2_beta']])

mu_betas <- matrix(NA, 50, 10)
p <- matrix(NA, 50, 10)
for(k in 1:10){
  mu_betas[,k] <- apply(mu_beta_values[,grep(paste0("^beta_", k, "_[0-9]+$"),
                                             colnames(mu_beta_values))], 1, mean)
  p[,k] <- apply(p_values[,grep(paste0("^beta_", k, "_[0-9]+$"),
                                     colnames(mu_beta_values))], 1, mean)
  
}  


# final probabilities across datasets 
p_star <- apply(p, 2, mean)

z_star <- matrix(NA, 50, 10)
for(k in 1:10){z_star[,k] = ifelse(p[,k] > 0.5, 1, 0)}

bases = rep(1:10, 5)
coef = rep(coef, 5)

{
  pdf(file = "output_plots_sim1_m_5.pdf", width = 15, height = 7)
  par(mfrow=c(1,10), mar = c(5, 2, 2, 0.1), cex = 0.5) 
  for(k in 1:10){
    boxplot(mu_betas[,k]*z_star[,k], xlab = paste0("beta for base ", bases[k]), ylim = c(-2.5, 1.5))
    abline(a = coef[k], b = 0, lwd = 2, col = "red", lty = 2)
  }
  
}
dev.off()


# plot with same format as pedro
boxplot(t(t(mu_betas)*t(z_star)), horizontal = TRUE, yaxt="n")

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
abline(v = seq(-2,1.5,0.5), lty = 2)


Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


# mean curve and credible interval 
mu_beta_final <- apply(mu_betas, 2, mean)
mode_z <- apply(z_star, 2, Mode)


Xt <- seq(0, 1, length = 100)
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt), norder = 4, nbasis = 10)
B <- getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)



plot(x = seq(0, 1, length = 100), y = as.numeric(as.vector(mu_beta_final*mode_z)%*% t(B)), lwd = 2, col= "red", type = "l")




data_y <- do.call("rbind", data[['y']])

as.numeric(data_y)








plot(x = rep(seq(0, 1, length = 100), each = 5), y = as.numeric(data_y), xlab = "t", ylab = expression(g(t)))
lines(x = seq(0, 1, length = 100), y = as.numeric(mu_betas[1,]%*% t(B)), lwd = 2, col= "red")


list(expression(xi_1),
     expression(xi_2),
     expression(xi_3),
     expression(xi_4),
     expression(xi_5),
     expression(xi_6),
     expression(xi_7),
     expression(xi_8), 
     expression(xi_9),
     expression(xi_10))
     
expression(eta_k)