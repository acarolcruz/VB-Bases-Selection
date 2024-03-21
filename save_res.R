elbo_res <- matrix(NA, nrow = 50, ncol = 10)

for(sim in 1:50){
  for(case in 1:10){
    load(paste0("Simulation1p_info2/case",case,"/Sim", sim, "VBresults.RData"))
    elbo_res[sim,case] <- out$conv_elbo
  }
}

tempo <- matrix(NA, nrow = 50, ncol = 10)
for(case in 1:10){
  tempo[,case] <- scan(paste0("Simulation1p_info2/case",case,"/", "tempo.txt"))
}

rowSums(tempo)
colSums(tempo)

elbo_optimum <- apply(elbo_res, 1, function(x){which.max(x)})


res <- list()
res[['mu']] <- list()
res[['p']] <- list()
res[['sigma2_beta']] <- list() #just the diagonal of the Sigma_beta_i
res[['Sigma2_beta']] <- list() #whole matrix
res[["conv_elbo"]] <- c()
res[["delta2"]] <- c()
res[["delta1"]] <- c()
sigma_values <- list()

for(sim in 1:50){
  
  load(paste0("Simulation1p_info2/case",elbo_optimum[sim],"/Sim", sim, "VBresults.RData"))
  
  res[['mu']][[sim]] <- out$mu_beta
  res[['p']][[sim]] <- out$p
  
  j <- 1
  for(i in seq(1, 50, 10)){
    sigma_values[[j]] = out$Sigma_beta[,i:(i+9)]
    j <- j+1
  }
  res[['sigma2_beta']][[sim]] <- sapply(1:5, function(x){diag(sigma_values[[x]])})
  res[['conv_elbo']][sim] <- out$conv_elbo
  res[['delta2']][sim] <- out$delta2
  res[['delta1']][sim] <- out$delta1
  res[['Sigma2_beta']][[sim]] <- out$Sigma_beta
}

mu_beta_values <- do.call("rbind", res[['mu']])
p_values <- do.call("rbind", res[['p']])
#sigma2_beta_values <- do.call("rbind", res[['sigma2_beta']])

mu_betas <- matrix(NA, 50, 10)
p <- matrix(NA, 50, 10)
for(k in 1:10){
  mu_betas[,k] <- apply(mu_beta_values[,grep(paste0("^beta_", k, "_[0-9]+$"),
                                             colnames(mu_beta_values))], 1, mean)
  p[,k] <- apply(p_values[,grep(paste0("^beta_", k, "_[0-9]+$"),
                                colnames(mu_beta_values))], 1, mean)
  
} 


# final probabilities across datasets 
#p_star <- apply(p, 2, mean)

z_star <- matrix(NA, 50, 10)
for(k in 1:10){z_star[,k] = ifelse(p[,k] > 0.5, 1, 0)}

bases = rep(1:10, 5)
coef = rep(coef, 5)

{
  pdf(file = "output_plots_sim1_across_cases.pdf", width = 15, height = 7)
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
abline(v = c( -2.9, -0.9, -0.6, -0.5,  0.0, 0.7,  1.6,  2.4), lty = 2)
#c( -2.9, -0.9, -0.6, -0.5,  0.0, 0.7,  1.6,  2.4)
#c(-0.6, -0.5, -0.2, 0, 0.1, 0.7, 1.1, 1.8 )


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



plot(x = seq(0, 1, length = 100), y = as.numeric(as.vector(mu_beta_final*mode_z)%*% t(B)), lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t))


# sigma2
boxplot(res$delta2/(res$delta1 - 1), ylim = c(0.01, 0.019))
abline(a = 0.01, b = 0, lwd = 2, col = "red", lty = 2)



# credible band ----

gi_t <- t(t(mu_betas)*t(z_star))%*%t(B) # 50x100


ETI <- apply(gi_t, 2, function(x){quantile(x, c(0.025, 0.975))})

plot(x = seq(0, 1, length = 100), y = as.numeric(as.vector(mu_beta_final*mode_z)%*% t(B)), lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t))
lines(x = seq(0, 1, length = 100), y = ETI[1,], lwd = 1, lty = 2)
lines(x = seq(0, 1, length = 100), y = ETI[2,], lwd = 1, lty = 2)
legend("topright", lwd = 1, col = c("red", "black"), legend= c("Estimated curve", "Credible bands"))



# plot with some points (from one dataset)

data_y <- do.call("rbind", data[['y']])
as.numeric(data_y)
plot(x = rep(seq(0, 1, length = 100), each = 5), y = as.numeric(data_y), xlab = "t", ylab = expression(g(t)))
lines(x = seq(0, 1, length = 100), y = as.numeric(mu_betas[1,]%*% t(B)), lwd = 2, col= "red")


# alternative way of generating the credible band in VB ----

for(sim in 1:50){
  betas_sample <- lapply(c(1, 11, 21, 31, 41), function(x){MASS::mvrnorm(100, mu = mu_beta_values[sim,x:(x+9)],  Sigma = res$Sigma2_beta[[sim]][,x:(x+9)])})
  z_sample <- lapply(c(1, 11, 21, 31, 41), function(x){matrix(rbinom(10*100, 1, prob = p_values[sim, x:(x+9)]), ncol = 10, byrow = TRUE)})
  betas_s_merge <- do.call("rbind", betas_sample)
  z_s_merge <- do.call("rbind", z_sample)
}

LL <- apply((betas_s_merge*z_s_merge)%*%t(B) , 2, function(y){quantile(y, c(0.025,0.975))})[1,]
UL <- apply((betas_s_merge*z_s_merge)%*%t(B) , 2, function(y){quantile(y, c(0.025,0.975))})[2,]

plot(x = seq(0, 1, length = 100), y = as.numeric(as.vector(mu_beta_final*mode_z)%*% t(B)), lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t))
lines(x = seq(0, 1, length = 100), y = LL, lwd = 1, lty = 2)
lines(x = seq(0, 1, length = 100), y = UL, lwd = 1, lty = 2)
legend("topright", lwd = 1, col = c("red", "black"), legend= c("Estimated curve", "Credible bands"))


