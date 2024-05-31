elbo_res <- c()

case <- 1
nsim <- 1:100
scenario <-  "Simulation1_VB_gibbs"

output_plot <- "Simulation1_vb_1"

# With no random initialization ----

for(sim in nsim){
  load(paste0(scenario,"/case",case,"/Sim", sim, "VBresults.RData"))
  elbo_res[sim] <- out$conv_elbo
}

tempo <- scan(paste0(scenario,"/case",case,"/", "tempo.txt"))


mean(tempo)


res <- list()
res[['mu']] <- list()
res[['p']] <- list()
res[['sigma2_beta']] <- list() #just the diagonal of the Sigma_beta_i
res[['Sigma2_beta']] <- list() #whole matrix
res[["conv_elbo"]] <- c()
res[["delta2"]] <- c()
res[["delta1"]] <- c()
res[["lambda2"]] <- c()
res[["lambda1"]] <- c()
res[["w"]] <- c()
res[["y"]] <- list()
sigma_values <- list()

for(sim in nsim){
  
  load(paste0(scenario,"/case",case,"/Sim", sim, "VBresults.RData"))
  
  res[['mu']][[sim]] <- out$mu_beta
  res[['p']][[sim]] <- out$p
  j <- 1
  for(i in seq(1, 50, 10)){
    sigma_values[[j]] = out$Sigma_beta[,i:(i+9)] # previous way of saving Sigma matrices in one matrix
    j <- j+1
  }
  #res[['sigma2_beta']][[sim]] <- sapply(1:5, function(x){diag(out$Sigma_beta[,,x])})
  res[['sigma2_beta']][[sim]] <- sapply(1:5, function(x){diag(sigma_values[[x]])})
  res[['conv_elbo']][sim] <- out$conv_elbo
  res[['delta2']][sim] <- out$delta2
  res[['delta1']][sim] <- out$delta1
  res[["lambda1"]][sim] <- out$lambda1
  res[["lambda2"]][sim] <- out$lambda2
  res[['Sigma2_beta']][[sim]] <- sigma_values  # out$Sigma_beta #sigma_values #
  res[["w"]][sim] <- out$w
  res[['y']][[sim]] <- out$data
}

mu_beta_values <- do.call("rbind", res[['mu']])
p_values <- do.call("rbind", res[['p']])

# final probabilities across datasets
z_values <- matrix(NA, length(nsim), 50)
for(k in 1:50){z_values[,k] = ifelse(p_values[,k] > 0.5, 1, 0)}

# computing etak
eta <- matrix(NA, length(nsim), 10)

etakis <- (mu_beta_values*z_values)
for(k in 1:10){
  eta[,k] <- apply(etakis[,grep(paste0("^beta_", k, "_[0-9]+$"),
                                             colnames(mu_beta_values))], 1, mean)
} 

apply(eta, 2, mean)
apply(eta, 2, sd)

# plot with same format as Pedro's paper:
{
pdf(file = sprintf("%s_%d.pdf", output_plot, 1), width = 6.2, height = 7) #width = 6.25, height = 7.2
par(mfrow=c(1,1), mar = c(5, 4.1, 2, 0.1), cex = 1) 
boxplot(eta, horizontal = TRUE, yaxt="n") #, ylim = c(-2, 2) for fourier
axis(2, 1:10, labels=c(expression(hat(xi)[1]),
                       expression(hat(xi)[2]),
                       expression(hat(xi)[3]),
                       expression(hat(xi)[4]),
                       expression(hat(xi)[5]),
                       expression(hat(xi)[6]),
                       expression(hat(xi)[7]),
                       expression(hat(xi)[8]), 
                       expression(hat(xi)[9]),
                       expression(hat(xi)[10])), cex.axis=1.2)
abline(v = coef, lty = 2)
dev.off()
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


#b-splines
Xt <- seq(0, 1, length = 100)
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt), nbasis = 10)
B <- getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)

#fourier
Xt <- seq(0, 2*pi, length = 100)
basisBspline_Simulated_Data <- create.fourier.basis(range(Xt), nbasis = 10)
B <- getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)

# table with estimates
#new
data.frame(Coefficients = 1:10, True = coef, Estimate = round(apply(eta, 2, mean), 4), SD = round(apply(eta, 2, sd),4))

apply(eta, 2, mean)
apply(eta, 2, sd)


# sigma2
boxplot(res$delta2/(res$delta1 - 1))
abline(a = 0.2^2, b = 0, lwd = 2, col = "red", lty = 2)
abline(a = 0.01, b = 0, lwd = 2, col = "red", lty = 2)

# tau2
boxplot(res$lambda2/(res$lambda1 - 1))
abline(a = 700, b = 0, lwd = 2, col = "red", lty = 2)


# credible band ----

# Generating the credible band in VB (matches with Pedro's apprach in MCMC) ----
index_values <- list(seq(1,10), seq(11,20), seq(21,30), seq(31,40), seq(41,50))

# selecting a dataset to display the curve
id_sim = 20 

LL <- NULL
UL <- NULL
estimates <- NULL
for(i in 1:5){
  betas_sample <- MASS::mvrnorm(200, mu = res$mu[[id_sim]][index_values[[i]]],  Sigma = res$Sigma2_beta[[id_sim]][[i]])
  z_sample <- matrix(rbinom(10*200, 1, prob = res$p[[id_sim]][index_values[[i]]]), ncol = 10, byrow = TRUE)
  
  estimates <- cbind(estimates, (betas_sample*z_sample))
}

eta_plot <- matrix(NA, 200, 10)
for(k in 1:10){
  eta_plot[,k] <- apply(estimates[,grep(paste0("^beta_", k, "_[0-9]+$"),
                                colnames(mu_beta_values))], 1, mean)
}

apply(eta_plot, 2, mean)


curve_f <- list()
for(s in 1:200){
  curve_f[[s]] <- tapply(estimates[s,], rep(1:10,5), mean)%*%t(B)
}  

LL <- apply(do.call(rbind,curve_f),2,function(i)quantile(i,probs = c(0.025,0.975)))[1,]
UL <-  apply(do.call(rbind,curve_f),2,function(i)quantile(i,probs = c(0.025,0.975)))[2,]


library(scales)

{
  pdf(file = sprintf("%s_%d.pdf", output_plot, 2), width = 6.2, height = 7)
  par(mar = c(5, 4, 2, 0.1) + 0.1) 
  #fourier: B = B[,-1] remove intercept
  #plot(x = Xt, y = c(cos(Xt) + sin(2*Xt)), lwd = 2, type = "l", ylab = expression(g[t]), xlab = expression(t), ylim = c(-3, 3), col = "blue") 
  #Bsplines
  plot(x = Xt, y = coef%*%t(B), lwd = 2, type = "l", ylab = expression(g[t]), xlab = expression(t), ylim = c(-3, 3), col = "blue") 
  lines(x = Xt, y = apply(eta_plot, 2, mean)%*%t(B), lwd = 2, lty = 2, col= "red")
  # plot points from the curves
  for(i in 1:5){
    points(x = Xt, y = res$y[[id_sim]][[i]], cex = 0.8, pch = 16, col = alpha("black", 0.2))
  }
  lines(x = Xt, y = LL, lwd = 2, lty = 3)
  lines(x = Xt, y = UL, lwd = 2, lty = 3)
  legend("topright", lwd = 2, lty = c(2, 1, 3), col = c("red", "blue", "black"), legend= c("Estimated curve", "True curve", "Credible band"))
    dev.off()
}


