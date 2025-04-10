elbo_res <- c()

case <- 1
nsim <- 1:100
scenario <-  "/Users/carol/Documents/Phd STATISTICS - Western/Research/VB-Bases-Selection/Results/Simulations results used in paper/Simulation3_corr_w"
K <- 10
m <- 5

output_plot <- "Simulation1_1"

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
  # j <- 1
  # for(i in seq(1, 50, 10)){
  #   sigma_values[[j]] = out$Sigma_beta[,i:(i+9)] # previous way of saving Sigma matrices in one matrix
  #   j <- j+1
  # }
  res[['sigma2_beta']][[sim]] <- sapply(1:5, function(x){diag(out$Sigma_beta[,,x])})
  #res[['sigma2_beta']][[sim]] <- sapply(1:5, function(x){diag(sigma_values[[x]])})
  res[['conv_elbo']][sim] <- out$conv_elbo
  res[['delta2']][sim] <- out$delta2
  res[['delta1']][sim] <- out$delta1
  res[["lambda1"]][sim] <- out$lambda1
  res[["lambda2"]][sim] <- out$lambda2
  res[['Sigma2_beta']][[sim]] <- out$Sigma_beta  # out$Sigma_beta #sigma_values #
  res[["w"]][sim] <- out$w
  res[['y']][[sim]] <- out$data
}

mu_beta_values <- do.call("rbind", res[['mu']])
p_values <- do.call("rbind", res[['p']])

# final probabilities across datasets
z_values <- matrix(NA, length(nsim), K*m)
for(k in 1:(K*m)){z_values[,k] = ifelse(p_values[,k] > 0.5, 1, 0)}

# TP and FP per curve
seq_values <- lapply(c(seq(1, m*K, K)), function(x){seq(x,x+K-1)})
TP <- (6*max(nsim))/colSums(sapply(1:m, function(i){apply(z_values[,seq_values[[i]]], 2, sum)}))
FP <- colSums(sapply(1:m, function(i){apply(z_values[,seq_values[[i]]], 2, sum)})[c(2,5,9,10),])/(4*max(nsim))

#per dataset and curve, then compute the mean over datasets (simulated scenario 1 - Bsplines)
seq_values <- lapply(c(seq(1, m*K, K)), function(x){seq(x,x+K-1)})
true <- c(1,0,1,1,0,1,1,1,0,0)
Sens <- sapply(1:m, function(i){mean(rowSums(z_values[,seq_values[[i]]][,c(1,3,4,6,7,8)])/6)})
Spec <- sapply(1:m, function(i){mean(1 - rowSums(z_values[,seq_values[[i]]][,c(2,5,9,10)])/4)})
Acc <- sapply(1:m, function(i){mean(rowSums(z_values[,seq_values[[i]]] == matrix(rep(true,100), ncol = 10, byrow = T))/10)})

#per dataset and curve, then compute the mean over datasets (simulated scenario 3 - Fourier)
seq_values <- lapply(c(seq(1, m*K, K)), function(x){seq(x,x+K-1)})
true <- c(0,1,1,0,0,0,0,0,0,0)
Sens <- sapply(1:m, function(i){mean(rowSums(z_values[,seq_values[[i]]][,c(2,3)])/2)})
Spec <- sapply(1:m, function(i){mean(1 - rowSums(z_values[,seq_values[[i]]][,c(1,4,5,6,7,8,9,10)])/8)})
Acc <- sapply(1:m, function(i){mean(rowSums(z_values[,seq_values[[i]]] == matrix(rep(true,100), ncol = 10, byrow = T))/10)})

# computing etak
eta <- matrix(NA, length(nsim), K)

etakis <- (mu_beta_values*z_values)
for(k in 1:K){
  eta[,k] <- apply(etakis[,grep(paste0("^beta_", k, "_[0-9]+$"),
                                             colnames(mu_beta_values))], 1, mean)
} 

apply(eta, 2, mean)
apply(eta, 2, sd)

z_final <- matrix(NA, length(nsim), K)
for(k in 1:K){
  z_final[,k] <- apply(z_values[,grep(paste0("^beta_", k, "_[0-9]+$"),
                                colnames(mu_beta_values))], 1, Mode)
} 

#per dataset and mean curve, then compute the mean over the datasets(simulated scenario 1 - Bsplines)
true <- c(1,0,1,1,0,1,1,1,0,0)
Sens <- mean(rowSums(z_final[,c(1,3,4,6,7,8)])/6)
Spec <- mean(1 - rowSums(z_final[,c(2,5,9,10)])/4)
Acc <- mean(rowSums(z_final == matrix(rep(true,100), ncol = 10, byrow = T))/10)

#per dataset and mean curve curve, then compute the mean over datasets (simulated scenario 3 - Fourier)
true <- c(0,1,1,0,0,0,0,0,0,0)
Sens <- mean(rowSums(z_final[,c(2,3)])/2)
Spec <- mean(1 - rowSums(z_final[,c(1,4,5,6,7,8,9,10)])/8)
Acc <- mean(rowSums(z_final == matrix(rep(true,100), ncol = 10, byrow = T))/10)

# plot with same format as Pedro's paper:
{
pdf(file = sprintf("%s_%d.pdf", output_plot, 1), width = 6.2, height = 7) #width = 6.25, height = 7.2
par(mfrow=c(1,1), mar = c(5, 4.1, 2, 0.1), cex = 1) 
boxplot(eta, horizontal = TRUE, yaxt="n") #, ylim = c(-2, 2) for fourier
axis(2, 1:K, labels=c(expression(hat(xi)[1]),
                       expression(hat(xi)[2]),
                       expression(hat(xi)[3]),
                       expression(hat(xi)[4]),
                       expression(hat(xi)[5]),
                       expression(hat(xi)[6]),
                       expression(hat(xi)[7]),
                       expression(hat(xi)[8]), 
                       expression(hat(xi)[9]),
                       expression(hat(xi)[10]),
                       expression(hat(xi)[11]),
                       expression(hat(xi)[12]),
                       expression(hat(xi)[13]),
                       expression(hat(xi)[14]),
                       expression(hat(xi)[15])), cex.axis=1.2)
#abline(v = coef, lty = 2)
abline(v = seq(-2,2, by = 0.5), lty = 2)
dev.off()
}

# Boxplot of w
boxplot(res$w)
mean(res$w)
median(res$w)
quantile(res$w, c(0.25,0.75))
sd(res$w)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


#b-splines
Xt <- seq(0, 1, length = 100)
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt), nbasis = K)
B <- getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)

#fourier
Xt <- seq(0, 2*pi, length = 100)
basisBspline_Simulated_Data <- create.fourier.basis(range(Xt), nbasis = K)
B <- getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)

# table with estimates
# bsliples
data.frame(Coefficients = 1:K, True = coef, Estimate = round(apply(eta, 2, mean), 4), SD = round(apply(eta, 2, sd),4))

#fourier
data.frame(Coefficients = 1:K, Estimate = round(apply(eta, 2, mean), 4), SD = round(apply(eta, 2, sd),4))


# sigma2
boxplot(res$delta2/(res$delta1 - 1))

# tau2
boxplot(res$lambda2/(res$lambda1 - 1))


# credible band ----

# Generating the credible band in VB (matches with Pedro's apprach in MCMC) ----
index_values <- lapply(c(seq(1, m*K, K)), function(x){seq(x,x+K-1)})

# selecting a dataset to display the curve
id_sim <- 20 

LL <- NULL
UL <- NULL
estimates <- NULL
for(i in 1:m){
  betas_sample <- MASS::mvrnorm(200, mu = res$mu[[id_sim]][index_values[[i]]],  Sigma = res$Sigma2_beta[[id_sim]][,,i])
  z_sample <- matrix(rbinom(K*200, 1, prob = res$p[[id_sim]][index_values[[i]]]), ncol = K, byrow = TRUE)
  
  estimates <- cbind(estimates, (betas_sample*z_sample))
}

eta_plot <- matrix(NA, 200, K)
for(k in 1:K){
  eta_plot[,k] <- apply(estimates[,grep(paste0("^beta_", k, "_[0-9]+$"),
                                colnames(mu_beta_values))], 1, mean)
}

apply(eta_plot, 2, mean)


curve_f <- list()
for(s in 1:200){
  curve_f[[s]] <- tapply(estimates[s,], rep(1:K, m), mean)%*%t(B)
}  

LL <- apply(do.call(rbind,curve_f), 2, function(i)quantile(i, probs = c(0.025, 0.975)))[1,]
UL <-  apply(do.call(rbind,curve_f), 2, function(i)quantile(i, probs = c(0.025, 0.975)))[2,]


library(scales)

{
  pdf(file = sprintf("%s_%d.pdf", output_plot, 2), width = 6.2, height = 7)
  par(mar = c(5, 4, 2, 0.1) + 0.1) 
  #fourier: B = B[,-1] remove intercept
  plot(x = Xt, y = c(cos(Xt) + sin(2*Xt)), lwd = 2, type = "l", ylab = expression(g[t]), xlab = expression(t), ylim = c(-3, 3), col = "blue") 
  #Bsplines
  #plot(x = Xt, y = coef%*%t(B), lwd = 2, type = "l", ylab = expression(g[t]), xlab = expression(t), ylim = c(-3, 3), col = "blue") 
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


