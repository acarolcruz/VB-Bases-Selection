# Real data application
source("sim_data_corr.R")
source("elbo_formulas_corr.R")
source("aux_fun_corr.R")
source("VB_corr_vem.R")

# Motorcycle data
library(monomvn)
library(glmnet)
library(faraway)

data(mcycle)
plot(accel~jitter(mcycle$time),mcycle, ylab = "Head acceleration (in units of g)", xlab = "time in ms")
y = list(mcycle$accel)
Xt <- jitter(mcycle$times)

# Regression splines:
regression_splines <- lm(y[[1]] ~ B[[1]] - 1)
summary(regression_splines) #MSE = 23.16**2 (536.39)

#bs(Xt,20)

K <- 30
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt),norder = 4, nbasis = K)
B <- lapply(1, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

initial_values <- list(p = rep(1, K), delta2 = 95484.5, lambda2 = 1000, w = 10)
system.time(out_mcycle_K30 <- vb_bs_corr(y = y, B = B, m = 1, mu_ki = 0.5, lambda_1 = 1e-6, lambda_2 = 1e-6, delta_1 = 100, delta_2 = 99*30, maxIter = 1000, K = K, initial_values, convergence_threshold = 0.001, Xt = Xt))

#   user  system elapsed 
# 71.736   0.428  72.643 

# Compute R2 vem
seq_values <- lapply(c(seq(1, 1*K, K)), function(x){seq(x,x+K-1)})
yhat <- as.numeric((out_mcycle_K30$mu_beta[seq_values[[i]]]*ifelse(out_mcycle_K30$p[seq_values[[i]]] > 0.50, 1, 0))%*%t(B[[i]]))
RSS_vb = sum((yhat- y[[1]])^2)
TSS = sum((y[[1]] - mean(y[[1]]))^2) 
R2a_vb = 1 - (RSS_vb/(length(y[[1]])-sum(ifelse(out_mcycle_K30$p[seq_values[[i]]] > 0.50, 1, 0))))/(TSS/(length(y[[1]])-1))
R2a_vb
#0.7936501

K <- 20
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt),norder = 4, nbasis = K)
B <- lapply(1, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

initial_values <- list(p = rep(1, K), delta2 = 40000, lambda2 = 100, w = 10)#91786.5, before: 92000 delta1 = 100, delta2=99*30
system.time(out_mcycle_K20 <- vb_bs_corr(y = y, B = B, m = 1, mu_ki = 0.5, lambda_1 = 1e-6, lambda_2 = 1e-6, delta_1 = 0.1, delta_2 = 0.1, maxIter = 1000, K = K, initial_values, convergence_threshold = 0.001, Xt = Xt))
#user  system elapsed 
# 47.745   0.137  48.066

seq_values <- lapply(c(seq(1, 1*K, K)), function(x){seq(x,x+K-1)})
yhat <- as.numeric((out_mcycle_K20$mu_beta[seq_values[[i]]]*ifelse(out_mcycle_K20$p[seq_values[[i]]] > 0.50, 1, 0))%*%t(B[[i]]))
RSS_vb = sum((yhat- y[[1]])^2)
TSS = sum((y[[1]] - mean(y[[1]]))^2) 
R2a_vb = 1 - (RSS_vb/(length(y[[1]])-sum(ifelse(out_mcycle_K20$p[seq_values[[i]]] > 0.50, 1, 0))))/(TSS/(length(y[[1]])-1))
#0.78729

K <- 15
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt),norder = 4, nbasis = K)
B <- lapply(1, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

initial_values <- list(p = rep(1, K), delta2 = 91517, lambda2 = 100, w = 10)
system.time(out_mcycle_K15 <- vb_bs_corr(y = y, B = B, m = 1, mu_ki = 0.5, lambda_1 = 1e-6, lambda_2 = 1e-6, delta_1 = 100, delta_2 = 99*30, maxIter = 1000, K = K, initial_values, convergence_threshold = 0.001, Xt = Xt))

#   user  system elapsed 
# 39.756   0.160  40.021 

seq_values <- lapply(c(seq(1, 1*K, K)), function(x){seq(x,x+K-1)})
yhat <- as.numeric((out_mcycle_K15$mu_beta[seq_values[[i]]]*ifelse(out_mcycle_K15$p[seq_values[[i]]] > 0.50, 1, 0))%*%t(B[[i]]))
RSS_vb = sum((yhat- y[[1]])^2)
TSS = sum((y[[1]] - mean(y[[1]]))^2) 
R2a_vb = 1 - (RSS_vb/(length(y[[1]])-sum(ifelse(out_mcycle_K15$p[seq_values[[i]]] > 0.50, 1, 0))))/(TSS/(length(y[[1]])-1))
#0.7854

# independent case ----
source("expected_values.R")
source("elbo_formulas.R")
source('elbo.R')
source('VB_bs_main.R')
system.time(out_mcycle_K20_ind <- vb_bs(y = y, B = B, m = 1, mu_ki = 0.5, lambda_1 = 1e-6, lambda_2 = 1e-6, delta_1 = 100, delta_2 = 99*30, maxIter = 1000, K = K, initial_values, convergence_threshold = 0.001))
# 0.7875
#   user  system elapsed 
# 1.557   0.050   1.634 

seq_values <- lapply(c(seq(1, 1*K, K)), function(x){seq(x,x+K-1)})
yhat <- as.numeric((out_mcycle_K20_ind$mu_beta[seq_values[[i]]]*ifelse(out_mcycle_K20_ind$p[seq_values[[i]]] > 0.50, 1, 0))%*%t(B[[i]]))
RSS_vb = sum((yhat- y[[1]])^2)
TSS = sum((y[[1]] - mean(y[[1]]))^2) 
R2a_vb = 1 - (RSS_vb/(length(y[[1]])-sum(ifelse(out_mcycle_K20_ind$p[seq_values[[i]]] > 0.50, 1, 0))))/(TSS/(length(y[[1]])-1))

# Bayesian lasso ----
system.time(BL <- blasso(X = B[[1]], y = y[[1]], T = 10000, thin = 1, 
                         beta=rep(-1,1*20),lambda2=1,s2=1,
                         RJ = F, rao.s2=F, icept = F, normalize=F))

# user  system elapsed 
# 0.275   0.009   0.406 

beta_BL <-  apply(BL$beta,2, mean)
estimativa_funcionalBL <- B[[1]]%*%beta_BL

RSS <- sum((y[[1]]-estimativa_funcionalBL)^2)
TSS <- sum((y[[1]] - mean(y[[1]]))^2)
K_model <- sum(abs(beta_BL)>0)
R2_adj_BL <- 1 - mean((RSS/(length(y[[1]])-K_model))/(TSS/(length(y[[1]]-1))))
# 0.7457

# LASSO ----
regularizador=seq(0.001,1,length=9)

result_LASSO <- NULL
R2_adj_LASSO_mat <- c()
time_lasso <- c()
for(lambda in 1:length(regularizador)){
  #setTxtProgressBar(pb,lambda)
  
  lambdaLIVROtexto <- regularizador[lambda] #parametrização de livro texto.
  time_lasso <- c(time_lasso, system.time(LASSO <- glmnet(B[[1]], y[[1]], alpha = 1,lambda = lambdaLIVROtexto/(2*nrow(B[[1]])), standardize = F, intercept = F)))
  result_LASSO[[lambda]] <- LASSO
  betaLASSO <- as.numeric(coef(LASSO)[-1])
  estimativa_funcional_lasso <- B[[1]]%*%betaLASSO
  
  RSS=sum((y[[1]] - estimativa_funcional)^2)
  TSS=sum((y[[1]] - mean(y[[1]]))^2)
  K_model <- sum(abs(betaLASSO)>0)
  R2_adj_LASSO_mat[lambda] <- 1-mean((RSS/(length(y[[1]])-K_model))/(TSS/(length(y[[1]])-1)))
}
mean(R2_adj_LASSO_mat)
lambda = which.max(R2_adj_LASSO_mat)
# 0.7715

# to plot LASSO's results
LASSO <- glmnet(B[[1]], y[[1]], alpha = 1,lambda = lambdaLIVROtexto/(2*nrow(B[[1]])), standardize = F, intercept = F)
betaLASSO <- as.numeric(coef(LASSO)[-1])
estimativa_funcional_lasso <- B[[1]]%*%betaLASSO

# time 
mean(time_lasso[seq(from = 3, to = 45, 5)])
# 0.001444444

# Credible bands - old way
LL = sapply(1:20, function(i){qnorm(0.025, mean = as.numeric(out_mcycle_K20$mu_beta[i]), sd = sqrt(diag(out_mcycle_K20$Sigma_beta[,,1])[i]))*qbinom(0.025, 1, prob = out_mcycle_K20$p[i])})%*%t(B[[1]])

UL = sapply(1:20, function(i){qnorm(0.975, mean = as.numeric(out_mcycle_K20$mu_beta[i]), sd = sqrt(diag(out_mcycle_K20$Sigma_beta[,,1])[i]))*qbinom(0.025, 1, prob = out_mcycle_K20$p[i])})%*%t(B[[1]])

LL <- NULL
UL <- NULL
betas_sample <- MASS::mvrnorm(200, mu = out_mcycle_K20$mu_beta,  Sigma = out_mcycle_K20$Sigma_beta[,,1])
z_sample <- matrix(rbinom(20*200, 1, prob = out_mcycle_K20$p), ncol = 20, byrow = TRUE)
curve <- (betas_sample*z_sample)%*%t(B[[1]])
LL <- apply(curve, 2, function(x){quantile(x, 0.25)})
UL <- apply(curve, 2, function(x){quantile(x, 0.975)})

# Credible band as in Pedro's article ----
LL <- NULL
UL <- NULL
estimates <- NULL

betas_sample <- MASS::mvrnorm(200, mu = out_mcycle_K20$mu_beta,  Sigma = out_mcycle_K20$Sigma_beta[,,1])
z_sample <- matrix(rbinom(20*200, 1, prob =out_mcycle_K20$p), ncol = 20, byrow = TRUE)

estimates <- cbind(estimates, (betas_sample*z_sample))

curve_f <- list()
for(s in 1:200){
  curve_f[[s]] <- estimates[s,]%*%t(B[[1]])
}  

LL <- apply(do.call(rbind,curve_f),2,function(i)quantile(i,probs = c(0.025,0.975)))[1,]
UL <-  apply(do.call(rbind,curve_f),2,function(i)quantile(i,probs = c(0.025,0.975)))[2,]


#finer grid
Xt_finer <- seq(2.43,57.63,length = 200)

K <- 20
basisbspline_finer <- create.bspline.basis(range(Xt_finer),norder = 4, nbasis = K)
B_finer <- lapply(1, function(x){getbasismatrix(Xt_finer, basisbspline_finer, nderiv = 0)})


# comparing VEM with smoothing splines and regression splines
{
  pdf(file = sprintf("%s_%d.pdf", "plot_real_data", 1), width = 12.5, height = 7.2)
  seq_values <- lapply(c(seq(1, 1*K, K)), function(x){seq(x,x+K-1)})
  par(mfrow=c(1,1), mar = c(5, 4.1, 2, 0.1), cex = 1) 
  for(i in 1:1){
    plot(y[[i]]~Xt,  cex = 0.8, pch = 16, col = alpha("black", 0.2), ylab = expression(g[t]), xlab = expression(t))
    lines(Xt_finer, y = (out_mcycle_K20$mu_beta*ifelse(out_mcycle_K20$p > 0.50, 1, 0))%*%t(B_finer[[i]]), lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t))
    #lines(predict(smooth.spline(Xt, y[[1]]), Xt))
    lines(Xt, fitted(regression_splines), lwd = 2, col = "blue")
    lines(predict(smooth.spline(Xt, y[[1]]), Xt), lwd = 2, col = "orange")
    lines(Xt, LL, col = "black", lty = 2)
    lines(Xt, UL, col = "black", lty = 2)
    legend("topright", legend = c("Proposed method", "Regression splines", "Smoothing splines", "Credible Band"), col = c("red", "blue", "orange", "black"), lty = c(1, 1, 1, 2), bty = "n")
    print(ifelse(out_mcycle_K20$p[seq_values[[i]]] > 0.50, 1, 0))
  }
}
dev.off()

# comparing VEM with Lasso and Bayesian Lasso
library(scales)
{
  pdf(file = sprintf("%s_%d.pdf", "plot_real_data_pres", 1), width = 12.5, height = 7.2)
  seq_values <- lapply(c(seq(1, 1*K, K)), function(x){seq(x,x+K-1)})
  par(mfrow=c(1,1), mar = c(5, 4.1, 2, 0.1), cex = 1) 
  for(i in 1:1){
    plot(y[[i]]~Xt,  cex = 1, pch = 16, col = alpha("black", 0.3), ylab = "Head acceleration (units of g)", xlab = "Time in ms")
    lines(Xt_finer, y = (out_mcycle_K20$mu_beta*ifelse(out_mcycle_K20$p > 0.50, 1, 0))%*%t(B_finer[[i]]), lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t))
    #lines(predict(smooth.spline(Xt, y[[1]]), Xt))
    lines(Xt, fitted(regression_splines), lwd = 2, col = "green")
    lines(Xt, estimativa_funcional_lasso, lwd = 2, col = "blue")
    lines(Xt, estimativa_funcionalBL, lwd = 2, col = "orange")
    #abline(h=0, col = "grey")
    lines(Xt, LL, col = "black", lty = 2)
    lines(Xt, UL, col = "black", lty = 2)
    legend("topright", legend = c("Proposed method", "Regression splines", "LASSO", "Bayesian LASSO", "Credible Band"), col = c("red", "green", "blue", "orange", "black"), lty = c(1, 1, 1, 1, 2), bty = "n")
    print(ifelse(out_mcycle_K20$p[seq_values[[i]]] > 0.50, 1, 0))
  }
}
dev.off()



# Compute R2 - smoothing splines
model <- smooth.spline(Xt, y[[1]])
yhat_smooth <- predict(model, x = Xt)$y
(sum((y[[1]]-yhat_smooth)^2))/(model$df)

RSS_s = sum((yhat_smooth - y[[1]])^2)
TSS = sum((y[[1]] - mean(y[[1]]))^2) 
R2a_s = 1 - (RSS_s/(length(y[[1]])-model$df))/(TSS/(length(y[[1]])-1))
R2a_s
# Compute R2 - regression splines
yhat_r <- fitted(regression_splines)
RSS_r <- sum((yhat_r - y[[1]])^2)
R2a_s = 1 - (RSS_r/(length(y[[1]])-20))/(TSS/(length(y[[1]])-1))
R2a_s
# Compute R2 vb
yhat_vb <- as.numeric((out_mcycle_K20$mu_beta[seq_values[[i]]]*ifelse(out_mcycle_K20$p[seq_values[[i]]] > 0.50, 1, 0))%*%t(B[[i]]))
RSS_vb = sum((yhat_vb - y[[1]])^2)
R2a_vb = 1 - (RSS_s/(length(y[[1]])-sum(ifelse(out_mcycle_K20$p[seq_values[[i]]] > 0.50, 1, 0))))/(TSS/(length(y[[1]])-1))
R2a_vb

K <- 15
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt),norder = 4, nbasis = K)
B <- lapply(1, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

initial_values <- list(p = rep(1, K),  delta2 = 90480, lambda2 = 10000, w = 10)
out_mcycle_K15 <- vb_bs_corr(y=y, B = B, m = 1, mu_ki = 0.5, lambda_1 = 10^(-10), lambda_2 = 10^(-10), delta_1 = 100, delta_2 = 99*30, maxIter = 1000, K = K, initial_values, convergence_threshold = 0.001, Xt = Xt)


seq_values <- lapply(c(seq(1, 1*K, K)), function(x){seq(x,x+K-1)})
par(mfrow=c(1,1), mar = c(5, 2, 2, 0.1), cex = 1) 
for(i in 1:1){
  plot(y[[i]]~Xt, ylab = expression(g[t]), xlab = expression(t))
  lines(Xt, y = as.numeric((out_mcycle_K15$mu_beta[seq_values[[i]]]*ifelse(out_mcycle_K15$p[seq_values[[i]]] > 0.50, 1, 0))%*%t(B[[i]])), lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t))
  print(ifelse(out_mcycle_K15$p[seq_values[[i]]] > 0.50, 1, 0))
}


K <- 10
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt), norder = 4, nbasis = K)
B <- lapply(1, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

initial_values <- list(p = rep(1, K),  delta2 = 89000, lambda2 = 10000, w = 1e3)
out_mcycle_K10 <- vb_bs_corr(y=y, B = B, m = 1, mu_ki = 0.5, lambda_1 = 1e-6, lambda_2 = 1e-6, delta_1 = 100, delta_2 = 99*30, maxIter = 1000, K = K, initial_values, convergence_threshold = 0.001, Xt = Xt)


K <- 5
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt), norder = 4, nbasis = K)
B <- lapply(1, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

initial_values <- list(p = rep(1, K),  delta2 = 91200, lambda2 = 10000, w = 1e3)
out_mcycle_K5 <- vb_bs_corr(y=y, B = B, m = 1, mu_ki = 0.5, lambda_1 = 10^(-10), lambda_2 = 10^(-10), delta_1 = 100, delta_2 = 99*30, maxIter = 1000, K = K, initial_values, convergence_threshold = 0.001, Xt = Xt)



seq_values <- lapply(c(seq(1, 1*K, K)), function(x){seq(x,x+K-1)})
par(mfrow=c(1,1), mar = c(5, 2, 2, 0.1), cex = 1) 
for(i in 1:1){
  plot(y[[i]]~Xt, ylab = expression(g[t]), xlab = expression(t))
  lines(Xt, y = as.numeric((out_mcycle_K15$mu_beta[seq_values[[i]]]*ifelse(out_mcycle_K15$p[seq_values[[i]]] > 0.50, 1, 0))%*%t(B[[i]])), lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t))
  print(ifelse(out_mcycle_K15$p[seq_values[[i]]] > 0.50, 1, 0))
}




# Second dataset - 5 curves

library(fda)
data = CanadianWeather
y = lapply(c(8,9,17,21,23), function(x){as.numeric(data$dailyAv[100:300,x,2])})
Xt <- seq(1:length(y[[1]]))
sds = lapply(1:5, function(x){sd(y[[x]])})
y_new <- lapply(1:5, function(x){y[[x]]/sds[[x]]})
cum_y <- lapply(1:5, function(x){cumsum(y_new[[x]])})
plot(y_new[[1]]~Xt, type = "l", ylim = c(0,8))
lines(y_new[[2]]~Xt, type = "l", col = "gray")
lines(y_new[[3]]~Xt, type = "l", col = "red")
lines(y_new[[4]]~Xt, type = "l", col = "blue")
lines(y_new[[5]]~Xt, type = "l", col = "green")

plot(Xt, cumsum(y_new[[1]]))

K <- 50
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt),norder = 4, nbasis = K)
fourier_basis_data <- create.fourier.basis(range(Xt), nbasis = K)
B <- lapply(1:5, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)}) #basis
B <- lapply(1:5, function(x){getbasismatrix(Xt, fourier_basis_data, nderiv = 0)[,-1]}) #fourier

# regression spline
library(splines)
regression_splines <- lm(y[[1]] ~ bs(Xt,30))
summary(regression_splines) #MSE = (0.885^2)

result <- smooth.spline(Xt, y_new[[1]])
yfdPar <- fdPar(basisBspline_Simulated_Data, Lfdobj=2, lambda = 1e-2)
yfd <- smooth.basis(Xt, y_new[[1]], yfdPar)$fd 
yhat <- eval.fd(Xt, yfd)

R <- smooth.basis(Xt, y_new[[1]], yfdPar)$penmat 
p <- smooth.basis(Xt, y_new[[1]], yfdPar)$df # cost in df
s_df <- length(Xt)-p # the higher s_df the higher lambda
sigma2hat <- (sum((Xt-yhat)^2))/(s_df)

plot(y_new[[1]] ~ Xt, col=gray(0.75)) 
lines(Xt, yhat, col="red", type="l")

K <- 10
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt),norder = 4, nbasis = K)
fourier_basis_data <- create.fourier.basis(range(Xt), nbasis = K)
B <- lapply(1:5, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)}) #basis
B <- lapply(1:5, function(x){getbasismatrix(Xt, fourier_basis_data, nderiv = 0)}) #fourier


initial_values <- list(p = rep(1, K),  delta2 = 704, lambda2 = 1000, w = 10) #500 for K=50, 473 for k= 30, 456 for K=20, 448 for K=15, 430 for K=10, 
out_weather_K20 <- vb_bs_corr(y=y_new, B = B, m = 5, mu_ki = 0.5, lambda_1 = 1e-6, lambda_2 = 1e-6, delta_1 = 100, delta_2 = 99*0.5, maxIter = 1000, K = K, initial_values, convergence_threshold = 0.001, Xt = Xt)

seq_values <- lapply(c(seq(1, 1*K, K)), function(x){seq(x,x+K-1)})
par(mfrow=c(1,1), mar = c(5, 2, 2, 0.1), cex = 1) 
for(i in 1:5){
  plot(y_new[[i]]~Xt, ylab = expression(g[t]), xlab = expression(t))
  lines(Xt, y = as.numeric((out_weather_K20$mu_beta[seq_values[[i]]]*ifelse(out_weather_K20$p[seq_values[[i]]] > 0.50, 1, 0))%*%t(B[[i]])), lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t))
  print(ifelse(out_weather_K20$p[seq_values[[i]]] > 0.50, 1, 0))
}


# Third dataset

boy.ht <- t(as.matrix(growth$hgtm))
girl.ht <- t(as.matrix(growth$hgtf))


Xt <- as.numeric(colnames(boy.ht))

plot(girl.ht[10,] ~ Xt)
lines(Xt, boy.ht[10,])

y <- lapply(1:10, function(i) as.numeric(boy.ht[i,]))
Xt <- growth$age
#Xt <- seq(0,1,len=31)
sds <- lapply(1:10, function(x){sd(y[[x]])})
y_growth <- lapply(1:10, function(x){y[[x]]/sds[[x]]})

(y[[1]]-mean(y[[1]]))%*%t(y[[1]]-mean(y[[1]]))


library(splines)
regression_splines <- lm(y_growth[[1]] ~ bs(Xt, 20, degree = 5))
summary(regression_splines) #MSE = 0.05

result <- smooth.spline(Xt, y_growth[[1]])


yfdPar <- fdPar(basisBspline_Simulated_Data, Lfdobj=4, lambda = 1e-2)
yfd <- smooth.basis(Xt, y_growth[[1]], yfdPar)$fd 
yhat <- eval.fd(Xt, yfd)

R <- smooth.basis(Xt, y_growth[[1]], yfdPar)$penmat 
p <- smooth.basis(Xt, y_growth[[1]], yfdPar)$df # cost in df
s_df <- length(Xt)-p # the higher s_df the higher lambda
sigma2hat <- (sum((Xt-yhat)^2))/(s_df)



K <- 10
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt), nbasis = K, norder = 6)
B <- lapply(1:10, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

initial_values <- list(p = rep(1, K),  delta2 = 16063.4, lambda2 = 1000, w = 2)
out_growth_K20 <- vb_bs_corr(y = y_growth, B = B, m = 10, mu_ki = 0.5, lambda_1 = 1e-6, lambda_2 = 1e-6, delta_1 = 0.1, delta_2 = 0.1, maxIter = 1000, K = K, initial_values, convergence_threshold = 0.001, Xt = Xt)

K <- 10
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt), nbasis = K, norder = 4)
B <- lapply(1:10, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

initial_values <- list(p = rep(1, K),  delta2 = 14, lambda2 = 10000, w = 10)
out_growth_K10 <- vb_bs_corr(y = y_growth, B = B, m = 10, mu_ki = 0.5, lambda_1 = 10^(-10), lambda_2 = 10^(-10), delta_1 = 100, delta_2 = 99*0.8, maxIter = 1000, K = K, initial_values, convergence_threshold = 0.001, Xt = Xt)


K <- 10
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt), nbasis = K, norder = 6)
B <- lapply(1:10, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

initial_values <- list(p = rep(1, K),  delta2 = 14, lambda2 = 10000, w = 2)
out_growth_K10 <- vb_bs_corr(y = y_growth, B = B, m = 10, mu_ki = 0.5, lambda_1 = 10^(-10), lambda_2 = 10^(-10), delta_1 = 100, delta_2 = 99*0.8, maxIter = 1000, K = K, initial_values, convergence_threshold = 0.001, Xt = Xt)



seq_values <- lapply(c(seq(1, 10*K, K)), function(x){seq(x,x+K-1)})
par(mfrow=c(1,1), mar = c(5, 2, 2, 0.1), cex = 1) 
for(i in 1:10){
  plot(y[[i]]~Xt, ylab = expression(g[t]), xlab = expression(t))
  lines(Xt, y = as.numeric((out_growth_K20$mu_beta[seq_values[[i]]]*ifelse(out_growth_K20$p[seq_values[[i]]] > 0.50, 1, 0))%*%t(B[[i]]))*sds[[i]], lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t))
  print(ifelse(out_growth_K20$p[seq_values[[i]]] > 0.50, 1, 0))
}



source("sim_data.R")
source("expected_values.R")
source("elbo_formulas.R")
source("elbo.R")
source("VB_bs_main.R")

K <- 6
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt), norder = 4, nbasis = K)
B <- list(getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0))

initial_values <- list(p = rep(1, K),  delta2 = 40000, lambda2 = 100)
out <- vb_bs(y = y, B = B, m = 5, mu_ki = 0.4, lambda_1 = 10^-10, lambda_2 = 10^-10, delta_1 = 10^-10, delta_2 = 10^-10, maxIter = 1000, K = K, initial_values = initial_values, convergence_threshold = 0.001)


mu_betas <- matrix(NA, 1, K)
p <- matrix(NA, 1, K)
for(k in 1:K){
  mu_betas[k] <- mean(out$mu_beta[grep(paste0("^beta_", k, "_[0-9]+$"),
                                       names(out$mu_beta))])
  p[k] <- mean(out$p[grep(paste0("^beta_", k, "_[0-9]+$"),
                          names(out$mu_beta))])
  
} 

z <- c()
for(k in 1:K){z[k] = ifelse(p[k] > 0.50, 1, 0)}


plot(y[[1]]~Xt, ylab = expression(g[t]), xlab = expression(t))
lines(Xt, y = as.numeric((out$mu_beta[1:K]*ifelse(out$p[1:K] > 0.50, 1, 0))%*%t(B[[1]])), lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t))

plot(y[[1]]~seq(1:length(y[[1]])), col = "gray", ylim = c(0.1, 7.7))
lines(y[[2]]~seq(1:length(y[[1]])), col = "gray")
lines(y[[3]]~seq(1:length(y[[1]])), col = "gray")
lines(y[[4]]~seq(1:length(y[[1]])), col = "gray")
lines(y[[5]]~seq(1:length(y[[1]])), col = "gray")
lines(Xt, y = as.numeric(as.vector(mu_betas*z)%*%t(B[[1]])), lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t), ylim = c(0.1, 7.7))

# Co2 data - Gaussian Processes for Machine Learning p 119


library(splines)
Xt <- c(1,2,3,5,6,7,9,10,seq(11,22,len = 12))
regression_splines <- lm(CO2_data$CO2[1:20] ~ bs(Xt,10))
summary(regression_splines) #MSE = 2.25

plot(Xt, CO2_data$CO2[1:20])

K <- 6
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt), norder = 4, nbasis = K)
B <- list(getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0))

initial_values <- list(p = rep(1, K),  delta2 = 11, lambda2 = 1000, w = 2)
out <- vb_bs_corr(y = list(CO2_data$CO2[1:20]), B = B, m = 1, mu_ki = 0.5, lambda_1 = 1e-6, lambda_2 = 1e-6, delta_1 = 100, delta_2 = 99*10, maxIter = 1000, K = K, initial_values, convergence_threshold = 0.001, Xt = Xt)

seq_values <- lapply(c(seq(1, 1*K, K)), function(x){seq(x,x+K-1)})
par(mfrow=c(1,1), mar = c(5, 2, 2, 0.1), cex = 1) 
for(i in 1:1){
  plot(list(CO2_data$CO2[1:20])[[i]]~Xt, ylab = expression(g[t]), xlab = expression(t))
  lines(Xt, y = as.numeric((out$mu_beta[seq_values[[i]]]*ifelse(out$p[seq_values[[i]]] > 0.50, 1, 0))%*%t(B[[i]])), lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t))
  print(ifelse(out$p[seq_values[[i]]] > 0.50, 1, 0))
}

# Presentation ----

plot(basisBspline_Simulated_Data, ylim = c(0,2), knots = FALSE)

{
pdf(file = sprintf("%s_%d.pdf", "real_data", 1), width = 12.5, height = 2.5)
par(mfrow=c(1,1), mar = c(2, 2, 1, 1), cex = 0.8) 
plot(basisBspline_Simulated_Data, ylim = c(0,1), knots = FALSE, col = rep(c("grey", "blue", "grey"), c(5, 7, 8)), lty = rep(c(2, 1, 2), c(5, 7, 8)))
}
dev.off()

{
pdf(file = sprintf("%s_%d.pdf", "real_data", 2), width = 12.5, height = 4.5)
par(mfrow=c(1,1), mar = c(0, 2, 1, 1), cex = 0.8) 
plot(y[[1]]~Xt,  cex = 1, pch = 16, col = alpha("black", 0.5), ylab = expression(g[t]), xlab = expression(t))
lines(Xt, y = (out_mcycle_K20$mu_beta*ifelse(out_mcycle_K20$p > 0.50, 1, 0))%*%t(B[[1]]), lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t))
}
dev.off()

{
  pdf(file = sprintf("%s_%d.pdf", "real_data", 3), width = 12.5, height = 2.5)
  par(mfrow=c(1,1), mar = c(2, 2, 1, 1), cex = 0.8) 
  plot(basisBspline_Simulated_Data, ylim = c(0,1), knots = FALSE, col = "blue", lty = 1)
}
dev.off()

{
  pdf(file = sprintf("%s_%d.pdf", "real_data", 4), width = 12.5, height = 4.5)
  par(mfrow=c(1,1), mar = c(0, 2, 1, 1), cex = 0.8) 
  plot(y[[1]]~Xt,  cex = 1, pch = 16, col = alpha("black", 0.5), ylab = expression(g[t]), xlab = expression(t))
}
dev.off()


points(2.5 + y[[1]]/50~Xt, ylab = expression(g[t]), xlab = expression(t))
lines(Xt, 2.5+ (out_mcycle_K20$mu_beta*ifelse(out_mcycle_K20$p > 0.50, 1, 0))%*%t(B[[1]])/50, col = 3, lwd = 3, lty = 1, xlab = "time in ms", ylab = "g(t)")

