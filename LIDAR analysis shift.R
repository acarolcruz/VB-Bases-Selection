# Real data application - LIDAR
# Adjusting for the intercept

source("elbo_formulas_corr.R")
source("aux_fun_corr.R")
source("VB_corr_vem.R")

# Motorcycle data
library(monomvn)
library(glmnet)
library(faraway)
library(fda)

lidar <- read.table("/Users/carol/Downloads/SemiPar/data/lidar.txt", header = TRUE)
plot(logratio~range, lidar, ylab = "logratio", xlab = "range")
y_orig <- list(lidar$logratio)
Xt <- lidar$range

intcp <- mean(y_orig[[1]][1:115]) #until before 500 as in the other analysis


y <- list(y_orig[[1]] - intcp)


# Regression splines:
K <- 50
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt), norder = 4, nbasis = K)
B <- lapply(1, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

regression_splines <- lm(y[[1]] ~ B[[1]] - 1)
summary(regression_splines) #MSE = 0.08127**2 0.0066 

K <- 30
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt), norder = 3, nbasis = K)
B <- lapply(1, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

delta2_ini <- (((length(Xt) + K + 2*100)/2) - 1)*0.01
initial_values <- list(p = rep(1, K), delta2 = delta2_ini , lambda2 = 100, w = 10)
system.time(out_lidar_K30 <- vb_bs_corr(y = y, B = B, m = 1, mu_ki = 0.5, lambda_1 = 1e-6, lambda_2 = 1e-6, delta_1 = 100, delta_2 = 99*0.01, maxIter = 1000, K = K, initial_values, convergence_threshold = 0.001, Xt = Xt))

#user  system elapsed 
#157.008   1.065 158.803 

# Compute R2 vem
seq_values <- lapply(c(seq(1, 1*K, K)), function(x){seq(x,x+K-1)})
yhat <- as.numeric((out_lidar_K30$mu_beta[seq_values[[1]]]*ifelse(out_lidar_K30$p[seq_values[[1]]] > 0.50, 1, 0))%*%t(B[[1]]))
RSS_vb <- sum((yhat - y[[1]])^2)
TSS <- sum((y[[1]] - mean(y[[1]]))^2) 
R2a_vb <- 1 - (RSS_vb/(length(y[[1]])-sum(ifelse(out_lidar_K30$p[seq_values[[1]]] > 0.50, 1, 0))))/(TSS/(length(y[[1]])-1))
R2a_vb
#0.9177

plot(y_orig[[1]]~Xt,  cex = 0.8, pch = 16, col = "grey", ylab = expression(g[t]), xlab = expression(t))
lines(Xt, y = (out_lidar_K30$mu_beta*ifelse(out_lidar_K30$p > 0.50, 1, 0))%*%t(B[[1]]) + intcp, lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t))

K <- 20
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt),norder = 4, nbasis = K)
B <- lapply(1, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

delta2_ini <- (((length(Xt) + K + 2*100)/2) - 1)*0.01
initial_values <- list(p = rep(1, K), delta2 = delta2_ini, lambda2 = 100, w = 10) 
system.time(out_lidar_K20 <- vb_bs_corr(y = y, B = B, m = 1, mu_ki = 0.5, lambda_1 = 1e-6, lambda_2 = 1e-6, delta_1 = 100, delta_2 = 99*0.01, maxIter = 1000, K = K, initial_values, convergence_threshold = 0.001, Xt = Xt))
#user  system elapsed 
#108.680   1.089 139.580 

seq_values <- lapply(c(seq(1, 1*K, K)), function(x){seq(x,x+K-1)})
yhat <- as.numeric((out_lidar_K20$mu_beta[seq_values[[1]]]*ifelse(out_lidar_K20$p[seq_values[[1]]] > 0.50, 1, 0))%*%t(B[[1]]))
RSS_vb = sum((yhat - y[[1]])^2)
TSS = sum((y[[1]] - mean(y[[1]]))^2) 
R2a_vb = 1 - (RSS_vb/(length(y[[1]])-sum(ifelse(out_lidar_K20$p[seq_values[[1]]] > 0.50, 1, 0))))/(TSS/(length(y[[1]])-1))
R2a_vb
#0.9220

plot(y_orig[[1]]~Xt,  cex = 0.8, pch = 16, col = "grey", ylab = expression(g[t]), xlab = expression(t))
lines(Xt, y = (out_lidar_K20$mu_beta*ifelse(out_lidar_K20$p > 0.50, 1, 0))%*%t(B[[1]]) + intcp, lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t))

K <- 15
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt),norder = 4, nbasis = K)
B <- lapply(1, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

delta2_ini <- (((length(Xt) + K + 2*100)/2) - 1)*0.01
initial_values <- list(p = rep(1, K), delta2 = delta2_ini, lambda2 = 100, w = 10)
system.time(out_lidar_K15 <- vb_bs_corr(y = y, B = B, m = 1, mu_ki = 0.5, lambda_1 = 1e-6, lambda_2 = 1e-6, delta_1 = 100, delta_2 = 99*0.01, maxIter = 1000, K = K, initial_values, convergence_threshold = 0.001, Xt = Xt))

#user  system elapsed 
# 97.718   0.677  98.964  

seq_values <- lapply(c(seq(1, 1*K, K)), function(x){seq(x,x+K-1)})
yhat <- as.numeric((out_lidar_K15$mu_beta[seq_values[[1]]]*ifelse(out_lidar_K15$p[seq_values[[1]]] > 0.50, 1, 0))%*%t(B[[1]]))
RSS_vb <- sum((yhat - y[[1]])^2)
TSS <- sum((y[[1]] - mean(y[[1]]))^2) 
R2a_vb <- 1 - (RSS_vb/(length(y[[1]])-sum(ifelse(out_lidar_K15$p[seq_values[[1]]] > 0.50, 1, 0))))/(TSS/(length(y[[1]])-1))
R2a_vb
#0.9230

plot(y_orig[[1]]~Xt,  cex = 0.8, pch = 16, col = "grey", ylab = expression(g[t]), xlab = expression(t))
lines(Xt, y = (out_lidar_K15$mu_beta*ifelse(out_lidar_K15$p > 0.50, 1, 0))%*%t(B[[1]]) + intcp, lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t))

K <- 10
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt),norder = 4, nbasis = K)
B <- lapply(1, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

delta2_ini <- (((length(Xt) + K + 2*100)/2) - 1)*0.01
initial_values <- list(p = rep(1, K), delta2 = delta2_ini, lambda2 = 100, w = 10)
system.time(out_lidar_K10 <- vb_bs_corr(y = y, B = B, m = 1, mu_ki = 0.5, lambda_1 = 1e-6, lambda_2 = 1e-6, delta_1 = 100, delta_2 = 99*0.01, maxIter = 1000, K = K, initial_values, convergence_threshold = 0.001, Xt = Xt))
#user  system elapsed 
#102.719   0.579 103.457 

seq_values <- lapply(c(seq(1, 1*K, K)), function(x){seq(x,x+K-1)})
yhat <- as.numeric((out_lidar_K10$mu_beta[seq_values[[1]]]*ifelse(out_lidar_K10$p[seq_values[[1]]] > 0.50, 1, 0))%*%t(B[[1]]))
RSS_vb <- sum((yhat- y[[1]])^2)
TSS <- sum((y[[1]] - mean(y[[1]]))^2) 
R2a_vb <- 1 - (RSS_vb/(length(y[[1]])-sum(ifelse(out_lidar_K10$p[seq_values[[1]]] > 0.50, 1, 0))))/(TSS/(length(y[[1]])-1))
R2a_vb
#0.9172


plot(y_orig[[1]]~Xt,  cex = 0.8, pch = 16, col = "grey", ylab = expression(g[t]), xlab = expression(t))
lines(Xt, y = ((out_lidar_K10$mu_beta*ifelse(out_lidar_K10$p > 0.50, 1, 0))%*%t(B[[1]])) + intcp, lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t))

K <- 6
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt),norder = 4, nbasis = K)
B <- lapply(1, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

delta2_ini <- (((length(Xt) + K + 2*100)/2) - 1)*0.01
initial_values <- list(p = rep(1, K), delta2 = delta2_ini, lambda2 = 100, w = 10)
system.time(out_lidar_K6 <- vb_bs_corr(y = y, B = B, m = 1, mu_ki = 0.5, lambda_1 = 1e-6, lambda_2 = 1e-6, delta_1 = 100, delta_2 = 99*0.01, maxIter = 1000, K = K, initial_values, convergence_threshold = 0.0001, Xt = Xt))
#user  system elapsed 
#89.462   0.480  90.093 

seq_values <- lapply(c(seq(1, 1*K, K)), function(x){seq(x,x+K-1)})
yhat <- as.numeric((out_lidar_K6$mu_beta[seq_values[[1]]]*ifelse(out_lidar_K6$p[seq_values[[1]]] > 0.50, 1, 0))%*%t(B[[1]]))
RSS_vb <- sum((yhat- y[[1]])^2)
TSS <- sum((y[[1]] - mean(y[[1]]))^2) 
R2a_vb <- 1 - (RSS_vb/(length(y[[1]])-sum(ifelse(out_lidar_K6$p[seq_values[[1]]] > 0.50, 1, 0))))/(TSS/(length(y[[1]])-1))
R2a_vb
#0.9003

# gcv ----

gcv <- function(i, out, y, B, seq_values){
  p_i <- out$p[seq_values[[i]]]
  Zi_hat <- ifelse(out$p[seq_values[[i]]] > 0.50, 1, 0)
  
  EG_i <- p_i*t(B[[i]])
  
  S_k <- E_inv_sigma2(out$delta1, out$delta2)*B[[i]]%*%diag(Zi_hat)%*%out$Sigma_beta[,,i]%*%EG_i%*%solve(calPsi(Xt, Xt, w = out$w))
  
  S_k <- E_inv_sigma2(out$delta1, out$delta2)*B[[i]]%*%diag(Zi_hat)%*%solve(E_inv_sigma2(out$delta1, out$delta2)*(diag(E_inv_tau2(out$lambda1, out$lambda2), K) + EG_i%*%solve(calPsi(Xt, Xt, w = out$w))%*%t(EG_i)))%*%EG_i%*%solve(calPsi(Xt, Xt, out$w))
  
  yhat <- as.numeric((out$mu_beta[seq_values[[i]]]*ifelse(out$p[seq_values[[i]]] > 0.50, 1, 0))%*%t(B[[i]]))
  RSS_vb = sum((yhat- y[[i]])^2)
  
  numerador <- t(y[[i]] - B[[i]]%*%(Zi_hat*out$mu_beta[seq_values[[i]]]))%*%(y[[i]] - B[[i]]%*%diag(Zi_hat)%*%out$mu_beta[seq_values[[i]]])
  
  GCV_k <- (1/length(Xt))*numerador/((1 - (1/length(Xt))*sum(diag(S_k)))^2) 
  return(GCV_k)
}  

K <- 6
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt),norder = 4, nbasis = K)
B <- lapply(1, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})
seq_values <- lapply(c(seq(1, 1*K, K)), function(x){seq(x,x+K-1)})

gcv_lidar6 <- gcv(i=1,out=out_lidar_K6,y=y,B=B, seq_values = seq_values)

K <- 10
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt),norder = 4, nbasis = K)
B <- lapply(1, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})
seq_values <- lapply(c(seq(1, 1*K, K)), function(x){seq(x,x+K-1)})

gcv_lidar10 <- gcv(i=1,out=out_lidar_K10,y=y,B=B, seq_values = seq_values)


K <- 15
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt),norder = 4, nbasis = K)
B <- lapply(1, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})
seq_values <- lapply(c(seq(1, 1*K, K)), function(x){seq(x,x+K-1)})

gcv_lidar15 <- gcv(i=1,out=out_lidar_K15,y=y,B=B, seq_values = seq_values)

K <- 20
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt),norder = 4, nbasis = K)
B <- lapply(1, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})
seq_values <- lapply(c(seq(1, 1*K, K)), function(x){seq(x,x+K-1)})

gcv_lidar20 <- gcv(i=1,out=out_lidar_K20,y=y,B=B, seq_values = seq_values)

K <- 30
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt),norder = 4, nbasis = K)
B <- lapply(1, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})
seq_values <- lapply(c(seq(1, 1*K, K)), function(x){seq(x,x+K-1)})

gcv_lidar30 <- gcv(i=1,out=out_lidar_K30,y=y,B=B, seq_values = seq_values)

# optimal K = 10 (smaller gcv)
gcv_lidar6
gcv_lidar10 
gcv_lidar15 
gcv_lidar20 
gcv_lidar30 

GCV_v <- c(gcv_lidar6,
           gcv_lidar10, 
           gcv_lidar15,
           gcv_lidar20, 
           gcv_lidar30)


gcv_res <- cbind(K = c(6,10,15,20,30), gcv = GCV_v)
which.min(gcv_res[,2])


# Comparative methods
K <- 15
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt), norder = 4, nbasis = K)
B <- lapply(1, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

#smoothing splines:
loglam  = seq(1,9,0.25)
nlam    = length(loglam)
dfsave  = rep(NA,nlam)
gcvsave = rep(NA,nlam)

for (ilam in 1:nlam) {
  cat(paste('log10 lambda = ',loglam[ilam], '\n '))
  lambda   = 10^loglam[ilam]
  fdParobj = fdPar(basisBspline_Simulated_Data, 2, lambda)
  smoothlist = smooth.basis(Xt, y[[1]],
                            fdParobj)
  dfsave[ilam]  = smoothlist$df
  gcvsave[ilam] = sum(smoothlist$gcv)
}

which.min(gcvsave) # log10(lambda) = 1.25 => lambda = 10^1.25 

lambda   = 10^loglam[which.min(gcvsave)]
fdParobj = fdPar(basisBspline_Simulated_Data, 2, lambda)
ss.fit = smooth.basis(Xt, y[[1]], fdParobj)
ss.fd = ss.fit$fd

y_hat_ss <- eval.fd(Xt, ss.fd)
RSS_r <- sum((y_hat_ss - y[[1]])^2)
TSS = sum((y[[1]] - mean(y[[1]]))^2) 
R2a_ss = 1 - (RSS_r/(length(y[[1]])-ss.fit$df))/(TSS/(length(y[[1]])-1))
R2a_ss
#0.9210
# Regression splines:

regression_splines <- lm(y[[1]] ~ B[[1]] - 1)

RSS <- sum((y[[1]]-fitted(regression_splines))^2)
TSS <- sum((y[[1]]- mean(y[[1]]))^2)
R2a_rs <- 1-(RSS/(length(y[[1]])-K))/(TSS/(length(y[[1]])-1))
#0.9205

# Bayesian lasso ----

#system.time(BL <- blasso(X = B[[1]], y = y[[1]], T = 10000, thin = 1, 
#                         beta=rep(-1,1*20),lambda2=1,s2=1,
#                         RJ = F, rao.s2=F, icept = T, normalize=F))
# 2 chains
df.blasso <- NULL

BL1 <- blasso(X = B[[1]], y = y[[1]], T = 10000, thin = 1, 
              beta=rep(-1, K),lambda2 = 1, s2 = 1,
              RJ = F, rao.s2=F, icept = T, normalize=F)
BL2 <- blasso(X = B[[1]], y = y[[1]], T = 10000, thin = 1, 
              beta=rep(1, K),lambda2 = 5, s2 = 5,
              RJ = F, rao.s2 = F, icept = T, normalize = F)

burnIn=5000
df.blasso1=as.data.frame(cbind(lambda=BL1$lambda2,s2=BL1$s2,beta0 = BL1$mu, beta=BL1$beta))
df.blasso1=df.blasso1[(burnIn+1):nrow(df.blasso1),]
df.blasso2=as.data.frame(cbind(lambda=BL2$lambda2,s2=BL2$s2,beta0 = BL1$mu, beta=BL2$beta))
df.blasso2=df.blasso2[(burnIn+1):nrow(df.blasso2),]

thin=50
df.blasso1 = df.blasso1[rep(c(T,rep(F,thin-1)),length.out=dim(df.blasso1)[1]),]
df.blasso2 = df.blasso2[rep(c(T,rep(F,thin-1)),length.out=dim(df.blasso2)[1]),]
df.blasso = rbind(df.blasso1,df.blasso2)


beta_BL=apply(df.blasso, 2, mean)
estimativa_funcionalBL <- cbind(rep(1,length(Xt)),B[[1]])%*%as.matrix(beta_BL[-c(1,2)])

RSS <- sum((y[[1]]-estimativa_funcionalBL)^2)
TSS <- sum((y[[1]]- mean(y[[1]]))^2)
K_model <- sum(abs(beta_BL[-c(1,2)])>0)
R2_adj_BL <- 1-(RSS/(length(y[[1]])-K_model))/(TSS/(length(y[[1]])-1))
# 0.9188

# LASSO ----
regularizador <- seq(0.001, 1, length = 9)

result_LASSO <- NULL
R2_adj_LASSO_mat <- c()
time_lasso <- c()
for(lambda in 1:length(regularizador)){
  #setTxtProgressBar(pb,lambda)
  
  lambdaLIVROtexto <- regularizador[lambda] #parametrização de livro texto.
  time_lasso <- c(time_lasso, system.time(LASSO <- glmnet(B[[1]], y[[1]], alpha = 1,lambda = lambdaLIVROtexto/(2*nrow(B[[1]])), standardize = T, intercept = F)))
  result_LASSO[[lambda]] <- LASSO
  betaLASSO <- as.numeric(coef(LASSO)[-1])
  estimativa_funcional_lasso <- B[[1]]%*%betaLASSO
  
  RSS=sum((y[[1]] - estimativa_funcional_lasso)^2)
  TSS=sum((y[[1]] - mean(y[[1]]))^2)
  K_model <- sum(abs(betaLASSO)>0)
  R2_adj_LASSO_mat[lambda] <- 1-mean((RSS/(length(y[[1]])-K_model))/(TSS/(length(y[[1]])-1)))
}
mean(R2_adj_LASSO_mat)
lambda <- which.max(R2_adj_LASSO_mat)
# 0.9213


# to plot LASSO's results
LASSO <- glmnet(B[[1]], y[[1]], alpha = 1, lambda = regularizador[which.max(R2_adj_LASSO_mat)]/(2*nrow(B[[1]])), standardize = F, intercept = F)
betaLASSO <- as.numeric(coef(LASSO)[-1])
estimativa_funcional_lasso <- B[[1]]%*%betaLASSO

# time 
mean(time_lasso[seq(from = 3, to = 45, 5)])
# 0.0023

# Credible band as in Pedro's article ----
LL <- NULL
UL <- NULL
estimates <- NULL
out <- out_lidar_K15
betas_sample <- MASS::mvrnorm(200, mu = out$mu_beta,  Sigma = out$Sigma_beta[,,1])
z_sample <- matrix(rbinom(15*200, 1, prob = out$p), ncol = 15, byrow = TRUE)

estimates <- cbind(estimates, (betas_sample*z_sample))

curve_f <- list()
for(s in 1:200){
  curve_f[[s]] <- estimates[s,]%*%t(B[[1]])
}  

LL <- apply(do.call(rbind,curve_f),2,function(i)quantile(i,probs = c(0.025,0.975)))[1,]
UL <-  apply(do.call(rbind,curve_f),2,function(i)quantile(i,probs = c(0.025,0.975)))[2,]

R2_res <- data.frame(cbind(method = c('VEM', 'RS', 'SS', 'Lasso', 'BLasso'), 
                           R2 = round(c(R2a_vb, R2a_rs, R2a_ss, mean(R2_adj_LASSO_mat), R2_adj_BL), 4))) 
save(y_orig, y, intcp, B,
     R2_res,
     gcv_res, 
     out_lidar_K15,
     regression_splines,
     ss.fd,
     beta_BL,
     betaLASSO, file = 'LIDAR_results.RData')

# comparing VEM with smoothing splines and regression splines
{
  pdf(file = sprintf("%s_%d.pdf", "plot_lidar_new", 1), width = 12.5, height = 7.2)
  seq_values <- lapply(c(seq(1, 1*K, K)), function(x){seq(x,x+K-1)})
  par(mfrow=c(1,1), mar = c(5, 4.1, 2, 0.1), cex = 1) 
  for(i in 1:1){
    plot(y_orig[[i]]~Xt,  cex = 0.8, pch = 16, col = "grey", ylab = expression(g[t]), xlab = expression(t))
    #(out_lidar_K20$mu_beta*ifelse(out_lidar_K20$p > 0.50, 1, 0))%*%t(B[[i]])
    lines(Xt, y = (out$mu_beta*out$p)%*%t(B[[i]]) + intcp, lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t))
    lines(Xt, fitted(regression_splines) + intcp, lwd = 2, col = "blue", lty = 2)
    lines(Xt, predict(smooth.spline(Xt, y[[1]]), Xt)$y + intcp, lwd = 2, col = "orange",  lty = 3)
    #lines(Xt, LL, col = "black", lty = 2)
    #lines(Xt, UL, col = "black", lty = 2)
    legend("topright", legend = c("Proposed method", "Regression splines", "Smoothing splines"), col = c("red", "blue", "orange"), lty = c(1, 1, 1), bty = "n")
    print(ifelse(out$p[seq_values[[i]]] > 0.50, 1, 0))
  }
}
dev.off()


library(scales)
{
  pdf(file = sprintf("%s_%d.pdf", "plot_real_data_lidar_new", 1), width = 12.5, height = 7.2)
  seq_values <- lapply(c(seq(1, 1*K, K)), function(x){seq(x,x+K-1)})
  par(mfrow=c(1,1), mar = c(5, 4.1, 2, 0.1), cex = 1.2) 
  for(i in 1:1){
    plot(y_orig[[i]]~Xt,  cex = 0.8, pch = 16, col = 'grey', ylab = "logratio", xlab = "range")
    lines(Xt, y = (out$mu_beta*ifelse(out$p > 0.50, 1, 0))%*%t(B[[i]]) + intcp, lwd = 2, col= "red", ylab = expression(g[t]), xlab = expression(t))
    lines(Xt, B[[1]]%*%ss.fd$coefs + intcp, lwd = 2, col = "purple")
    lines(Xt, B[[1]]%*%regression_splines$coefficients + intcp, lwd = 2, lty = 4, col = "green")
    lines(Xt, estimativa_funcional_lasso + intcp, lwd = 2, col = "blue",lty = 2 )
    lines(Xt, estimativa_funcionalBL + intcp, lwd = 3, col = "orange", lty = 3)
    #abline(h=0, col = "grey")
    #polygon(c(Xt, rev(Xt)), c(LL, rev(UL)), col = "#00000020", border = NA)
    #lines(Xt, LL, col = "black", lty = 2)
    #lines(Xt, UL, col = "black", lty = 2)
    legend("topright", legend = c("Proposed method", "LASSO", 'Regression splines', 'Bayesian LASSO', 'Smoothing splines'), col = c("red", 'blue','green', 'orange', 'purple'), lwd = c(2, 2, 2, 2, 2), lty = c(1, 1, 1, 1, 1), bty = "n")
    print(ifelse(out$p[seq_values[[i]]] > 0.50, 1, 0))
  }
}
dev.off()

# Presentation ----

which(ifelse(out_lidar_K10_old$p > 0.50, 1, 0) == 1)
colors = rep("blue", K)
colors[which(ifelse(out_lidar_K10_old$p > 0.50, 1, 0) == 0)] <- "black"
{
  pdf(file = sprintf("%s_%d.pdf", "lidar_basis_vem", 1), width = 12.5, height = 2.5)
  par(mfrow=c(1,1), mar = c(5, 4.1, 2, 0.1), cex = 0.8) 
  plot(basisBspline_Simulated_Data, ylim = c(0,1), knots = FALSE, col = colors,
       lty = rep(c(2, 1), c(5,5)), lwd = rep(c(1.5, 2), c(5,5)))
}
dev.off()

# testing putting an intercept

K <- 15
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt),norder = 4, nbasis = K)
B <- lapply(1, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

B_inter <- B
B_inter[[1]] <- cbind(c(rep(1, 107), rep(0, length(Xt) - 107)), B[[1]])

K <- ncol(B_inter[[1]])

list(y - mean(y[[1]]))

delta2_ini <- (((length(Xt) + K + 2*100)/2) - 1)*0.01
initial_values <- list(p = rep(1, K), delta2 = delta2_ini, lambda2 = 100, w = 10) 
system.time(out_lidar_K15_inter <- vb_bs_corr(y =y, B = B_inter, m = 1, mu_ki = 0.5, lambda_1 = 1e-6, lambda_2 = 1e-6, delta_1 = 100, delta_2 = 99*0.01, maxIter = 1000, K = K, initial_values, convergence_threshold = 0.001, Xt = Xt))

seq_values <- lapply(c(seq(1, 1*K, K)), function(x){seq(x,x+K-1)})
yhat <- as.numeric((out_lidar_K15_inter$mu_beta[seq_values[[1]]]*ifelse(out_lidar_K15_inter$p[seq_values[[1]]] > 0.50, 1, 0))%*%t(B_inter[[1]]))
RSS_vb = sum((yhat - y[[1]])^2)
TSS = sum((y[[1]] - mean(y[[1]]))^2) 
R2a_vb = 1 - (RSS_vb/(length(y[[1]])-sum(ifelse(out_lidar_K15_inter$p[seq_values[[1]]] > 0.50, 1, 0))))/(TSS/(length(y[[1]])-1))
R2a_vb
#0.9059

plot(y[[i]]~Xt,  cex = 0.8, pch = 16, col = "grey", ylab = expression(g[t]), xlab = expression(t))
lines(Xt, y = (out_lidar_K15_inter$mu_beta*ifelse(out_lidar_K15_inter$p > 0.50, 1, 0))%*%t(B_inter[[i]]), lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t))



K <- 15
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt),norder = 4, nbasis = K)
B <- lapply(1, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

B_inter <- B
B_inter[[1]] <- cbind(c(rep(1, 107), rep(0, length(Xt) - 107)), B[[1]])

K <- ncol(B_inter[[1]])

y_test <- y[[1]] + 2

delta2_ini <- (((length(Xt) + K + 2*100)/2) - 1)*0.01
initial_values <- list(p = c(rep(0, 7), rep(1, 8)), delta2 = delta2_ini, lambda2 = 100, w = 10) 
system.time(out_lidar_K15<- vb_bs_corr(y = list(y_test), B = B, m = 1, mu_ki = 0.5, lambda_1 = 1e-6, lambda_2 = 1e-6, delta_1 = 100, delta_2 = 99*0.01, maxIter = 1000, K = K, initial_values, convergence_threshold = 0.001, Xt = Xt))

seq_values <- lapply(c(seq(1, 1*K, K)), function(x){seq(x,x+K-1)})
plot(y_test~Xt,  cex = 0.8, pch = 16, col = "grey", ylab = expression(g[t]), xlab = expression(t))
lines(Xt, y = (out_lidar_K15$mu_beta*ifelse(out_lidar_K15$p > 0.50, 1, 0))%*%t(B[[i]]), lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t))

plot(basisBspline_Simulated_Data, ylim = c(0,1), knots = FALSE)



plot(Xt[1:115], y_orig[[1]][1:115])

y_flat <- y_orig[[1]][1:115]

mean(y_flat)
Xt_flat <- Xt[1:115]

K <- 10
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt),norder = 4, nbasis = K)
B <- lapply(1, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})


delta2_ini <- (((length(Xt) + K + 2*100)/2) - 1)*0.01
initial_values <- list(p = rep(1, K), delta2 = delta2_ini, lambda2 = 100, w = 10) 
system.time(out_lidar_K15<- vb_bs_corr(y = list(y_orig[[1]] + abs(mean(y_flat))), B = B, m = 1, mu_ki = 0.5, lambda_1 = 1e-6, lambda_2 = 1e-6, delta_1 = 100, delta_2 = 99*0.01, maxIter = 1000, K = K, initial_values, convergence_threshold = 0.001, Xt = Xt))

seq_values <- lapply(c(seq(1, 1*K, K)), function(x){seq(x,x+K-1)})
plot(y_orig[[1]] ~ Xt,  cex = 0.8, pch = 16, col = "grey", ylab = expression(g[t]), xlab = expression(t))
lines(Xt, y = (out_lidar_K15$mu_beta*ifelse(out_lidar_K15$p > 0.50, 1, 0))%*%t(B[[1]]) - abs(mean(y_flat)), lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t))

y2 <-c(y_flat - mean(y_flat), y[[1]][116:221])
plot(Xt, y2)

# Mixed of ordem 2 and ordem 4 B splines

K <- 10
basisBspline_Simulated_Data1 <- create.bspline.basis(range(Xt),norder = 4, nbasis = K)
B1 <- lapply(1, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data1, nderiv = 0)})

K <- 20
basisBspline_Simulated_Data2 <- create.bspline.basis(range(Xt), norder = 2, nbasis = K)
B2 <- lapply(1, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data2, nderiv = 0)})

B = list(cbind(B2[[1]], B1[[1]]))
K <- ncol(B[[1]])

delta2_ini <- (((length(Xt) + K + 2*100)/2) - 1)*0.01
initial_values <- list(p = rep(1, K), delta2 = delta2_ini, lambda2 = 100, w = 10) 
system.time(out_lidar_K20_mixed<- vb_bs_corr(y = list(y[[1]] - mean(y[[1]])), B = B, m = 1, mu_ki = 0.5, lambda_1 = 1e-6, lambda_2 = 1e-6, delta_1 = 100, delta_2 = 99*0.01, maxIter = 1000, K = K, initial_values, convergence_threshold = 0.001, Xt = Xt))

seq_values <- lapply(c(seq(1, 1*K, K)), function(x){seq(x,x+K-1)})
plot(y[[1]]~Xt,  cex = 0.8, pch = 16, col = "grey", ylab = expression(g[t]), xlab = expression(t))
lines(Xt, y = (out_lidar_K20_mixed$mu_beta*ifelse(out_lidar_K20_mixed$p > 0.50, 1, 0))%*%t(B[[1]]) + mean(y[[1]]), lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t))




plot(basisBspline_Simulated_Data1, ylim = c(0,1), knots = FALSE)
lines(Xt, basisBspline_Simulated_Data2, ylim = c(0,1), knots = FALSE)

plot(Xt, B[[1]][,1], type = "l")
for(i in 2:30){lines(Xt, B[[1]][,i], type = "l")}
