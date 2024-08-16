# Weather data
# Real data application
source("elbo_formulas_corr.R")
source("aux_fun_corr.R")
source("VB_corr_vem.R")


# Second dataset - 6 curves

library(fda)
library(invgamma)
data <- CanadianWeather
cities <- c('Montreal','Quebec', 'Arvida', 'Bagottville', 'Sherbrooke', "Vancouver")
y = lapply(1:6, function(x){as.numeric(data$dailyAv[1:365,cities[x],1])}) #temperature
Xt <- 1:length(y[[1]])
sds = lapply(1:6, function(x){sd(y[[x]])})
y_new <- lapply(1:6, function(x){y[[x]]/sds[[x]]})

# Raw curves ----
{
  pdf(file = sprintf("%s_%d.pdf", "plot_weather_raw", 1), width = 12.5, height = 7.2)
  plot(y_new[[1]]*sds[[1]]~Xt, type = "l", col = "orange", xlab = "Day", ylab = "Temperature", ylim = c(-19.3, 22.5))
  lines(y_new[[2]]*sds[[2]]~Xt, type = "l", col = "purple")
  lines(y_new[[3]]*sds[[3]]~Xt, type = "l", col = "green")
  lines(y_new[[4]]*sds[[4]]~Xt, type = "l", col = "blue")
  lines(y_new[[5]]*sds[[5]]~Xt, type = "l", col = "red")
  lines(y_new[[6]]*sds[[6]]~Xt, type = "l", col = "black")
  legend("topright", legend = cities, col = c("orange", "purple", "green","blue", "red", "black"), lty = c(1, 1, 1, 1, 1, 1), bty = "n")
  dev.off()
}

apply(matrix(unlist(lapply(1:6, function(i){c(min(y[[i]]), max(y[[i]]))})),ncol = 2, byrow = T), 2, max)
apply(matrix(unlist(lapply(1:6, function(i){c(min(y[[i]]), max(y[[i]]))})),ncol = 2, byrow = T), 2, min)

# regression spline (to get estimated mse)
library(splines) 
regression_splines <- lm(y_new[[1]] ~ bs(Xt,50))
summary(regression_splines) #MSE = (0.0566^2)

# Temperature - whole year ---- 
K <- 50
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt),norder = 4, nbasis = K)
B <- lapply(1:6, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)}) 

p_initial <- rep(1, K)
initial_values <- list(p = rep(p_initial, 1), delta2 = 4, lambda2 = 1000, w = 100)
time_tempK50_year <- system.time(out_tempK50_year <- vb_bs_corr(y = y_new, B = B, m = 6, mu_ki = 0.5, lambda_1 = 0.5, lambda_2 = 0.5, delta_1 = 10, delta_2 = 9*0.01, maxIter = 1000, K = K, initial_values = initial_values, convergence_threshold = 0.001, Xt = Xt))

K <- 30
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt),norder = 4, nbasis = K)
B <- lapply(1:6, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

p_initial <- rep(1, K) 
initial_values <- list(p = rep(p_initial, 2), delta2 = 4, lambda2 = 1000, w = 100)
time_tempK30_year <- system.time(out_tempK30_year <- vb_bs_corr(y = y_new, B = B, m = 6, mu_ki = 0.5, lambda_1 = 0.5, lambda_2 = 0.5, delta_1 = 10, delta_2 = 9*0.01, maxIter = 1000, K = K, initial_values = initial_values, convergence_threshold = 0.001, Xt = Xt))

K <- 20
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt),norder = 4, nbasis = K)
B <- lapply(1:6, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

p_initial <- rep(1, K) 
initial_values <- list(p = rep(p_initial, 2), delta2 = 3.5, lambda2 = 1000, w = 100)
time_tempK20_year <- system.time(out_tempK20_year <- vb_bs_corr(y = y_new, B = B, m = 6, mu_ki = 0.5, lambda_1 = 0.5, lambda_2 = 0.5, delta_1 = 10, delta_2 = 9*0.01, maxIter = 1000, K = K, initial_values = initial_values, convergence_threshold = 0.001, Xt = Xt))

K <- 15
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt),norder = 4, nbasis = K)
B <- lapply(1:6, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

p_initial <- rep(1, K)
initial_values <- list(p = rep(p_initial, 1), delta2 = 3.4, lambda2 = 1000, w = 100)
time_tempK15_year <- system.time(out_tempK15_year <- vb_bs_corr(y = y_new, B = B, m = 6, mu_ki = 0.5, lambda_1 = 0.5, lambda_2 = 0.5, delta_1 = 10, delta_2 = 9*0.01, maxIter = 1000, K = K, initial_values = initial_values, convergence_threshold = 0.001, Xt = Xt))

K <- 10
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt),norder = 4, nbasis = K)
B <- lapply(1:6, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

p_initial <- rep(1, K)
initial_values <- list(p = rep(p_initial, 1), delta2 = 3.39, lambda2 = 1000, w = 100)
time_tempK10_year <- system.time(out_tempK10_year <- vb_bs_corr(y = y_new, B = B, m = 6, mu_ki = 0.5, lambda_1 = 0.5, lambda_2 = 0.5, delta_1 = 10, delta_2 = 9*0.01, maxIter = 1000, K = K, initial_values = initial_values, convergence_threshold = 0.001, Xt = Xt))

# Some hyperparamters considered for sigma2
#10,5*0.05
#40,2
#10,9*0.005
#2 and 1*0.005
#5, 2

# GCV evaluation

K <- 10
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt),norder = 4, nbasis = K)
B <- lapply(1:6, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

seq_values <- lapply(c(seq(1, 6*K, K)), function(x){seq(x,x+K-1)})

out <- out_tempK10_year
gcv_temp10 <- sapply(1:6, function(i){gcv(i=i,out=out,y=y_new,B=B, seq_values = seq_values)})

prop_temp10 <- round(sapply(1:6, function(i){sum(ifelse(out$p[seq_values[[i]]] > 0.50, 1, 0))/K}),4)

K <- 15
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt),norder = 4, nbasis = K)
B <- lapply(1:6, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

seq_values <- lapply(c(seq(1, 6*K, K)), function(x){seq(x,x+K-1)})

out <- out_tempK15_year
gcv_temp15 <- sapply(1:6, function(i){gcv(i=i,out=out,y=y_new,B=B, seq_values = seq_values)})

prop_temp15 <- round(sapply(1:6, function(i){sum(ifelse(out$p[seq_values[[i]]] > 0.50, 1, 0))/K}),4)

K <- 20
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt),norder = 4, nbasis = K)
B <- lapply(1:6, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

seq_values <- lapply(c(seq(1, 6*K, K)), function(x){seq(x,x+K-1)})

out <- out_tempK20_year
gcv_temp20 <- sapply(1:6, function(i){gcv(i=i,out=out,y=y_new,B=B, seq_values = seq_values)})

prop_temp20 <- round(sapply(1:6, function(i){sum(ifelse(out$p[seq_values[[i]]] > 0.50, 1, 0))/K}),4)

K <- 30
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt),norder = 4, nbasis = K)
B <- lapply(1:6, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

seq_values <- lapply(c(seq(1, 6*K, K)), function(x){seq(x,x+K-1)})

out <- out_tempK30_year
gcv_temp30 <- sapply(1:6, function(i){gcv(i=i,out=out,y=y_new,B=B, seq_values = seq_values)})

prop_temp30 <- round(sapply(1:6, function(i){sum(ifelse(out$p[seq_values[[i]]] > 0.50, 1, 0))/K}),4)

K <- 50
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt),norder = 4, nbasis = K)
B <- lapply(1:6, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

seq_values <- lapply(c(seq(1, 6*K, K)), function(x){seq(x,x+K-1)})

out <- out_tempK50_year
gcv_temp50 <- sapply(1:6, function(i){gcv(i=i,out=out,y=y_new,B=B, seq_values = seq_values)})

prop_temp50 <- round(sapply(1:6, function(i){sum(ifelse(out$p[seq_values[[i]]] > 0.50, 1, 0))/K}),4)

optimal_K <- apply(data.frame(gcv_temp10, gcv_temp15, gcv_temp20, gcv_temp30, gcv_temp50), 1, function(x){c(10,15,20,30,50)[which.min(x)]})

data.frame(curve = c(1,2,3,4,5,6), K_10 = gcv_temp10, prop_10 = prop_temp10, K_15 = gcv_temp15, prop_15 = prop_temp15, K_20 = gcv_temp20, prop_20 = prop_temp20, K_30 = gcv_temp30, prop_30 = prop_temp30, K_50 = gcv_temp50, prop_50 = prop_temp50, optimal_K_gcv = optimal_K)

# Comparative analysis

# Comparison with other methods (using optimal K)----

# Compute R2 ---- 

K <- 30
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt), norder = 4, nbasis = K)
B_all <- lapply(1:6, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

K <- 20
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt), norder = 4, nbasis = K)
B_all[[1]] <- getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)
B_all[[2]] <- getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)

K_all <- optimal_K

# regression splines
R2a_rs <- sapply(1:6, function(i){
  regression_splines <- lm(y_new[[i]] ~ B_all[[i]] - 1)
  yhat_r <- fitted(regression_splines)
  RSS_r <- sum((yhat_r - y_new[[i]])^2)
  TSS = sum((y_new[[1]] - mean(y_new[[1]]))^2) 
  R2a_rs = 1 - (RSS_r/(length(y_new[[i]])-K))/(TSS/(length(y_new[[i]])-1))
  R2a_rs
})

# LASSO 
estimativa_funcional_LASSO <- NULL
regularizador <- seq(0.001, 1, length = 9)
R2_lasso <- c()
for(i in 1:6){
  result_LASSO <- NULL
  R2_adj_LASSO_mat <- c()
  time_lasso <- c()
  for(lamb in 1:length(regularizador)){
    lambdaLIVROtexto <- regularizador[lamb] #parametrização de livro texto.
    time_lasso <- c(time_lasso, system.time(LASSO <- glmnet(B_all[[i]], y_new[[i]], alpha = 1,lambda = lambdaLIVROtexto/(2*nrow(B_all[[i]])), standardize = F, intercept = F)))
    result_LASSO[[lamb]] <- LASSO
    betaLASSO <- as.numeric(coef(LASSO)[-1])
    estimativa_funcional_lasso <- B_all[[i]]%*%betaLASSO
    
    RSS=sum((y_new[[i]] - estimativa_funcional_lasso)^2)
    TSS=sum((y_new[[i]] - mean(y_new[[i]]))^2)
    K_model <- sum(abs(betaLASSO)>0)
    R2_adj_LASSO_mat[lamb] <- 1-mean((RSS/(length(y_new[[i]])-K_model))/(TSS/(length(y_new[[i]])-1)))
  }
  lamb_op <- which.max(R2_adj_LASSO_mat)
  LASSO <- glmnet(B_all[[i]], y_new[[i]], alpha = 1, lambda = regularizador[lamb_op]/(2*nrow(B_all[[i]])), standardize = F, intercept = F)
  betaLASSO <- as.numeric(coef(LASSO)[-1])
  estimativa_funcional_lasso <- B_all[[i]]%*%betaLASSO
  estimativa_funcional_LASSO <- cbind(estimativa_funcional_LASSO, estimativa_funcional_lasso)
  #print(lamb_op)
  R2_lasso <- c(R2_lasso, R2_adj_LASSO_mat[which.max(R2_adj_LASSO_mat)])
  R2_lasso
  
}




# Bayesian LASSO ----
# 2 chains
df.blasso <- NULL
for(i in 1:6){
  BL1=blasso(X = B_all[[i]], y = y_new[[i]], T = 10000, thin = 1, 
             beta=rep(-1, K_all[i]),lambda2 = 1, s2 = 1,
             RJ = F, rao.s2=F, icept = T, normalize=F)
  BL2=blasso(X = B_all[[i]], y = y_new[[i]], T = 10000, thin = 1, 
             beta=rep(1, K_all[i]),lambda2 = 5, s2 = 5,
             RJ = F, rao.s2 = F, icept = T, normalize = F)
  
  burnIn=5000
  df.blasso1=as.data.frame(cbind(lambda=BL1$lambda2,s2=BL1$s2,beta0 = BL1$mu, beta=BL1$beta))
  df.blasso1=df.blasso1[(burnIn+1):nrow(df.blasso1),]
  df.blasso2=as.data.frame(cbind(lambda=BL2$lambda2,s2=BL2$s2,beta0 = BL1$mu, beta=BL2$beta))
  df.blasso2=df.blasso2[(burnIn+1):nrow(df.blasso2),]
  
  thin=50
  df.blasso1 = df.blasso1[rep(c(T,rep(F,thin-1)),length.out=dim(df.blasso1)[1]),]
  df.blasso2 = df.blasso2[rep(c(T,rep(F,thin-1)),length.out=dim(df.blasso2)[1]),]
  df.blasso[[i]] = rbind(df.blasso1,df.blasso2)
}

estimativa_BLASSO <- NULL
R2_adj_BL=vector()
beta_BL <- NULL
for(i in 1:6){
  beta_BL[[i]]=apply(df.blasso[[i]], 2, mean)
  estimativa_funcionalBL=cbind(rep(1,365),B_all[[i]])%*%as.matrix(beta_BL[[i]][-c(1,2)])
  estimativa_BLASSO <- cbind(estimativa_BLASSO, estimativa_funcionalBL)
  
  RSS=sum((y_new[[i]]-estimativa_funcionalBL)^2)
  TSS=sum((y_new[[i]]- mean(y_new[[i]]))^2)
  K_model=sum(abs(beta_BL[[i]][-c(1,2,3)])>0)
  R2_adj_BL[i] = 1-(RSS/(length(y_new[[i]])-K_model))/(TSS/(length(y_new[[i]])-1))
}

# VEM 

out_all <- list(out_tempK20_year, out_tempK20_year, out_tempK30_year, out_tempK30_year, out_tempK30_year, out_tempK30_year)

R2_adj_VEM <- sapply(1:6, function(i){
  out <- out_all[[i]]
  seq_values <- lapply(c(seq(1, 6*K_all[i], K_all[i])), function(x){seq(x,x+K_all[i]-1)})
  estimative_VEM <- (out$mu_beta[seq_values[[i]]]*ifelse(out$p[seq_values[[i]]] > 0.50, 1, 0))%*%t(B_all[[i]])*sds[[i]]
  
  RSS <- sum((y[[i]] - estimative_VEM)^2)
  TSS <- sum((y[[i]] - mean(y[[i]]))^2)
  K_model <- sum(ifelse(out$p[seq_values[[i]]] > 0.50, 1, 0))
  R2_adj_VEM <- 1 - (RSS/(length(y[[i]])-K_model))/(TSS/(length(y[[i]]-1)))
  R2_adj_VEM
})
#0.9955794 0.9963907 0.9944341 0.9964815 0.9933751 0.9966489

# R2 table ----

data.frame(Curve = c(1,2,3,4,5,6), VEM = round(R2_adj_VEM,4), RS = round(R2a_rs,4), BLASSO = round(R2_adj_BL,4), LASSO = round(R2_lasso,4))

# Plot - VEM and others

# CB for VEM according to the optimal K

K <- 30
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt), norder = 4, nbasis = K)
B <- lapply(1:6, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

CB_VEM <- lapply(1:6, function(i){cb_vem(m = 6, K = 30, B = B, i = i, out_tempK30_year)})

K <- 20
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt), norder = 4, nbasis = K)
B <- lapply(1:6, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

CB_VEM[[1]] <- cb_vem(m = 6, K = 20, B = B, i = 1, out_tempK20_year)
CB_VEM[[2]] <- cb_vem(m = 6, K = 20, B = B, i = 2, out_tempK20_year)

# Plot ----
color = c("orange", "purple", "green", "blue")
color = c("orange", "purple", "green","blue", "red", "black")
{
  pdf(file = sprintf("%s_%d.pdf", "plot_temp_mean_curves", 1), width = 12.5, height = 7.2)
  #plot(y_new[[1]]*sds[[1]]~Xt, type = "l", col = "gray", xlab = "Day", ylab = "Temperature", ylim = c(-19.3, 20.6))
  plot(y_new[[1]]*sds[[1]]~Xt, type = "l", col = "orange", xlab = "Day", ylab = "Temperature", ylim = c(-19.3, 22.5))
  lines(y_new[[2]]*sds[[2]]~Xt, type = "l", col = "purple")
  lines(y_new[[3]]*sds[[3]]~Xt, type = "l", col = "green")
  lines(y_new[[4]]*sds[[4]]~Xt, type = "l", col = "blue")
  lines(y_new[[5]]*sds[[5]]~Xt, type = "l", col = "red")
  lines(y_new[[6]]*sds[[6]]~Xt, type = "l", col = "black")
  for(i in 1:6){
    #lines(y_new[[i]]*sds[[i]]~Xt, type = "l", col = "gray")
    out <- out_all[[i]]
    seq_values <- lapply(c(seq(1, 6*K_all[i], K_all[i])), function(x){seq(x,x+K_all[i]-1)})
    lines(Xt, y = (out$mu_beta[seq_values[[i]]]*ifelse(out$p[seq_values[[i]]] > 0.50, 1, 0))%*%t(B_all[[i]])*sds[[i]], lwd = 2, col=color[i] , type = "l", ylab = expression(g[t]), xlab = expression(t))
  }
  
  legend("topright", legend = cities, col = c("orange", "purple", "green","blue", "red", "yellow"), lty = c(1, 1, 1, 1, 1, 1), bty = "n")
  
  
  dev.off()
}

library(scales)
{
  pdf(file = sprintf("%s_%d.pdf", "plot_temp_all_VEM", 1), width = 10, height = 5)
  #seq_values <- lapply(c(seq(1, 4*K, K)), function(x){seq(x,x+K-1)})
  par(mfrow=c(1,1), mar = c(5, 4.1, 2, 0.1), cex = 1) 
  for(i in 1:6){
    out <- out_all[[i]]
    seq_values <- lapply(c(seq(1, 6*K_all[i], K_all[i])), function(x){seq(x,x+K_all[i]-1)})
    plot(y[[i]]~Xt,  cex = 0.5, pch = 16, col = alpha("black", 0.3), ylab = "Temperature", xlab = "day", ylim = c(-19.3, 22.5), main = cities[i])
    #, ylim = c(-19.3, 20.6)
    lines(Xt, y = (out$mu_beta[seq_values[[i]]]*ifelse(out$p[seq_values[[i]]] > 0.50, 1, 0))%*%t(B_all[[i]])*sds[[i]], lwd = 2, col= "black", type = "l", ylab = expression(g[t]), xlab = expression(t))
    #lines(predict(smooth.spline(Xt, y[[i]]), Xt), lwd = 2, col = "blue")
    #lines(Xt, estimativa_funcional_LASSO[,i]*sds[[i]], lwd = 2, col = "blue", lty = 2)
    #lines(Xt, estimativa_BLASSO[,i]*sds[[i]] , lwd = 2, col = "green")
    #lines(Xt, fitted(lm(y_new[[i]] ~ B[[i]] - 1))*sds[[i]], lwd = 2, col = "red", lty = 3)
    #lines(Xt, CB_RS[[i]][,1], col = "red", lty = 2)
    #lines(Xt, CB_RS[[i]][,2], col = "red", lty = 2)
    #lines(Xt, CB_SS[[i]][,1], col = "blue", lty = 2)
    #lines(Xt, CB_SS[[i]][,2], col = "blue", lty = 2)
    # lines(Xt, CB_VEM[[i]][,1]*sds[[i]], col = "black", lty = 2)
    #  lines(Xt, CB_VEM[[i]][,2]*sds[[i]], col = "black", lty = 2)
    polygon(c(Xt, rev(Xt)), c(CB_VEM[[i]][,1]*sds[[i]], rev(CB_VEM[[i]][,2]*sds[[i]])), col = "#00000030", border = NA)
  }
  #print(ifelse(out_mcycle_K20$p[seq_values[[i]]] > 0.50, 1, 0))
  dev.off()   
}

{
  pdf(file = sprintf("%s_%d.pdf", "plot_temp_all", 1), width = 10, height = 5)
  par(mfrow=c(1,1), mar = c(5, 4.1, 2, 0.1), cex = 1) 
  #plot(y[[i]]~Xt,  cex = 1, pch = 16, col = alpha("black", 0.3), ylab = "Temperature", xlab = "day")
  
  for(i in 1:6){
    plot(y[[i]]~Xt,  cex = 0.5, pch = 16, col = alpha("black", 0.3), ylab = "Temperature", xlab = "day", ylim = c(-19.3, 22.5), main = cities[i])
    seq_values <- lapply(c(seq(1, 6*K_all[i], K_all[i])), function(x){seq(x,x+K_all[i]-1)})
    out <- out_all[[i]]
    lines(Xt, y = (out$mu_beta[seq_values[[i]]]*ifelse(out$p[seq_values[[i]]] > 0.50, 1, 0))%*%t(B_all[[i]])*sds[[i]], lwd = 2, col= "black", type = "l", ylab = "Temperature", xlab = "day", ylim=c(-19.3, 20.6))
    
    #lines(predict(smooth.spline(Xt, y[[i]]), Xt), lwd = 2, col = "blue")
    lines(Xt, estimativa_funcional_LASSO[,i]*sds[[i]], lwd = 2, col = "blue", lty = 2)
    lines(Xt, estimativa_BLASSO[,i]*sds[[i]], lwd = 2, col = "green")
    lines(Xt, fitted(lm(y_new[[i]] ~ B_all[[i]] - 1))*sds[[i]], lwd = 2, col = "red", lty = 2)
    #lines(Xt, y[[i]],  cex = 1, pch = 16, col = col = alpha("black", 0.3), ylab = "Temperature", xlab = "day")
    #lines(Xt, CB_RS[[i]][,1], col = "red", lty = 2)
    #lines(Xt, CB_RS[[i]][,2], col = "red", lty = 2)
    #lines(Xt, CB_SS[[i]][,1], col = "blue", lty = 2)
    #lines(Xt, CB_SS[[i]][,2], col = "blue", lty = 2)
    #polygon(c(Xt, rev(Xt)), c(CB_VEM[[i]][,1]*sds[[i]], rev(CB_VEM[[i]][,2]*sds[[i]])), col = "#00000020", border = NA)
    legend("topright", legend = c("Proposed method", "Regression splines", "LASSO", "Bayesian LASSO"), col = c("black", "red", "blue","green"), lty = c(1, 1, 1, 1), bty = "n")
    print(ifelse(out$p[seq_values[[i]]] > 0.50, 1, 0))
  }
  dev.off()   
}
RS_fit <- fitted(lm(y_new[[i]] ~ B_all[[i]] - 1))
summer_t = 150:250
{
  pdf(file = sprintf("%s_%d.pdf", "plot_temp_all_summer", 1), width = 10, height = 5)
  par(mfrow=c(1,1), mar = c(5, 4.1, 2, 0.1), cex = 1) 
  #plot(y[[i]]~Xt,  cex = 1, pch = 16, col = alpha("black", 0.3), ylab = "Temperature", xlab = "day")
  
  for(i in 1:6){
    RS_fit <- fitted(lm(y_new[[i]] ~ B_all[[i]] - 1))
    #plot(Xt[summer_t],y[[i]][summer_t],  cex = 0.5, pch = 16, col = alpha("black", 0.3), ylab = "Temperature", xlab = "day", ylim = c(10, 25), main = cities[i])
    seq_values <- lapply(c(seq(1, 6*K_all[i], K_all[i])), function(x){seq(x,x+K_all[i]-1)})
    out <- out_all[[i]]
    plot(Xt[summer_t], y = (out$mu_beta[seq_values[[i]]]*ifelse(out$p[seq_values[[i]]] > 0.50, 1, 0))%*%t(B_all[[i]][summer_t,])*sds[[i]], lwd = 2, col= "black", type = "l", ylab = "Temperature", xlab = "day", ylim=c(10, 22), main = cities[i])
    
    #lines(predict(smooth.spline(Xt, y[[i]]), Xt), lwd = 2, col = "blue")
    lines(Xt[summer_t], estimativa_funcional_LASSO[summer_t,i]*sds[[i]], lwd = 2, col = "blue", lty = 4)
    lines(Xt[summer_t], estimativa_BLASSO[summer_t,i]*sds[[i]], lwd = 2, col = "green")
    lines(Xt[summer_t], RS_fit[summer_t]*sds[[i]], lwd = 2, col = "red", lty = 2)
    #points(y[[i]]~Xt)
    #points(Xt[zoom_t], y[[i]][zoom_t])
    legend("topright", legend = c("Proposed method", "Regression splines", "LASSO", "Bayesian LASSO"), col = c("black", "red", "blue","green"), lty = c(1, 1, 1, 1), bty = "n")
    #print(ifelse(out_mcycle_K20$p[seq_values[[i]]] > 0.50, 1, 0))
  }
  
  dev.off()
}  



