# computing GCV
library(faraway)
data(mcycle)
plot(accel~jitter(mcycle$times,amount=0.9),mcycle)
y = list(mcycle$accel)
Xt <- jitter(mcycle$times)

K <- 20
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt),norder = 4, nbasis = K)
B <- lapply(1, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)})

seq_values <- lapply(c(seq(1, 1*K, K)), function(x){seq(x,x+K-1)})

out <- out_mcycle_K20

p_i <- out$p[seq_values[[i]]]
Zi_hat <- ifelse(out$p[seq_values[[i]]] > 0.50, 1, 0)
  
EG_i <- p_i*t(B[[1]])

S_k <- E_inv_sigma2(out$delta1, out$delta2)*B[[i]]%*%diag(Zi_hat)%*%out$Sigma_beta[,,i]%*%EG_i%*%solve(calPsi(Xt, Xt, w = out$w))

S_k <- E_inv_sigma2(out$delta1, out$delta2)*B[[i]]%*%diag(Zi_hat)%*%solve(E_inv_sigma2(out$delta1, out$delta2)*(diag(E_inv_tau2(out$lambda1, out$lambda2), K) + EG_i%*%solve(calPsi(Xt, Xt, w = out$w))%*%t(EG_i)))%*%EG_i%*%solve(calPsi(Xt, Xt, out$w))

seq_values <- lapply(c(seq(1, 1*K, K)), function(x){seq(x,x+K-1)})
yhat <- as.numeric((out$mu_beta[seq_values[[i]]]*ifelse(out$p[seq_values[[i]]] > 0.50, 1, 0))%*%t(B[[i]]))
RSS_vb = sum((yhat- y[[1]])^2)

numerador <- t(y[[i]] - B[[i]]%*%(Zi_hat*out$mu_beta[seq_values[[i]]]))%*%(y[[i]] - B[[i]]%*%diag(Zi_hat)%*%out$mu_beta[seq_values[[i]]])

GCV_k <- (1/133)*numerador/((1 - (1/133)*sum(diag(S_k)))^2) 

# K = 30 519.8586
# K = 20 523.6492
# K = 15 528.0622
