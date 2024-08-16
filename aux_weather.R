# aux for weather data

cb_rs <- function(y,B,i){
  fitted_rs <- fitted(lm(y[[i]] ~ B[[i]] - 1))
  sigma2hat <- (sum((y[[i]]-fitted_rs)^2)/(length(Xt)-ncol(B[[i]])))
  
  #Sigma2hat <-  (y[[i]]-fitted_rs)%*%t((y[[i]]-fitted_rs))/(length(y[[i]]) - ncol(B[[i]]))
  #varcoef <- solve(t(B[[i]])%*%solve(Sigma2hat)%*%B[[i]])
  varcoef <- sigma2hat*solve(t(B[[i]])%*%B[[i]])
  varfhat <- B[[i]]%*%varcoef%*%t(B[[i]]) # var (f hat) = var(betaB)= B%*%var(beta)%*%t(B)
  sd_rs_fit <- sqrt(diag(varfhat))
  UB <- fitted_rs + 1.96*sd_rs_fit
  LB <- fitted_rs - 1.96*sd_rs_fit
  
  return(data.frame(LB = LB, UB = UB, fitted_rs))
}

cb_ss <- function(y,i,lambda = 10^(-4.5)){
  
  Bspline.basis <- create.bspline.basis(range = c(min(Xt), max(Xt)),
                                        norder = 4, breaks = sort(unique(Xt)))
  fdParobj <- fdPar(Bspline.basis, Lfdobj=2, lambda=lambda)
  lambda_opt <- lambda2gcv(log10(lambda), Xt, y[[i]], fdParobj, wtvec=rep(1,length(Xt)))
  yfdPar <- fdPar(Bspline.basis, Lfdobj=2, lambda=lambda_opt)
  yfd <- smooth.basis(Xt, y[[i]], yfdPar)$fd 
  yhat <- eval.fd(Xt, yfd)
  
  #loglamout <- plotGCVRMSE.fd(-6, -3, 0.25, x, argvals = Xt, y, fdParobj)
  
  B_ss <- getbasismatrix(Xt, Bspline.basis, nderiv = 0)
  R <- smooth.basis(Xt, y[[i]], yfdPar)$penmat 
  p <- smooth.basis(Xt, y[[i]], yfdPar)$df #effective number of parameters 
  s_df <- length(y[[i]])-p
  sigma2hat <- sum((y[[i]]-yhat)^2)/s_df
  
  A_star <- solve(t(B_ss)%*%B_ss + lambda_opt*R)%*%t(B_ss) 
  varcoef <- sigma2hat*A_star%*%t(A_star)
  varfhat <- B_ss%*%varcoef%*%t(B_ss)
  sd_ss_fit <- sqrt(diag(varfhat))
  UB <- yhat + 1.96*sd_ss_fit
  LB <- yhat - 1.96*sd_ss_fit
  
  return(data.frame(LB = LB, UB = UB, yhat = yhat))
}

cb_vem <- function(m, K, i, out, B, cl_l = 0.025, cl_u = 0.975){
  seq_values <- lapply(c(seq(1, m*K, K)), function(x){seq(x,x+K-1)})
  
  LL <- NULL
  UL <- NULL
  estimates <- NULL
  
  betas_sample <- MASS::mvrnorm(200, mu = out$mu_beta[seq_values[[i]]],  Sigma = out$Sigma_beta[,,i])
  z_sample <- matrix(rbinom(K*200, 1, prob = out$p[seq_values[[i]]]), ncol = K, byrow = TRUE)
  
  estimates <- cbind(estimates, (betas_sample*z_sample))
  
  curve_f <- list()
  for(s in 1:200){
    curve_f[[s]] <- estimates[s,]%*%t(B[[i]])
  }  
  
  LL <- apply(do.call(rbind,curve_f), 2, function(i)quantile(i,probs = c(cl_l, cl_u)))[1,]
  UL <- apply(do.call(rbind,curve_f), 2, function(i)quantile(i,probs = c(cl_l, cl_u)))[2,]
  
  return(data.frame(LL,UL))
}


Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


res_gibbs <- function(folder, K, m, maxIter = 10000, thin = 50){
  #burn-in and thinning
  load(paste0(folder,"chain",1, "/Gibbsresults.RData"))
  index <- seq(from = (maxIter/2)+1, to = maxIter, by = thin)
  res1 <- list(pars_mat = out[[1]][index,], z_mat = out[[2]][index,])
  save(res1, file = paste0(folder,"chain",1, "/IterGibbs_sample.RData"))
  
  #burn-in and thinning
  load(paste0(folder, "chain", 2, "/Gibbsresults.RData"))
  index <- seq(from = (maxIter/2)+1, to = maxIter, by = thin)
  res1 <- list(pars_mat = out[[1]][index,], z_mat = out[[2]][index,])
  save(res1, file = paste0(folder, "chain", 2, "/IterGibbs_sample.RData"))
  
  # chain1 ----
  #res 1: sample parameters after burnin 50%
  load(paste0(folder,"chain",1,"/IterGibbs_sample.RData"))
  
  out_chain1 <- res1
  
  # chain2 ----
  
  #res 1: sample parameters after burnin 50% and thinning k = 25
  load(paste0(folder,"chain",2, "/IterGibbs_sample.RData"))
  
  out_chain2 <- res1
  
  # generating pos results per sim ----
  
  out <- list(pars = rbind(out_chain1$pars_mat, out_chain2$pars_mat), z_mat = rbind(out_chain1$z_mat, out_chain2$z_mat))
  
  # posterior means and modes
  out2 <- list(pars = apply(out$pars, 2, mean),
               z_mat = apply(out$z_mat, 2, Mode))
  
  out_final <- list(betas = out2$pars[grep("^beta_",names(out2[[1]]))], thetas = out2$pars[grep("^theta_",names(out2[[1]]))], sigma2 = out2$pars[grep("^sigma2",names(out2[[1]]))], tau2 = out2$pars[grep("^tau2",names(out2[[1]]))], Z = out2$z_mat)
  
  for(g in names(out_final)){
    write(out_final[[g]], file = paste0(folder,  g, ".txt"), ncolumns = length(out_final[[g]]))
  }
  
  # Run this code after checking convergence and merging chain results
  beta <- c()
  z <- c()
  
  beta <- scan(paste0(folder, "betas.txt"))
  z <- scan(paste0(folder, "Z.txt"))
  
  
  tempo <- c()
  for(chain in 1:2){
    tempo[chain] <- scan(paste0(folder,"chain",1, "/tempo.txt"))
  }
  
  
  betas <- out$pars[,grep("^beta_",colnames(out[[1]]))]
  zs <- out$z_mat
  
  etakis <- (betas*zs)
  
  # credible band for MCMC
  seq_values <- lapply(c(seq(1, m*K, K)), function(x){seq(x,x+K-1)})
  
  curves <- NULL
  for(i in 1:m){
    curve_f <- c()
    for(s in 1:200){
      curve_f <- rbind(curve_f, (betas*zs)[s,seq_values[[i]]]%*%t(B[[i]]))
    } 
    curves[[i]] <- curve_f
  }
  
  eta_plot <- apply(etakis, 2, mean)
  
  CB_gibbs <<- lapply(1:m, function(i){apply(curves[[i]],2,function(x)quantile(x,probs = c(0.025,0.975)))})
  
  return(list(eta = eta_plot, z = z, tempo = tempo))
}


gcv <- function(i, out, y, B, seq_values){
  p_i <- out$p[seq_values[[i]]]
  Zi_hat <- ifelse(out$p[seq_values[[i]]] > 0.50, 1, 0)
  
  EG_i <- p_i*t(B[[1]])
  
  S_k <- E_inv_sigma2(out$delta1, out$delta2)*B[[1]]%*%diag(Zi_hat)%*%out$Sigma_beta[,,i]%*%EG_i%*%solve(calPsi(Xt, Xt, w = out$w))
  
  S_k <- E_inv_sigma2(out$delta1, out$delta2)*B[[1]]%*%diag(Zi_hat)%*%solve(E_inv_sigma2(out$delta1, out$delta2)*(diag(E_inv_tau2(out$lambda1, out$lambda2), K) + EG_i%*%solve(calPsi(Xt, Xt, w = out$w))%*%t(EG_i)))%*%EG_i%*%solve(calPsi(Xt, Xt, out$w))
  
  yhat <- as.numeric((out$mu_beta[seq_values[[i]]]*ifelse(out$p[seq_values[[i]]] > 0.50, 1, 0))%*%t(B[[1]]))
  RSS_vb = sum((yhat- y[[1]])^2)
  
  numerador <- t(y[[i]] - B[[1]]%*%(Zi_hat*out$mu_beta[seq_values[[i]]]))%*%(y[[i]] - B[[1]]%*%diag(Zi_hat)%*%out$mu_beta[seq_values[[i]]])
  
  GCV_k <- (1/length(Xt))*numerador/((1 - (1/length(Xt))*sum(diag(S_k)))^2) 
  return(GCV_k)
}  
