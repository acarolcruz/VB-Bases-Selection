#Auxiliary functions for VB with correlated errors



# Expected values -----
E_inv_tau2 <- function(lambda1, lambda2){
  lambda1/lambda2
}


E_inv_sigma2 <- function(delta1, delta2){
  delta1/delta2
}

E_log_sigma2 <- function(delta1, delta2){
  log(delta2) - digamma(delta1)
}

E_log_tau2 <- function(lambda1, lambda2){
  log(lambda2) - digamma(lambda1)
}


E_square_gi <- function(B, i, y, mu, Sigma, p, iter, psi){
  p_i <- p[iter - 1, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu))]
  mu_bi <- mu[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu))]
  Sigmai <- Sigma[,,i]
  #var_bi <- diag(Sigma[,,i])
  #var_beta_z <- as.numeric(t(var_bi)%*%(p_i*(1-p_i)*t(B[[i]]^2)) + t(var_bi)%*%((p_i*t(B[[i]]))^2) + t(mu_bi^2)%*%(p_i*(1-p_i)*t(B[[i]]^2)))
  
  var_z <- diag(p_i*(1-p_i))
  
  Var_zi_betai <- var_z*Sigmai + var_z*(mu_bi%*%t(mu_bi)) + Sigmai*(p_i%*%t(p_i)) 
  
  #t(mu_bi)%*%(p_i*(1-p_i)*t(B[[i]]^2))%*%mu_bi + t(p_i*t(B[[i]]))%*%Sigmai%*%(p_i*t(B[[i]])) + sum(diag(t(Sigmai)%*%(p_i*(1-p_i)*t(B[[i]]^2))))
  
  
  #t(p_i*t(B[[i]]))%*%mu_bi
  #if(det(psi) < 4.387555e-48) browser()
  res <- t(y[[i]] - B[[i]]%*%(mu_bi*p_i))%*%solve(psi)%*%(y[[i]] - B[[i]]%*%(mu_bi*p_i)) + sum(diag(solve(psi)%*%(B[[i]]%*%Var_zi_betai%*%t(B[[i]])))) 
  return(res)
}

E_sumk_beta_i_corr <- function(i, mu, Sigma, iter){
  mu_bi <- mu[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu))]
  var_bi <- diag(Sigma[,,i])
  
  return(sum(var_bi + mu_bi^2))
}


E_square_beta_i <- function(z, i, p, mu, Sigma, B, y, k, K, iter, psi){
  
  ni <- unlist(lapply(y,length))
  
  #Selecting indices in i that are not k, i.e, i1, ...., ik-1
  #ind_inotk <- grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu))[-which(grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu)) == grep(paste0("^beta_",k,"_",i,"$"), colnames(mu)))] 

  mu_bi <- mu[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu))]
  Sigmai <- Sigma[,,i]
  
  # vector with probabilities when k-th position is z
  pi_notk <- p
  pi_notk[k] <- z
  
  var_zi_notk <- pi_notk*(1-pi_notk)
  var_zi_notk[k] <- 0
  
  Var_zi_notk <- diag(var_zi_notk)
  
  #var_zi_betai_notk <- t((pi_notk*t(B[[i]]))^2)%*%diag(Sigmai) + t((var_zi_notk*t(B[[i]]))^2)%*%diag(Sigmai) + t((var_zi_notk*t(B[[i]]))^2)%*%mu_bi
  
  #var_zi_betai_notk <- as.numeric(t(diag(Sigmai))%*%(pi_notk*t(B[[i]])) + t(diag(Sigmai))%*%(var_zi_notk*(t(B[[i]])^2)) + t(mu_bi^2)%*%(var_zi_notk*(t(B[[i]])^2)))
  
  Var_zi_betai_notk <- Var_zi_notk*Sigmai + Var_zi_notk*(mu_bi%*%t(mu_bi)) + Sigmai*(pi_notk%*%t(pi_notk))
  
  #is.singular.matrix
  #pacote matrixcalc
  
  #print(is.singular.matrix(Var_zi_betai_notk))
  
  #diag(diag(Var_zi_betai_notk))
  
  #Esquare_i_beta <- t(y[[i]] - t(pi_notk*t(B[[i]]))%*%mu_bi)%*%solve(psi)%*%(y[[i]] - t(pi_notk*t(B[[i]]))%*%mu_bi) + sum(diag(solve(psi)%*%(B[[i]]%*%Var_zi_betai_notk%*%t(B[[i]]))))
  
  Esquare_i_beta <- t(y[[i]] - t(pi_notk*t(B[[i]]))%*%mu_bi)%*%solve(psi)%*%(y[[i]] - t(pi_notk*t(B[[i]]))%*%mu_bi) + sum(diag(solve(psi)%*%(B[[i]]%*%Var_zi_betai_notk%*%t(B[[i]]))))
  
  #sum(diag(t(Var_zi_betai_notk)%*%Sigmai)) + t(pi_notk)%*%Sigmai%*%pi_notk
  
  #um(sapply(1:ni[i], function(j){sum((var_bi_nok*p_i_nok*(1-p_i_nok) + var_bi_nok*(p_i_nok^2) + p_i_nok*(1-p_i_nok)*(mu_bi_nok^2))*B[[i]][j, -k]^2)})) + sum(var_bik*(z^2)*(B[[i]][, k]^2))
  
  return(Esquare_i_beta)
  
}  


E_log_theta_ki <- function(a1_ki, a2_ki){
  
  return(digamma(a1_ki) - digamma(a1_ki + a2_ki))
}


E_log_theta_ki_c <- function(a1_ki, a2_ki){
  
  return(digamma(a2_ki) - digamma(a1_ki + a2_ki))
}

# Elbo -----


elbo_corr <- function(y, ni, B, m, K, iter, delta_1, delta_2, lambda_1, lambda_2, delta1_q, delta2_values, mu_beta_values, lambda1_q, lambda2_values, a1_values, a2_values, Sigma_beta, p_values, mu_ki, psi){
  
  res <- sum(sapply(1:m, function(i){E_log_like_i_corr(y = y, B = B, ni = ni, i = i, iter = iter , delta1_q = delta1_q, delta2_values = delta2_values, mu_beta_values = mu_beta_values, Sigma_beta = Sigma_beta, p_values = p_values, psi = psi)})) + sum(sapply(1:m, function(i){diff_z_i(i = i, iter = iter , K = K, p_values = p_values, mu_beta_values = mu_beta_values, a1_values = a1_values, a2_values = a2_values)})) + sum(sapply(1:m, function(i){diff_beta_i(i = i, iter = iter, K = K, delta1_q = delta1_q, delta2_values = delta2_values, lambda1_q = lambda1_q, lambda2_values = lambda2_values, mu_beta_values = mu_beta_values, Sigma_beta = Sigma_beta)})) + sum(sapply(1:m, function(i){diff_theta_i(i = i, iter = iter, K = K, a1_values = a1_values, a2_values = a2_values, mu_beta_values = mu_beta_values, mu_ki = mu_ki)})) + diff_sigma2(iter = iter, delta_1 = delta_1, delta_2 = delta_2, delta1_q = delta1_q, delta2_values = delta2_values) + diff_tau2(iter = iter, lambda_1 = lambda_1, lambda_2 = lambda_2, lambda1_q = lambda1_q, lambda2_values = lambda2_values)
  
  return(res)
}

check_convergence <- function(elbo_c, elbo_prev, convergence_threshold) {
  if(is.null(elbo_prev) == TRUE) {
    return(FALSE)
  }
  else{
    dif <- elbo_c - elbo_prev
    if(abs(dif)  <= convergence_threshold) return(TRUE)
    else return(FALSE)
  }
}  


calPsi <- function(tt, s, w=1) {
  Psi <- matrix(rep(0, length(tt)*length(s)), nrow=length(tt))
  for (i in 1:nrow(Psi)) {
    for (j in 1:ncol(Psi)) {
      #Psi[i,j] <- exp(-(2/w) * (abs(tt[i] - s[j]) / (max(tt) - min(tt)))) #Gabriel's
      Psi[i,j] <- exp(- w * abs(tt[i] - s[j]) / (max(tt) - min(tt))) #Ronaldo's
      #Psi[i,j] <- exp(-w^2*(tt[i] - s[j])^2) #squared-exponential
    
    }
  }
  return(Psi)
}



# Creating and computing correlation matrix based on a Gaussian Process
calCov <- function(tt, s, sigma=1, w=1) {
  Cov <- sigma^2 * calPsi(tt, s, w)
  return(Cov)
}


# EM ----

# FUN <- function(w, Xt = Xt, n_i = n_i, out = out, Emui = Emui, var_beta_z = var_beta_z){
#   psi <- calPsi(Xt, Xt, w = w)
#   -(-0.5*n_i*(log(2*pi) + E_log_sigma2(out$delta1, out$delta2)) - 0.5*log(det(psi)) - 0.5*(out$delta1/out$delta2)*(t(y[[i]]-Emui)%*%solve(psi)%*%(y[[i]]-Emui) + sum(diag(solve(psi)%*%diag(var_beta_z[1,])))))
# }

# FUN_grad <- function(w=w){
#   Emui <- t(pzi*t(B[[i]]))%*%mu_betai
#   A <- (solve(psi)%*%dev_psi(Xt, Xt, w=w)%*%solve(psi))
#   var_beta_z <- t(sigma_betai + mu_betai^2)%*%(pzi*(1-pzi)*t(B[[i]]^2)) + t(sigma_betai)%*%(pzi*(1-pzi)*t(B[[i]]^2)) + t(mu_betai^2)%*%(pzi*(1-pzi)*t(B[[i]]^2))
#   -(-0.5*sum(diag(solve(psi)%*%dev_psi(Xt, Xt, w=w))) 
#     -0.5*(delta1/delta2)*(t(y[[i]]-Emui)%*%A%*%(y[[i]]-Emui) - sum(diag(A%*%diag(var_beta_z[1,])))))
# }

# FUN_grad <- function(w=w, Xt, out){
#   psi <- calPsi(Xt, Xt, w = w)
#   
#   sum(sapply(1:m, function(i){
#     p_i <- out$p[grep(paste0("^beta_[0-9]+_",i,"$"), names(out$mu_beta))]
#     mu_bi <- out$mu_beta[grep(paste0("^beta_[0-9]+_",i,"$"), names(out$mu_beta))]
#     Sigmai <- out$Sigma_beta[,,i]
#     var_bi <- diag(Sigmai)
#   
#     var_z <- diag(p_i*(1-p_i))
#     
#     Var_zi_betai <- var_z*Sigmai + var_z*(mu_bi%*%t(mu_bi))  + Sigmai*(p_i%*%t(p_i))
# 
#     (-0.5*sum(diag(solve(psi)%*%dev_psi(Xt, Xt, w=w))) + 0.5*E_inv_sigma2(out$delta1, out$delta2)*(t(y[[i]] - t(p_i*t(B[[i]]))%*%mu_bi)%*%solve(psi)%*%dev_psi(Xt, Xt, w=w)%*%solve(psi)%*%(y[[i]] - t(p_i*t(B[[i]]))%*%mu_bi) -0.5*E_inv_sigma2(out$delta1, out$delta2)*sum(diag(-solve(psi)%*%dev_psi(Xt, Xt, w=w)%*%solve(psi)%*%(B[[i]]%*%Var_zi_betai%*%t(B[[i]]))))))}))
# 
# }

dev_psi <- function(tt, s, w=1) {
  dev_Psi <- matrix(rep(0, length(tt)*length(s)), nrow=length(tt))
  for (i in 1:nrow(dev_Psi)) {
    for (j in 1:ncol(dev_Psi)) {
      #dev_Psi[i,j] <- (2*abs(tt[i]-s[j])/(w^2))*exp(-2*abs(tt[i]-s[j])/w) # gabriel's
      dev_Psi[i,j] <- (-abs(tt[i]-s[j]) / (max(tt) - min(tt)))*exp(-w*abs(tt[i]-s[j]) / (max(tt) - min(tt))) # ronaldo's
      #dev_Psi[i,j] <- (-abs(tt[i]-s[j]))*exp(log(w))*exp(-exp(log(w))*abs(tt[i]-s[j])) # ronaldo's (2)
      #dev_Psi[i,j] <- -2*w*(tt[i] - s[j])^2 * exp(-w^2*(tt[i] - s[j])^2) #squared-exponential 
    }
  }
  return(dev_Psi)
}

# checar essa derivada
dev2_psi <- function(tt, s, w=1) {
  dev2_Psi <- matrix(rep(0, length(tt)*length(s)), nrow=length(tt))
  for (i in 1:nrow(dev2_Psi)) {
    for (j in 1:ncol(dev2_Psi)) {
      #dev2_Psi[i,j] <- ((4*(abs((tt[i]-s[j]))^2)/(w^4)) - 4*(abs(tt[i]-s[j])/(w^3)))*exp(-2*abs(tt[i]-s[j])/w) # gabriel's
      dev2_Psi[i,j] <- (abs(tt[i]-s[j]) / (max(tt) - min(tt)))^2*exp(-w*abs(tt[i]-s[j]) / (max(tt) - min(tt))) # ronaldo's
      #dev2_Psi[i,j] <- (-2*exp(-w^2*(tt[i] - s[j])^2)*(tt[i] - s[j])^2)*(1 + 2*w^2*(tt[i] - s[j])^2) #squared-exponential 
      #dev2_Psi[i,j] <- exp(-exp(log(w))*abs(tt[i]-s[j]))*((abs(tt[i]-s[j])*exp(log(w)))^2 - abs(tt[i]-s[j])*exp(log(w))) # ronaldo's (2)
    }
  }
  return(dev2_Psi)
}

# if using uniroot: uncomment this code and comment the other dev_elbo
# dev_elbo <- function(w, m,  delta1_q, delta2_q, Xt,
#                      iter, y, B, p_values, mu_beta_values, Sigma_beta){
# 
# 
#   d_psi <- dev_psi(Xt, Xt, w)
#   d2_psi <- dev2_psi(Xt, Xt, w)
#   psi <- calPsi(Xt, Xt, w)
# 
#   A_matrix <- lapply(1:m, function(i){
#     p_i <- p_values[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))]
#     mu_bi <- mu_beta_values[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))]
#     Sigmai <- Sigma_beta[,,i]
# 
#     var_z <- diag(p_i*(1-p_i))
# 
#     Var_zi_betai <- var_z*Sigmai + var_z*(mu_bi%*%t(mu_bi)) + Sigmai*(p_i%*%t(p_i))
# 
#     (y[[i]] - t(p_i*t(B[[i]]))%*%mu_bi)%*%t(y[[i]] - t(p_i*t(B[[i]]))%*%mu_bi) + B[[i]]%*%Var_zi_betai%*%t(B[[i]])
#   })
# 
#   dev_E <- E_inv_sigma2(delta1 = delta1_q, delta2 = delta2_q)
#   s_psi <- solve(psi)
# 
#   dev1 <- -(m / 2) * sum(diag(s_psi %*% d_psi)) +
#     0.5 * dev_E * sum(diag(s_psi %*% d_psi %*% s_psi %*% Reduce('+', A_matrix)))
# 
# 
#   return(dev1)
# }
# 
elbo_omega <- function(w, y, Xt, B, ni, m, K, iter, delta_1, delta_2, lambda_1, lambda_2, delta1_q, delta2_values, mu_beta_values, lambda1_q, lambda2_values, a1_values, a2_values, Sigma_beta, p_values, mu_ki){
   psi <- calPsi(Xt, Xt, w)

   elbo_value <- elbo_corr(y = y, B = B, ni = ni, m = m, K = K, iter = iter, delta_1 = delta_1, delta_2 = delta_2, lambda_1 = lambda_1, lambda_2 = lambda_2, delta1_q = delta1_q, delta2_values = delta2_values, mu_beta_values = mu_beta_values, lambda1_q = lambda1_q, lambda2_values = lambda2_values, a1_values = a1_values, a2_values = a2_values, Sigma_beta = Sigma_beta, p_values = p_values, mu_ki = mu_ki, psi = psi)
   return(elbo_value)
}

dev_elbo <- function(w, y, Xt, B, ni, m, K, iter, delta_1, delta_2, lambda_1, lambda_2, delta1_q, delta2_values, mu_beta_values, lambda1_q, lambda2_values, a1_values, a2_values, Sigma_beta, p_values, mu_ki){
  delta2_q = delta2_values[iter]

  d_psi <- dev_psi(Xt, Xt, w)
  d2_psi <- dev2_psi(Xt, Xt, w)
  psi <- calPsi(Xt, Xt, w)

  A_matrix <- lapply(1:m, function(i){
    p_i <- p_values[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))]
    mu_bi <- mu_beta_values[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))]
    Sigmai <- Sigma_beta[,,i]

    var_z <- diag(p_i*(1-p_i))

    Var_zi_betai <- var_z*Sigmai + var_z*(mu_bi%*%t(mu_bi)) + Sigmai*(p_i%*%t(p_i))

    (y[[i]] - t(p_i*t(B[[i]]))%*%mu_bi)%*%t(y[[i]] - t(p_i*t(B[[i]]))%*%mu_bi) + B[[i]]%*%Var_zi_betai%*%t(B[[i]])
  })

  dev_E <- E_inv_sigma2(delta1 = delta1_q, delta2 = delta2_q)
  s_psi <- solve(psi)

  dev1 <- -(m / 2) * sum(diag(s_psi %*% d_psi)) +
    0.5 * dev_E * sum(diag(s_psi %*% d_psi %*% s_psi %*% Reduce('+', A_matrix)))


  return(dev1)
}

# 
# dev2_elbo <- function(w_initial, y, Xt, B, ni, m, K, iter, delta_1, delta_2, lambda_1, lambda_2, delta1_q, delta2_values, mu_beta_values, lambda1_q, lambda2_values, a1_values, a2_values, Sigma_beta, p_values, mu_ki){
#   
#   delta2_q = delta2_values[iter]
#   
#   d_psi <- dev_psi(Xt, Xt, w_initial)
#   d2_psi <- dev2_psi(Xt, Xt, w_initial)
#   psi <- calPsi(Xt, Xt, w_initial)
#   
#   A_matrix <- lapply(1:m, function(i){
#     p_i <- p_values[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))]
#     mu_bi <- mu_beta_values[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))]
#     Sigmai <- Sigma_beta[,,i]
#     
#     var_z <- diag(p_i*(1-p_i))
#     
#     Var_zi_betai <- var_z*Sigmai + var_z*(mu_bi%*%t(mu_bi)) + Sigmai*(p_i%*%t(p_i)) 
#     
#     (y[[i]] - t(p_i*t(B[[i]]))%*%mu_bi)%*%t(y[[i]] - t(p_i*t(B[[i]]))%*%mu_bi) + B[[i]]%*%Var_zi_betai%*%t(B[[i]])
#   })
#   
#   dev_E <- E_inv_sigma2(delta1 = delta1_q, delta2 = delta2_q)
#   s_psi <- solve(psi)
#   
#   dev2 <- -(m / 2) * sum(diag(
#     -s_psi %*% (d_psi %*% s_psi %*% d_psi + d2_psi))) + 
#     0.5 * dev_E * sum(diag(s_psi %*% ((d_psi) ^ 2 + d2_psi - d_psi %*% s_psi %*% d_psi) %*% (s_psi %*% Reduce('+', A_matrix))))
#   
#   
#   return(dev2)
# }


