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
  
  Var_zi_betai <- var_z*Sigmai + var_z*(mu_bi%*%t(mu_bi))  + Sigmai*(p_i%*%t(p_i)) 
  
  #t(mu_bi)%*%(p_i*(1-p_i)*t(B[[i]]^2))%*%mu_bi + t(p_i*t(B[[i]]))%*%Sigmai%*%(p_i*t(B[[i]])) + sum(diag(t(Sigmai)%*%(p_i*(1-p_i)*t(B[[i]]^2))))
  
  
  
  t(y[[i]] - t(p_i*t(B[[i]]))%*%mu_bi)%*%solve(psi)%*%(y[[i]] - t(p_i*t(B[[i]]))%*%mu_bi) + sum(diag(solve(psi)%*%(B[[i]]%*%Var_zi_betai%*%t(B[[i]])))) 
    
  #
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


calPsi <- function(t, s, w=1) {
  Psi <- matrix(rep(0, length(t)*length(s)), nrow=length(t))
  for (i in 1:nrow(Psi)) {
    for (j in 1:ncol(Psi)) {
      Psi[i,j] <- exp(-2*abs(t[i]-s[j])*(1/w))
      #Psi[i,j] <- exp(-2*w*abs(t[i]-s[j]))
    }
  }
  return(Psi)
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

dev_psi <- function(t, s, w=1) {
  dev_Psi <- matrix(rep(0, length(t)*length(s)), nrow=length(t))
  for (i in 1:nrow(dev_Psi)) {
    for (j in 1:ncol(dev_Psi)) {
      dev_Psi[i,j] <- (2*abs(t[i]-s[j])/(w^2))*exp(-2*abs(t[i]-s[j])/w)
      #dev_Psi[i,j] <- (-2*abs(t[i]-s[j]))*exp(-2*w*abs(t[i]-s[j]))
    }
  }
  return(dev_Psi)
}
