# Elbo calculation

E_log_like_i_corr <- function(y, ni, B, i, iter, delta1_q, delta2_values, mu_beta_values, Sigma_beta, p_values, psi){
  
  delta2_q <- delta2_values[iter]
  
  mu_i_q <- mu_beta_values[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))]
  sigma2_i_q <- diag(Sigma_beta[,,i])
  
  p_qi <- p_values[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))]
  
  
  res <- (-ni[i]/2)*(log(2*pi)+E_log_sigma2(delta1_q, delta2_q)) -0.5*log(det(psi)) -0.5*E_inv_sigma2(delta1_q, delta2_q)*E_square_gi(B = B, i = i, y = y, mu = mu_beta_values, Sigma = Sigma_beta, p = p_values, iter = iter, psi = psi)
  
  #cat("E_log_like_i:", res, "\n")
  return(res)
}


diff_z_i <- function(i, iter, K, p_values, mu_beta_values, a1_values, a2_values){
  
  p_qi <- p_values[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))]
  a1_qi <- a1_values[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))] 
  a2_qi <- a2_values[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))] 
  
  
  res <- c()
  for(k in 1:K){
    if((p_qi[k] == 1) | (p_qi[k] == 0)){
      res <- c(res, 0)
    } else{res <- c(res, p_qi[k]*(E_log_theta_ki(a1_qi[k], a2_qi[k]) - log(p_qi[k])) +(1-p_qi[k])*(E_log_theta_ki_c(a1_qi[k], a2_qi[k])-log(1-p_qi[k])))}
    
  }
  #cat("diff_z_i:", sum(res), "\n")
  return(sum(res))
}


diff_beta_i <- function(i, iter, K, delta1_q, delta2_values, lambda1_q, lambda2_values, mu_beta_values, Sigma_beta){
  
  delta2_q <- delta2_values[iter]
  
  lambda2_q <- lambda2_values[iter]
  
  mu_i_q <- mu_beta_values[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))]
  
  Sigma_i_q <- Sigma_beta[,,i]
  sigma2_i_q <- diag(Sigma_i_q)
  
  
  res <- -(K/2)*(E_log_sigma2(delta1_q, delta2_q) + E_log_tau2(lambda1_q, lambda2_q)) - sum(sapply(1:K, function(k){((sigma2_i_q[k] + mu_i_q[k]^2)/2)}))*E_inv_sigma2(delta1_q, delta2_q)*E_inv_tau2(lambda1_q, lambda2_q) -(-0.5*log(det(Sigma_i_q)) - 0.5*K)
  #cat("diff_beta_i:", res, "\n")
  return(res) 
}


diff_theta_i <- function(i, iter, K, a1_values, a2_values, mu_beta_values, mu_ki){
  
  a1_qi <- a1_values[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))] 
  a2_qi <- a2_values[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))] 
  
  res <- c()
  for(k in 1:K){
    res <- c(res, (mu_ki - a1_qi[k])*E_log_theta_ki(a1_qi[k], a2_qi[k]) + (1-mu_ki - a2_qi[k])*E_log_theta_ki_c(a1_qi[k], a2_qi[k]) + log(gamma(a1_qi)*gamma(a2_qi)))
  }
  sum(res)
  
  #cat("diff_theta_i :", sum(res), "\n")
  return(sum(res))
  
}

diff_sigma2 <- function(iter, delta_1, delta_2, delta1_q, delta2_values){
  
  delta2_q <- delta2_values[iter]
  
  res <-  - delta1_q*log(delta2_q)  + (delta1_q - delta_1)*E_log_sigma2(delta1_q, delta2_q) + (delta2_q - delta_2)*E_inv_sigma2(delta1_q, delta2_q)
  #cat("diff_sigma2_i :", res, "\n")
  return(res)
}

diff_tau2 <- function(iter, lambda_1, lambda_2, lambda1_q, lambda2_values){
  
  lambda2_q <- lambda2_values[iter]
  
  res <-  - lambda1_q*log(lambda2_q) + (lambda1_q - lambda_1)*E_log_tau2(lambda1_q, lambda2_q) + (lambda2_q - lambda_2)*E_inv_tau2(lambda1_q, lambda2_q)
  #cat("diff_tau2_i :", res, "\n")
  return(res)
}

