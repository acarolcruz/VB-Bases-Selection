elbo <- function(y, ni, B, m, K, iter, delta_1, delta_2, lambda_1, lambda_2, delta1_q, delta2_values, mu_beta_values, lambda1_q, lambda2_values, a1_values, a2_values, Sigma_beta, p_values, mu_ki){
  
 res <- sum(sapply(1:m, function(i){E_log_like_i(y = y, B = B, ni = ni, i = i, iter = iter , delta1_q = delta1_q, delta2_values = delta2_values, mu_beta_values = mu_beta_values, Sigma_beta = Sigma_beta, p_values = p_values)})) + sum(sapply(1:m, function(i){diff_z_i(i = i, iter = iter , K = K, p_values = p_values, mu_beta_values = mu_beta_values, a1_values = a1_values, a2_values = a2_values)})) + sum(sapply(1:m, function(i){diff_beta_i(i = i, iter = iter, K = K, delta1_q = delta1_q, delta2_values = delta2_values, lambda1_q = lambda1_q, lambda2_values = lambda2_values, mu_beta_values = mu_beta_values, Sigma_beta = Sigma_beta)})) + sum(sapply(1:m, function(i){diff_theta_i(i = i, iter = iter, K = K, a1_values = a1_values, a2_values = a2_values, mu_beta_values = mu_beta_values, mu_ki = mu_ki)})) + diff_sigma2(iter = iter, delta_1 = delta_1, delta_2 = delta_2, delta1_q = delta1_q, delta2_values = delta2_values) + diff_tau2(iter = iter, lambda_1 = lambda_1, lambda_2 = lambda_2, lambda1_q = lambda1_q, lambda2_values = lambda2_values)
 
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