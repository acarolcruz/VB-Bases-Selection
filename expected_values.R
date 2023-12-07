


E_inv_tau2 <- function(lambda1, lambda2){
  lambda1/lambda2
}


E_inv_sigma2 <- function(delta1, delta2){
  delta1/delta2
}

E_square_i <- function(B, i, y, mu, Sigma, p, iter){
  ni <- unlist(lapply(y,length))
  p_i <- p[iter - 1, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu))]
  mu_bi <- mu[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu))]
  var_bi <- diag(Sigma[[iter]][,grep(paste0("^beta_[0-9]+_",i,"$"),colnames(mu))])
  Vi <- sum(sapply(1:ni[i], function(j){sum((var_bi*p_i*(1-p_i) + var_bi*p_i^2 + p_i*(1-p_i)*mu_bi^2)*B[[i]][j,]^2)})) + sum(sapply(1:ni[i], function(j){(y[[i]][j] - sum(mu_bi*p_i*B[[i]][j,]))^2})) 
  
  return(Vi)
} 

E_sumk_beta_i <- function(i, mu, Sigma, iter){
  mu_bi <- mu[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu))]
  var_bi <- diag(Sigma[[iter]][,grep(paste0("^beta_[0-9]+_",i,"$"),colnames(mu))])
  
  return(sum(var_bi + mu_bi^2))
}


E_Ct_C <- function(z, i, p, mu, Sigma, B, y, k, K, iter){
  
  ni <- unlist(lapply(y,length))
  
  #Selecting indices in i that are not k, i.e, i1, ...., ik-1
  #ind_inotk <- grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu))[-which(grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu)) == grep(paste0("^beta_",k,"_",i,"$"), colnames(mu)))] 
  ind_inotk <- seq(1:K)[-k]
  
  p_i_nok <- p[ind_inotk]
  mu_bi_nok <- mu[iter, ind_inotk]
  var_bi_nok <- diag(Sigma[[iter]][,ind_inotk])
  
  mu_bik <- mu[iter, grep(paste0("^beta_",k,"_",i,"$"), colnames(mu))]
  p_ki <- p[k]
  var_bik <- diag(Sigma[[iter]][,grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu))])[k]
  
  CtC <- sum(sapply(1:ni[i], function(j){sum((var_bi_nok*p_i_nok*(1-p_i_nok) + var_bi_nok*p_i_nok^2 + p_i_nok*(1-p_i_nok)*mu_bi_nok^2)*B[[i]][j, -k]^2)})) + 
    sum(var_bik*p_ki^2*B[[i]][, k]^2) + sum(sapply(1:ni[i], function(j){(y[[i]][j] - sum(mu_bi_nok*p_i_nok*B[[i]][j,-k]) - z*mu_bik*B[[i]][j, k])^2})) 
  
  return(CtC)
  
}  
  
  
E_log_theta_ki <- function(i, k, a1, a2){
  a1_ki <- a1
  a2_ki <- a2
  
  return(digamma(a1_ki) - digamma(a1_ki + a2_ki))
}


E_log_theta_ki_c <- function(i, k, a1, a2){
  a1_ki <- a1
  a2_ki <- a2
  
  return(digamma(a2_ki) - digamma(a1_ki + a2_ki))
}

  