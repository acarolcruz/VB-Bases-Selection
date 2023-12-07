#' Implements VB for the bases selection in scalar-on-function functional regression
#' 
#' @param y list containing the values for the outcome variable
#' @param z vector for the bases selection 
#' @param Xt vector of time points (same for all curves)
#' @param B matrix with B-spline values
#' @param K number of bases
#' @param m number of curves
#' @param lambda_1 hyperparameter for the tau2 variance component of betaki coefficient
#' @param lambda_2 hyperparameter for the tau2 variance component  of betaki coefficient
#' @param delta_1 hyperparameter for the sigma2 variance component of betaki coefficient
#' @param delta_2 hyperparameter for the sigma2 variance component  of betaki coefficient
#' @param mu_ki hyperparameter for the distribution of theta
#' @returns parameters for the variational densities and ELBO
#' 


vb_bs <- function(y, B, m = 5, mu_ki = 1/2, lambda_1 = 10^(-10), lambda_2 = 10^(-10), delta_1 = 10^(-10), delta_2 = 10^(-10), maxIter = 1000, K = 10, initial_values){
  
  # Creating structure to save output
  mu_beta_values <- matrix(NA, maxIter, m*K) 
  colnames(mu_beta_values) <- as.vector(outer(paste0("beta_",1:K,"_"),1:m,paste0))
  Sigma_beta <- list()
  p_values <- matrix(NA, maxIter, m*K)
  lambda2_values <- c()
  delta2_values <- c()
  a1_values <-  matrix(NA, maxIter, m*K) 
  a2_values <- matrix(NA, maxIter, m*K) 
  
  L <- rep(-Inf, maxIter)  # Store the lower bounds
  
  p_values[1,] <- initial_values$p
  lambda2_values[1] <- initial_values$lambda2
  delta2_values[1] <- initial_values$delta2
  
  ni <- unlist(lapply(y,length))
  
  # is not being updated through variational inference
  delta1_q <- (sum(ni) + m*K + 2*delta_1)/2
  lambda1_q = (m*K + 2*lambda_1)/2
  
  for(iter in 2:maxIter){
    Sigma_beta[[iter]] <- matrix(NA, K, m*K) 
    
    for(i in 1:m){
      
      p_i <- p_values[iter - 1, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))]
      
      HHt <- list()
      Hy <- list()
      for(j in 1:ni[i]){
        H_ij <- p_i*B[[i]][j,]
        HHt[[j]] <- H_ij%*%t(H_ij)
        Hy[[j]] <- (H_ij)*y[[i]][j]
      }
      
      V_i <- solve(diag(1, K)*E_inv_tau2(lambda1 = lambda1_q, 
                                         lambda2 = lambda2_values[iter - 1]) + Reduce('+',HHt))
      
      #Update Sigma matrix for beta_i
      Sigma_bi_q <- (1/E_inv_sigma2(delta1 = delta1_q, 
                                    delta2 = delta2_values[iter - 1]))*V_i
      
      #Update mean of the beta_i 
      mu_bi_q <- E_inv_sigma2(delta1 = delta1_q, 
                              delta2 = delta2_values[iter - 1])*Sigma_bi_q%*%Reduce('+',Hy)
      
      mu_beta_values[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))] <- mu_bi_q
      Sigma_beta[[iter]][, grep(paste0("^beta_[0-9]+_",i,"$"),colnames(mu_beta_values))] <- Sigma_bi_q
      
    }
    
    #update delta2 and lambda2
    delta2_q <- (sum(sapply(1:m, function(i){E_square_i(B = B, i = i, y = y, mu = mu_beta_values, Sigma = Sigma_beta, p = p_values, iter = iter)})) + 
                   sum(sapply(1:m, function(i){E_sumk_beta_i(i = i, mu = mu_beta_values, Sigma = Sigma_beta,  iter = iter)}))*E_inv_tau2(lambda1 = lambda1_q,lambda2 = lambda2_values[iter - 1]) + 2*delta2_values[1])/2
    
    delta2_values[iter] <- delta2_q
    
    lambda2_q <- (sum(sapply(1:m, function(i){E_sumk_beta_i(i = i, mu = mu_beta_values, Sigma = Sigma_beta,  iter = iter)}))*E_inv_sigma2(delta1 = delta1_q, delta2 = delta2_values[iter]) + 2*lambda2_values[1])/2
    
    lambda2_values[iter] <- lambda2_q
    
    
    # #update a1 and a2
    # a1_q_vector <- c()
    # a2_q_vector <- c()
    # for(i in 1:m){
    #   p_i <- p_values[1, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))]
    #   for(k in 1:K){
    #     a1_q <- p_i[k] + mu_ki
    #     a2_q <- 2 - p_i[k] - mu_ki
    #     a1_q_vector <- c(a1_q_vector, a1_q) 
    #     a2_q_vector <- c(a2_q_vector, a2_q)
    #   }
    # }
    # 
    # a1_values[iter,] <- a1_q_vector
    # a2_values[iter,] <- a2_q_vector
    
    #Update p_ki
    
    for(i in 1:m){
      p_q <- c()
      a1_q_vector <- c()
      a2_q_vector <- c()
      
      # p values for the previous iteration
      p_i <- p_values[iter - 1, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))]
      
      # p value for the current iteration
      p_star <- p_values[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))]
      
      #ind_c indices that were not yet updated in iteration c for p
      ind_c <- is.na(p_star)
      p_star[ind_c] <- p_i[ind_c]
      
      for(k in 1:K){
        
        #update a1 for ki and a2 for ki
        a1_q <- p_i[k] + mu_ki
        a2_q <- 2 - p_i[k] - mu_ki
        
        rho_ki <- exp(sapply(0:1, function(z){-0.5*E_inv_sigma2(delta1 = delta1_q, delta2 = delta2_q)*E_Ct_C(z = z, i = i, p = p_star, mu = mu_beta_values, Sigma = Sigma_beta, B = B, y = y, k = k, K = K, iter = iter) + z*E_log_theta_ki(i = i, k = k, a1 = a1_q, a2 = a2_q) + (1-z)*E_log_theta_ki_c(i = i, k = k, a1 = a1_q, a2 = a2_q)}))
        
        p_ki_1 <- exp((-0.5*E_inv_sigma2(delta1 = delta1_q, delta2 = delta2_q)*E_Ct_C(z = 1, i = i, p = p_star, mu = mu_beta_values, Sigma = Sigma_beta, B = B, y = y, k = k, K = K, iter = iter) + E_log_theta_ki(i = i, k = k, a1 = a1_q, a2 = a2_q)))/sum(rho_ki)
        
        p_star[k] <- p_ki_1
        
        p_q <- c(p_q, p_ki_1)
        a1_q_vector <- c(a1_q_vector, a1_q)
        a2_q_vector <- c(a2_q_vector, a2_q)
        
      }
      
      p_values[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))] <- p_q
      a1_values[iter,grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))] <- a1_q_vector
      a2_values[iter,grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))] <- a2_q_vector
      
    }
    
    
  }
  
  return(list(mu_beta = mu_beta_values[maxIter,], Sigma_beta = Sigma_beta[[maxIter]], p = p_values[maxIter, ], lambda1 = lambda1_q,
              lambda2 = lambda2_values[maxIter], delta1 = delta1_q, delta2 = delta2_values[maxIter], a1 = a1_values[maxIter,], a2 = a2_values[maxIter,]))
  
}

