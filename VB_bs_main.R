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


vb_bs <- function(y, B, m = 5, mu_ki = 1/2, lambda_1 = 10^(-10), lambda_2 = 10^(-10), delta_1 = 10^(-10), delta_2 = 10^(-10), maxIter = 1000, K = 10, initial_values, convergence_threshold = 0.01){
  
  coef_star_f <- matrix(NA, maxIter, m*K)
  # Creating structure to save output
  mu_beta_values <- matrix(NA, maxIter, m*K) 
  colnames(mu_beta_values) <- as.vector(outer(paste0("beta_",1:K,"_"),1:m,paste0))
  Sigma_beta <- list()
  p_values <- matrix(NA, maxIter, m*K)
  lambda2_values <- c()
  delta2_values <- c()
  a1_values <-  matrix(NA, maxIter, m*K) 
  a2_values <- matrix(NA, maxIter, m*K) 
  
  p_values[1,] <- initial_values$p
  lambda2_values[1] <- initial_values$lambda2
  #delta2_values[1] <- max(sapply(1:m, function(i){var(y[[i]] - mean(y[[i]]))}))
  delta2_values[1] <- initial_values$delta2
  
  ni <- unlist(lapply(y,length))
  
  # is not being updated through variational inference
  delta1_q <- (sum(ni) + m*K + 2*delta_1)/2
  lambda1_q <- (m*K + 2*lambda_1)/2
  
  converged <- FALSE
  elbo_prev <- 0
  ELBO_values <- -Inf
  iter <- 1
  mu_prev <- initial_values$mu_0
  coef_star_prev <- initial_values$coef_0
  
  while(converged == FALSE & iter < maxIter){
    iter <- iter + 1
    
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
      
      V_i <- solve(diag(1, K)*E_inv_tau2(lambda1 = lambda1_q, lambda2 = lambda2_values[iter - 1]) + Reduce('+',HHt))
      
      #Update Sigma matrix for beta_i
      Sigma_bi_q <- (1/E_inv_sigma2(delta1 = delta1_q, delta2 = delta2_values[iter - 1]))*V_i
      
      #Update mean of the beta_i 
      mu_bi_q <- E_inv_sigma2(delta1 = delta1_q, delta2 = delta2_values[iter - 1])*Sigma_bi_q%*%Reduce('+',Hy)
      
      mu_beta_values[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))] <- mu_bi_q
      Sigma_beta[[iter]][, grep(paste0("^beta_[0-9]+_",i,"$"),colnames(mu_beta_values))] <- Sigma_bi_q
      
    }
    
    #update delta2 and lambda2
    delta2_q <- (sum(sapply(1:m, function(i){E_square_i(B = B, i = i, y = y, mu = mu_beta_values, Sigma = Sigma_beta, p = p_values, iter = iter)})) + sum(sapply(1:m, function(i){E_sumk_beta_i(i = i, mu = mu_beta_values, Sigma = Sigma_beta,  iter = iter)}))*E_inv_tau2(lambda1 = lambda1_q,lambda2 = lambda2_values[iter - 1]) + 2*delta2_values[1])/2
    
    delta2_values[iter] <- delta2_q
    
    lambda2_q <- (sum(sapply(1:m, function(i){E_sumk_beta_i(i = i, mu = mu_beta_values, Sigma = Sigma_beta,  iter = iter)}))*E_inv_sigma2(delta1 = delta1_q, delta2 = delta2_q) + 2*lambda2_values[1])/2
    
    lambda2_values[iter] <- lambda2_q
    

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
        a1_ki_q <- p_i[k] + mu_ki
        a2_ki_q <- 2 - p_i[k] - mu_ki
        
        log_rho_ki <- sapply(0:1, function(z){(-ni[i]/2)*E_log_sigma2(delta1 = delta1_q, delta2 = delta2_q)-0.5*E_inv_sigma2(delta1 = delta1_q, delta2 = delta2_q)*E_Ct_C(z = z, i = i, p = p_star, mu = mu_beta_values, Sigma = Sigma_beta, B = B, y = y, k = k, K = K, iter = iter) + z*E_log_theta_ki(a1_ki = a1_ki_q, a2_ki = a2_ki_q) + (1-z)*E_log_theta_ki_c(a1_ki = a1_ki_q, a2_ki = a2_ki_q)})
        
        if(sum(exp(log_rho_ki)) == 0){
          cat("sum pki = 0", "iter:", iter, "\n")
          p_ki_1 <- c(0,1)[which.max(log_rho_ki)]
        } else {
          p_ki_1 <- exp(((-ni[i]/2)*E_log_sigma2(delta1 = delta1_q, delta2 = delta2_q)-0.5*E_inv_sigma2(delta1 = delta1_q, delta2 = delta2_q)*E_Ct_C(z = 1, i = i, p = p_star, mu = mu_beta_values, Sigma = Sigma_beta, B = B, y = y, k = k, K = K, iter = iter) + E_log_theta_ki(a1_ki = a1_ki_q, a2_ki = a2_ki_q)))/sum(exp(log_rho_ki))
        }
        
        p_star[k] <- p_ki_1
        
        p_q <- c(p_q, p_ki_1)
        a1_q_vector <- c(a1_q_vector, a1_ki_q)
        a2_q_vector <- c(a2_q_vector, a2_ki_q)
        
      }
      
      p_values[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))] <- p_q
      a1_values[iter,grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))] <- a1_q_vector
      a2_values[iter,grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))] <- a2_q_vector
      
    }
    
    #elbo calculation

    elbo_c <- elbo(y = y, B = B, ni = ni, m = m, K = K, iter = iter, delta_1 = delta_1, delta_2 = delta_2,
                   lambda_1 = lambda_1, lambda_2 = lambda_2,
                   delta1_q = delta1_q, delta2_values = delta2_values,
                   mu_beta_values = mu_beta_values, lambda1_q = lambda1_q, lambda2_values = lambda2_values,
                   a1_values = a1_values, a2_values = a2_values, Sigma_beta = Sigma_beta, p_values = p_values, mu_ki = mu_ki)
  
    cat(elbo_c, "\n")
    
    ELBO_values <- c(ELBO_values, elbo_c)
    #converged_1 <- ifelse(abs(elbo_c - elbo_prev) <= convergence_threshold, TRUE, FALSE)
    #converged_2 <- ifelse(sum(abs(mu_prev - mu_beta_values[iter,])) <= convergence_threshold, TRUE, FALSE)
    #converged <- ifelse(converged_1 == FALSE & converged_2 == FALSE, FALSE, TRUE)
    converged <- check_convergence(elbo_c, elbo_prev, convergence_threshold)
    #cat(sum(abs(mu_prev - mu_beta_values[iter,])), "\n")  
    #cat(mu_beta_values[iter, 3], "\n") 
    elbo_prev <- elbo_c
    #mu_prev <- mu_beta_values[iter,]
    
    # convergence of Z*betas
    
    # coef_star_c <- c()
    # for(i in 1:m){
    #   mu_i_q <- mu_beta_values[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))]
    #   p_iq <- p_values[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))] 
    #   for(k in 1:K){
    #     z_star <- ifelse(p_iq[k] > 0.5, 1, 0)
    #     coef_star_c <- c(coef_star_c, z_star*mu_i_q[k])
    #   }
    # }
    
    #cat(coef_star_c, "\n")
    #print(as.vector(coef_star_c))
    #cat(sum(abs(coef_star_prev - coef_star_c)), "\n")  
    
    #coef_star_prev <- coef_star_c
    
    
    #coef_star_f[iter, ] <- coef_star_c
  }
  
  if(iter == maxIter & converged == FALSE){cat("The algorithm did not converge and the max iteration has been achieved")}

   if(converged == TRUE){
     return(list(mu_beta = mu_beta_values[iter,], Sigma_beta = Sigma_beta[[iter]], p = p_values[iter, ], lambda1 = lambda1_q,
               lambda2 = lambda2_values[iter], delta1 = delta1_q, delta2 = delta2_values[iter], a1 = a1_values[iter,], a2 = a2_values[iter,],
               conv_elbo = elbo_c, n_iterations = iter, elbo_values = ELBO_values))
   } else{
     #iter_best <- which.max(ELBO_values)
     return(list(mu_beta = mu_beta_values[maxIter,], Sigma_beta = Sigma_beta[[maxIter]], p = p_values[maxIter, ],
                 lambda1 = lambda1_q, lambda2 = lambda2_values[maxIter], delta1 = delta1_q, delta2 = delta2_values[maxIter],
                 a1 = a1_values[maxIter,], a2 = a2_values[maxIter,],
                 conv_elbo = ELBO_values[maxIter], n_iterations = iter, elbo_values = ELBO_values))}

  return(coef_star_f)
}

