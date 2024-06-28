# VEM 

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


vb_bs_corr <- function(y, B, m = 5, mu_ki = 1/2, lambda_1 = 10^(-10), lambda_2 = 10^(-10), delta_1 = 10^(-10), delta_2 = 10^(-10), maxIter = 1000, K = 10, initial_values, convergence_threshold = 0.01, Xt = seq(0, 1, length = times)){
  
  converged <- FALSE
  
  w_values <- c()
  w_values[1] <- initial_values$w
  
  w <- initial_values$w

  # Calculate the covariance matrix for initial value of w
  psi <- calPsi(Xt, Xt, w = initial_values$w)
  
  K <- ncol(B[[1]]) #for fourier bases, there is another columns for the intercept
  
  mu_beta_values <- matrix(NA, maxIter, m*K) 
  colnames(mu_beta_values) <- as.vector(outer(paste0("beta_",1:K,"_"),1:m,paste0))
  Sigma_beta <- array(NA, dim = c(K, K, m))
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
  
  elbo_prev <- 0
  ELBO_values <- -Inf
  iter <- 1

  while(converged == FALSE & iter < maxIter){
    
    iter <- iter + 1
    
    for(i in 1:m){
      
      p_i <- p_values[iter - 1, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))]
      
      EG_i <- p_i*t(B[[i]])
      EGGt <- EG_i%*%solve(psi)%*%t(EG_i)
      EyG <- t(y[[i]])%*%solve(psi)%*%t(EG_i)
      
      V_i <- solve(diag(1, K)*E_inv_tau2(lambda1 = lambda1_q, lambda2 = lambda2_values[iter - 1]) + EGGt)
      
      #Update Sigma matrix for beta_i
      Sigma_bi_q <- (1/E_inv_sigma2(delta1 = delta1_q, delta2 = delta2_values[iter - 1]))*V_i
      
      #Update mean of the beta_i 
      mu_bi_q <- (E_inv_sigma2(delta1 = delta1_q, delta2 = delta2_values[iter - 1])*EyG)%*%Sigma_bi_q
      
      mu_beta_values[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))] <- mu_bi_q
      Sigma_beta[,,i] <- Sigma_bi_q
      
    }
    
    #update delta2 and lambda2
    # sigma2 posterior parameter
    delta2_q <- (sum(sapply(1:m, function(i){E_square_gi(B = B, i = i, y = y, mu = mu_beta_values, Sigma = Sigma_beta, p = p_values, iter = iter, psi = psi)})) + sum(sapply(1:m, function(i){E_sumk_beta_i_corr(i = i, mu = mu_beta_values, Sigma = Sigma_beta,  iter = iter)}))*E_inv_tau2(lambda1 = lambda1_q,lambda2 = lambda2_values[iter - 1]) + 2*delta_2)/2
    
    delta2_values[iter] <- delta2_q
    
    # tau2 posterior parameter
    lambda2_q <- (sum(sapply(1:m, function(i){E_sumk_beta_i_corr(i = i, mu = mu_beta_values, Sigma = Sigma_beta,  iter = iter)}))*E_inv_sigma2(delta1 = delta1_q, delta2 = delta2_q) + 2*lambda_2)/2
    
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
        
        #a1_ki_q <- p_star[k] + mu_ki
        #a2_ki_q <- 1 - p_i[k] + mu_ki
        #a2_ki_q <- 2 - p_star[k] - mu_ki
        #+ log(2*pi))-0.5*log(det(psi))
        log_rho_ki <- sapply(0:1, function(z){(-ni[i]/2)*E_log_sigma2(delta1 = delta1_q, delta2 = delta2_q)-0.5*E_inv_sigma2(delta1 = delta1_q, delta2 = delta2_q)*E_square_beta_i(z = z, i = i, p = p_star, mu = mu_beta_values, Sigma = Sigma_beta, B = B, y = y, k = k, K = K, iter = iter, psi = psi) + z*E_log_theta_ki(a1_ki = a1_ki_q, a2_ki = a2_ki_q) - z*E_log_theta_ki_c(a1_ki = a1_ki_q, a2_ki = a2_ki_q) + E_log_theta_ki_c(a1_ki = a1_ki_q, a2_ki = a2_ki_q)})
        
        if(sum(exp(log_rho_ki)) == 0){
          cat("sum pki = 0", "iter:", iter, "\n")
          p_ki_1 <- c(0,1)[which.max(log_rho_ki)]
        } else {
          p_ki_1 <- exp(log_rho_ki[2])/sum(exp(log_rho_ki))
          
          #p_ki_1 <- exp(log_rho_ki[2])/(1 + exp(log_rho_ki[2]))
          
          #exp(((-ni[i]/2)*E_log_sigma2(delta1 = delta1_q, delta2 = delta2_q)-0.5*E_inv_sigma2(delta1 = delta1_q, delta2 = delta2_q)*E_square_beta_i(z = 1, i = i, p = p_star, mu = mu_beta_values, Sigma = Sigma_beta, B = B, y = y, k = k, K = K, iter = iter, psi = psi) + E_log_theta_ki(a1_ki = a1_ki_q, a2_ki = a2_ki_q)))/sum(exp(log_rho_ki))
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
    
    #res <- optimize(f = elbo_omega, y = y, Xt = Xt, B = B, ni = ni, m = m, K = K, iter = iter, delta_1 = delta_1, delta_2 = delta_2, lambda_1 = lambda_1, lambda_2 = lambda_2, delta1_q = delta1_q, delta2_values = delta2_values, mu_beta_values = mu_beta_values, lambda1_q = lambda1_q, lambda2_values = lambda2_values, a1_values = a1_values, a2_values = a2_values, Sigma_beta = Sigma_beta, p_values = p_values, mu_ki = mu_ki, interval = c(0, 1e10), maximum = TRUE)
    
    #w_c <- res$maximum
    
    
    #w_c <- uniroot(dev_elbo,  m=m, delta1_q = delta1_q, delta2_q = delta2_q, Xt = Xt, iter=iter, y=y, B=B, p_values=p_values, mu_beta_values=mu_beta_values, Sigma_beta= Sigma_beta, interval = c(0, 1e10), tol = convergence_threshold)$root
    
    w_c <- optim(initial_values$w, elbo_omega, dev_elbo, y = y, Xt = Xt, B = B, ni = ni, m = m, K = K, iter = iter, delta_1 = delta_1, delta_2 = delta_2, lambda_1 = lambda_1, lambda_2 = lambda_2, delta1_q = delta1_q, delta2_values = delta2_values, mu_beta_values = mu_beta_values, lambda1_q = lambda1_q, lambda2_values = lambda2_values, a1_values = a1_values, a2_values = a2_values, Sigma_beta = Sigma_beta, p_values = p_values, mu_ki = mu_ki, control = list(fnscale = -1), method = "L-BFGS-B", lower = 1e-10, upper = Inf)$par
    #1e-1
    w_values <- c(w_values, w_c)      
    
    w_prev <- w
    w <- w_c
  
    # Calculate the covariance matrix for initial value of w
    psi <- calPsi(Xt, Xt, w = w)
    
    #elbo_c <- res$objective
    elbo_c <- elbo_corr(y = y, B = B, ni = ni, m = m, K = K, iter = iter, delta_1 = delta_1, delta_2 = delta_2, lambda_1 = lambda_1, lambda_2 = lambda_2, delta1_q = delta1_q, delta2_values = delta2_values, mu_beta_values = mu_beta_values, lambda1_q = lambda1_q, lambda2_values = lambda2_values, a1_values = a1_values, a2_values = a2_values, Sigma_beta = Sigma_beta, p_values = p_values, mu_ki = mu_ki, psi = psi)
    
    ELBO_values <- c(ELBO_values, elbo_c)
    
    converged <- check_convergence(elbo_c, elbo_prev, convergence_threshold)

    elbo_prev <- elbo_c
    
  }  
  
  if(converged == TRUE){
    return(list(data = y, B = B, mu_beta = mu_beta_values[iter,], Sigma_beta = Sigma_beta, p = p_values[iter, ], lambda1 = lambda1_q,
                lambda2 = lambda2_values[iter], delta1 = delta1_q, delta2 = delta2_values[iter], a1 = a1_values[iter,], a2 = a2_values[iter,],
                conv_elbo = elbo_c, n_iterations = iter, elbo_values = ELBO_values, w = w))
  } else{
    cat("The algorithm did not converge and the max iteration has been achieved")
    #iter_best <- which.max(ELBO_values)
    #return(list(mu_beta = mu_beta_values[maxIter,], Sigma_beta = Sigma_beta, p = p_values[maxIter, ],
    #            lambda1 = lambda1_q, lambda2 = lambda2_values[maxIter], delta1 = delta1_q, delta2 = delta2_values[maxIter],
    #            a1 = a1_values[maxIter,], a2 = a2_values[maxIter,],
    #            conv_elbo = ELBO_values[maxIter], n_iterations = iter, elbo_values = ELBO_values))}
    
    #return(coef_star_f)
    return(list(data = y, B = B))
  }
}

