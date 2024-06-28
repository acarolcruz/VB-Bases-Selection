# Simulating data (same number of time points for all curves, same coefficients and the basis used for all curves)

library(fda)
library(dqrng)
library(mvnfast)

sim_data_corr <- function(m = 5, times = 100, ordem = 4, K = 10, coef = c(-2,0,1.5,1.5,0,-1,-0.5,-1,0,0), seed = 1234, sigma = 0.1, w = 0.7, basis_type = "B-splines"){
  
  # generating basis functions (same for all curves)
  #Xt = seq(0, 1, length = times)
  if(basis_type == "B-splines"){
    Xt <- seq(0, 1, length = times)
    #Xt = seq(0, 100, length = times)
    basis_data <- create.bspline.basis(range(Xt), norder = ordem, nbasis = K)
    B_simulated_data <- getbasismatrix(Xt, basis_data, nderiv = 0)
  } else if(basis_type == "Mixed"){
     # Xt <- seq(0, 1, length = times)
     # basis_data_fourier <- create.fourier.basis(range(Xt), nbasis = floor(K/2))
     # B_simulated_data_fourier <- getbasismatrix(Xt, basis_data_fourier, nderiv = 0)
     
     Xt <- seq(0, 1, length = times)
     #basis_data_fourier <- create.fourier.basis(range(Xt), dropind = 1, nbasis = round(K/2)) 
     basis_data_fourier <- create.bspline.basis(range(Xt), norder = ordem, nbasis = round(K/2))
     B_simulated_data_fourier <- getbasismatrix(Xt, basis_data_fourier, nderiv = 0)
      
     Xt <- seq(0, 1, length = times)
     basis_data_bs <- create.bspline.basis(range(Xt), norder = ordem, nbasis = round(K/2))
     B_simulated_data_bs <- getbasismatrix(Xt, basis_data_bs, nderiv = 0)
      
     B_simulated_data <- cbind(B_simulated_data_fourier, B_simulated_data_bs)
     colnames(B_simulated_data) <- as.vector(paste0("B_",1:K))
     
    #basis_data_fourier <- create.fourier.basis(range(Xt), nbasis = 7)
    #B_simulated_data_fourier <- getbasismatrix(Xt, basis_data_fourier, nderiv = 0)
    
    #basis_data_bs <- create.bspline.basis(range(Xt), norder = ordem, nbasis = 9)
    #B_simulated_data_bs <- getbasismatrix(Xt, basis_data_bs, nderiv = 0)
    
    #B_simulated_data <- list(B_simulated_data_fourier, B_simulated_data_bs)
    
  } else if (basis_type == "inverse"){
    
    Xt = seq(0, 2*pi, length = times)
    basis_data <- create.bspline.basis(range(Xt), nbasis = K)
    B_simulated_data <- getbasismatrix(Xt, basis_data, nderiv = 0)
    
  } else{
        Xt = seq(0, 2*pi, length = times)
        basis_data <- create.fourier.basis(range(Xt), nbasis = K, dropind=1)
        B_simulated_data <- getbasismatrix(Xt, basis_data, nderiv = 0)
  }
      
      
    
  
  
  
  # Calculate the covariance matrix
  Cov_i <- calCov(Xt, Xt, sigma = sigma, w=w)
  
  # Generate a number of functions from the process
  n.samples <- m
  et <- matrix(rep(0,length(Xt)*n.samples), ncol=n.samples)
  set.seed(seed)
  for (i in 1:n.samples) {
    # Each column represents a sample from a multivariate normal distribution
    # with zero mean and covariance Cov_i
    et[,i] <- mvrnorm(n=1, mu=rep(0, length(Xt)), Sigma=Cov_i)
  }
  
  
  
  # generating data for each curve
  #set.seed(seed)
  argvals <- lapply(1:m, function(s) Xt)
  if(basis_type == "B-splines"){
    y <- lapply(1:m, function(s) as.numeric(B_simulated_data %*% coef) + et[,s])
    B <- lapply(argvals, function(s) B_simulated_data)
  } else if(basis_type == "Mixed"){
    #y_f <- lapply(1:2, function(s) as.numeric(B_simulated_data[[1]] %*% coef[[1]]) + et[,s])
    #y_bs <- lapply(1:3, function(s) as.numeric(B_simulated_data[[2]] %*% coef[[2]]) + et[,s])
    
    #y <- c(y_f, y_bs)
    #B <- lapply(argvals, function(s) cbind(B_simulated_data_fourier, B_simulated_data_bs))
    
    y <- lapply(1:m, function(s) as.numeric(B_simulated_data %*% coef) + et[,s])
    B <- lapply(argvals, function(s) B_simulated_data)
  } else if(basis_type == "inverse"){
    y <- lapply(1:m, function(s) cos(Xt) + sin(2*Xt) + et[,s]) # same as Pedro's
    
    Xt <- seq(0, 2*pi, length = times)
    basis_data <- create.bspline.basis(range(Xt), norder = ordem, nbasis = K)
    B_simulated_data <- getbasismatrix(Xt, basis_data, nderiv = 0)
    
    B <- lapply(argvals, function(s) B_simulated_data)
  } else{
    y <- lapply(1:m, function(s) cos(Xt) + sin(2*Xt) + et[,s]) # same as Pedro's
    B <- lapply(argvals, function(s) B_simulated_data)
  }

  
  
  return(list(y = y, B = B, Xt = Xt))
  
}


# test: sim_data_corr(m = m, times = times, ordem = ordem, K = K, coef = coef, seed = seed + sim, sigma = sigma, w = w)
