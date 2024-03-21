# Simulating data (same number of time points for all curves, same coefficients and the basis used for all curves)

library(fda)
library(dqrng)
library(mvnfast)

sim_data_corr <- function(m = 5, times = 100, ordem = 4, K = 10, coef = c(-2,0,1.5,1.5,0,-1,-0.5,-1,0,0), seed = 1234, sigma = 0.1, w = 0.7){
  
  # generating basis functions (same for all curves)
  Xt = seq(0, 1, length = times)
  basisBspline_Simulated_Data = create.bspline.basis(range(Xt), norder = ordem, nbasis = K)
  B_simulated_data = getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)
  
  
  # Creating and computing correlation matrix based on a Gaussian Process
  calCov <- function(t, s, sigma=1, w=1) {
    Psi <- matrix(rep(0, length(t)*length(s)), nrow=length(t))
    for (i in 1:nrow(Psi)) {
      for (j in 1:ncol(Psi)) {
        Psi[i,j] <- (sigma^2)*exp(-2*abs(t[i]-s[j])*(1/w))
        #Psi[i,j] <- (sigma^2)*exp(-2*w*abs(t[i]-s[j]))
      }
    }
    return(Psi)
  }
  

  # Calculate the covariance matrix
  Cov_i <- calCov(Xt, Xt, sigma = sigma, w=w)
  
  # Generate a number of functions from the process
  n.samples <- m
  et <- matrix(rep(0,length(Xt)*n.samples), ncol=n.samples)
  for (i in 1:n.samples) {
    set.seed(seed+i)
    # Each column represents a sample from a multivariate normal distribution
    # with zero mean and covariance Cov_i
    et[,i] <- mvrnorm(n=1, mu=rep(0, length(Xt)), Sigma=Cov_i)
  }
  
  
  
  # generating data for each curve
  set.seed(seed)
  argvals <- lapply(1:m, function(s) Xt)
  y <- lapply(1:m, function(s) as.numeric(B_simulated_data %*% coef) + et[,s])
  B <- lapply(argvals, function(s) getbasismatrix(Xt, create.bspline.basis(range(Xt), norder = ordem, nbasis = K), nderiv = 0))
  

  
  
  return(list(y = y, B = B, Xt = Xt))
  
}


# test: sim_data_corr(m = m, times = times, ordem = ordem, K = K, coef = coef, seed = seed + sim, sigma = sigma, w = w)
