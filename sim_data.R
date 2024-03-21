# Simulating data (same number of time points for all curves, same coefficients and the basis used for all curves)

library(fda)
library(dqrng)
library(mvnfast)

sim_data <- function(m = 5, times = 100, ordem = 4, K = 10, coef = c(-2,0,1.5,1.5,0,-1,-0.5,-1,0,0), seed = 1234, sigma = 0.1){
  
  # generating basis functions (same for all curves)
  Xt = seq(0, 1, length = times)
  basisBspline_Simulated_Data = create.bspline.basis(range(Xt), norder = ordem, nbasis = K)
  B_simulated_data = getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)
  
  # generating data for each curve
  set.seed(seed)
  argvals <- lapply(1:m, function(s) Xt)
  y <- lapply(1:m, function(s) as.numeric(B_simulated_data %*% coef) + rnorm(length(argvals[[s]]), 0, sigma))
  B <- lapply(argvals, function(s) getbasismatrix(Xt, create.bspline.basis(range(Xt), norder = ordem, nbasis = K), nderiv = 0))
  
  # generating true value for the z (basis selection) (same for all curves)
  #z = rep(0, 10)
  #z[which(coef != 0)] <- 1
  

  return(list(y = y, B = B))
  
}


# test: sim_data()
