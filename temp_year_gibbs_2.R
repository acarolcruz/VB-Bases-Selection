# Gibbs for weather data (temperature)
source("gibbs_basis.R")
library(fda)
library(invgamma)
library(dqrng)
library(mvnfast)


data <- CanadianWeather
cities = c('Montreal','Quebec', 'Arvida', 'Bagottville', 'Sherbrooke')
y = lapply(2:5, function(x){as.numeric(data$dailyAv[1:365,cities[x],1])}) #temperature
Xt <- 1:length(y[[1]])
sds = lapply(1:4, function(x){sd(y[[x]])})
y_new <- lapply(1:4, function(x){y[[x]]/sds[[x]]})



# K = 30 ----

dir.create(paste0("Temperature_year_gibbs/K30","/chain", 2), showWarnings = FALSE)

K <- 30
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt), norder = 4, nbasis = K)
B <- lapply(1:4, function(x){getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)}) #basis

# chain 2
set.seed(1268)
z_initial <- rbinom(K, 1, 0.7)

initial_values2 <- list(tau2 = 5, sigma2 = 5, theta = rep(3/4,4*K), beta = rep(1,4*K), z = rep(z_initial, 4))

tt <- Sys.time()

out <- gibbs_basis(y = y_new, B = B, m = 4, mu_ki = 0.5, lambda1 = 0, lambda2 = 0, delta1 = 0, delta2 = 0, maxIter = 10000, K = 30, initial_values = initial_values2)

tt <- difftime(Sys.time(), tt, units = 'mins')

save(out, file = paste0("Temperature_year_gibbs/K30","/chain2/", "Gibbsresults.RData"))
write(tt, file = paste0("Temperature_year_gibbs/K30", "/chain2/", "tempo.txt"))
