# Simulated Data
Xt = seq(0,1,length=100)
basisBspline_Simulated_Data = create.bspline.basis(range(Xt),norder=4,nbasis=10)
B_simulated_data = getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv=0)
coef = c(-2,0,1.5,1.5,0,-1,-0.5,-1,0,0)

z = rep(0, 10)
z[which(coef != 0)] <- 1

set.seed(1234)
argvals = lapply(1:5,function(s) Xt)
raw_data =lapply(1:5,function(s) as.numeric(B_simulated_data%*%coef)+rnorm(length(argvals[[s]]), 0, 0.1))
y = raw_data
K =10

B = lapply(argvals,function(s) getbasismatrix(Xt,create.bspline.basis(range(Xt),norder=4,nbasis=K), nderiv=0))


#set.seed(1234)
#y = lapply(1:5,function(s) rowSums(t(coef*z_list[[s]]*t(B_simulated_data))) + rnorm(length(argvals[[s]]), 0, 0.1))



# Latent variable - initial value of the chain
z_vector = rep(z, 5)


# 

sum(sapply(1:5, function(x){(y[[x]] - mean(y[[x]]))^2}))/(100/2)

