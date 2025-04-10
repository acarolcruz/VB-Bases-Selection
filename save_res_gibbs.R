# Run this code after checking convergence and merging chain results
output_plot = "Simulation1_gibbs_miss_paper"
folder <- "Simulation1_gibbs_miss_paper"

betas <- matrix(NA, nrow = 100, ncol = 50)
z <- matrix(NA, nrow = 100, ncol = 50)

for(sim in 1:100){
    betas[sim,] <- scan(paste0(folder, "/Sim", sim, "betas.txt"))
    z[sim,] <- scan(paste0(folder,"/Sim", sim, "Z.txt"))
}

tempo <- matrix(NA, nrow = 2, ncol = 100)
for(chain in 1:2){
  tempo[chain,] <- scan(paste0(folder,"/chain", chain, "/tempo.txt"))
}

rowSums(tempo)
colSums(tempo)

betas_global <- matrix(NA, 100, 10)
z_global <- matrix(NA, 100, 10)
for(k in 1:10){
  betas_global[,k] <- apply(betas[,grep(paste0("^beta_", k, "_[0-9]+$"),
                                             colnames(res1$pars_mat))], 1, mean)
  z_global[, k] <- apply(z[,grep(paste0("^beta_", k, "_[0-9]+$"),
                                 colnames(res1$pars_mat))], 1, Mode)
} 


#bases = rep(1:10, 5)
#coef = rep(coef, 5)
eta <- matrix(NA, 100, 10)
#p <- matrix(NA, length(nsim), 10)

etakis <- (betas*z)
for(k in 1:10){
  eta[,k] <- apply(etakis[,grep(paste0("^beta_", k, "_[0-9]+$"),
                                colnames(mu_beta_values))], 1, mean)
} 



# plot with same format as pedro

boxplot(t(t(betas_global)*t(z_global)), horizontal = TRUE, yaxt="n")

axis(2, 1:10, labels=c(expression(hat(xi)[1]),
                       expression(hat(xi)[2]),
                       expression(hat(xi)[3]),
                       expression(hat(xi)[4]),
                       expression(hat(xi)[5]),
                       expression(hat(xi)[6]),
                       expression(hat(xi)[7]),
                       expression(hat(xi)[8]), 
                       expression(hat(xi)[9]),
                       expression(hat(xi)[10])), cex.axis=1.2)
abline(v = coef, lty = 2)
#c( -2.9, -0.9, -0.6, -0.5,  0.0, 0.7,  1.6,  2.4)
#c(-10.4,  -3.3,  -2.0,  -1.8,   0.0 , 2.4,   5.5,   8.4)

{
pdf(file = sprintf("%s_%d.pdf", output_plot, 1), width = 6.2, height = 7)
par(mfrow=c(1,1), mar = c(5, 4.1, 2, 0.1), cex = 1) 
boxplot(t(t(betas_global)*t(z_global)), horizontal = TRUE, yaxt="n")#, ylim = c(-2, 2) for fourier
#t(t(mu_betas)*t(z_star)
axis(2, 1:10, labels=c(expression(hat(xi)[1]),
                       expression(hat(xi)[2]),
                       expression(hat(xi)[3]),
                       expression(hat(xi)[4]),
                       expression(hat(xi)[5]),
                       expression(hat(xi)[6]),
                       expression(hat(xi)[7]),
                       expression(hat(xi)[8]), 
                       expression(hat(xi)[9]),
                       expression(hat(xi)[10])), cex.axis=1.2)
abline(v = sort(coef), lty = 2)
#c( -2.9, -0.9, -0.6, -0.5,  0.0, 0.7,  1.6,  2.4)
#c(-10.4,  -3.3,  -2.0,  -1.8,   0.0 , 2.4,   5.5,   8.4)
dev.off()
}


Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


# mean curve and credible interval 
beta_final <- apply(betas_global, 2, mean)
mode_z_gibbs <- apply(z_global, 2, Mode)





Xt <- seq(0, 1, length = 100)
basisBspline_Simulated_Data <- create.bspline.basis(range(Xt), norder = 4, nbasis = 10)
B <- getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv = 0)


# table with estimates

data.frame(Coefficients = 1:10, True = coef, Mean = as.vector(round(beta_final*mode_z_gibbs, 4)), SD = round(apply(t(t(betas_global)*t(z_global)), 2, sd),4))


# credible band for MCMC
gi_t_MCMC <- t(t(betas_global)*t(z_global))%*%t(B) # 50x100


ETI_MCMC <- apply(gi_t_MCMC, 2, function(x){quantile(x, c(0.025, 0.975))})


sim = 20

load(paste0(paste0(folder,"/chain1/Sim", sim, "IterGibbs_sample.RData")))
out_chain1 <- res1
load(paste0(paste0(folder,"/chain2/Sim", sim, "IterGibbs_sample.RData")))
out_chain2 <- res1

betas_sample <- rbind(out_chain1$pars_mat[,1:50], out_chain2$pars_mat[,1:50])
zs_sample <- rbind(out_chain1$z_mat, out_chain2$z_mat)

mode_z <- apply(zs_sample, 2, Mode)
mean_betas <- apply(betas_sample, 2, mean)

eta_plot <- matrix(NA, 200, 10)
for(k in 1:10){
  eta_plot[,k] <- apply((betas_sample*zs_sample)[,grep(paste0("^beta_", k, "_[0-9]+$"),
                                        colnames(betas_sample))], 1, mean)
}

apply(eta_plot, 2, mean)


curve_f <- list()
for(s in 1:200){
  curve_f[[s]] <- tapply((betas_sample*zs_sample)[s,], rep(1:10,5), mean)%*%t(B)
}  

LL <- apply(do.call(rbind,curve_f),2,function(i)quantile(i,probs = c(0.025,0.975)))[1,]
UL <-  apply(do.call(rbind,curve_f),2,function(i)quantile(i,probs = c(0.025,0.975)))[2,]



{
  pdf(file = sprintf("%s_%d.pdf", output_plot, 2), width = 6.2, height = 7)
  par(mar = c(5, 4, 2, 0.1) + 0.1) 
  plot(x = seq(0, 1, length = 100), y = coef%*%t(B), lwd = 2, type = "l", ylab = expression(g[t]), xlab = expression(t), ylim = c(-3, 3), col = "blue")
 lines(x = seq(0, 1, length = 100), y = apply(eta_plot, 2, mean)%*%t(B),lwd = 2, lty = 2, col= "red")
  lines(x = seq(0, 1, length = 100), y = LL, lwd = 2, lty = 3)
  lines(x = seq(0, 1, length = 100), y = UL, lwd = 2, lty = 3)
  #polygon(c(Xt, rev(Xt)), c(LL, rev(UL)), col = "#00000060", border = NA)
  legend("topright", lwd = 2, lty = c(2, 1, 3), col = c("red", "blue", "black"), legend= c("Estimated curve", "True curve", "Credible band"))
}
dev.off()

# comparing VB and MCMC (in terms of mean curve) for this to work you need to run the code in save_res.R first
plot(x = seq(0, 1, length = 100), y = as.numeric(as.vector(beta_final*mode_z_gibbs)%*% t(B)), lwd = 2, col= "red", type = "l", ylab = expression(g[t]), xlab = expression(t))
lines(x = seq(0, 1, length = 100), y = as.numeric(as.vector(mu_beta_final*mode_z)%*% t(B)), lwd = 2, col= "blue", type = "l", ylab = expression(g[t]), xlab = expression(t))
legend("topright", lwd = 1, col = c("red", "blue"), legend= c("MCMC", "VB"))

