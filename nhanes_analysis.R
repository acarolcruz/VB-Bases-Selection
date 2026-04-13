files <- list.files("R/", pattern = "\\.R$", full.names = TRUE)
invisible(lapply(files, source))

sds = lapply(1:8, function(x){sd(y_raw[[x]])})
y_new <- lapply(1:8, function(x){y_raw[[x]]/sds[[x]]})


library(splines)
regression_splines <- lm(y_new[[1]] ~ bs(Xt_5min,50))
summary(regression_splines)

0.09908^2


fit <- vem_fit(y = y_raw, Xt = Xt_5min, K = c(10, 20, 30), selection_metric = "mean",
               center = FALSE, scale = TRUE, delta_2 = 9*0.1,
               initial_values_fn = function(K, m) list(p = rep(1, m*K), delta2 = 4, lambda2 = 1000, w = 100))

run_comparisons <- function(y_input, y_orig, fit_obj, Xt) {
  m <- length(y_input)
  yhat_VEM <- yhat_RS <- yhat_SS <- yhat_LASSO <- yhat_BL <- vector("list", m)
  R2_VEM <- R2_RS <- R2_SS <- R2_LASSO <- R2_BL <- numeric(m)
  beta_LASSO_list <- beta_BL_list <- zhat_list <- vector("list", m)

  # Tracks global index across different K values
  current_start <- 1

  for (i in seq_len(m)) {
    nm      <- names(y_input)[i]
    y_std <- y_input[[i]]
    K       <- fit_obj$selected_K[i]
    s_sd <- as.numeric(fit$scaling_params$sds[i])

    basis_obj <- create.bspline.basis(range(Xt), norder = 4, nbasis = K)
    B_i       <- getbasismatrix(Xt, basis_obj, nderiv = 0)

    N     <- length(y_std)
    TSS   <- sum((y_std - mean(y_std))^2)

    idx_range <- current_start:(current_start + K - 1)
    zhat      <- ifelse(fit_obj$model$prob[idx_range] > 0.5, 1, 0)
    beta_vem  <- fit_obj$model$mu_beta[idx_range]

    current_start <- current_start + K # Increment for next curve

    yhat_v <- as.numeric(B_i %*% (beta_vem * zhat))
    RSS_v  <- sum((y_std - yhat_v)^2)
    df_v   <- max(1, sum(zhat))
    R2_VEM[i] <- 1 - (RSS_v / (N - df_v)) / (TSS / (N - 1))
    yhat_VEM[[i]] <- yhat_v*s_sd
    zhat_list[[i]] <- zhat

    rs      <- lm(y_std ~ B_i - 1)
    yhat_r  <- as.numeric(fitted(rs))
    RSS_r   <- sum((y_std - yhat_r)^2)
    R2_RS[i] <- 1 - (RSS_r / (N - K)) / (TSS / (N - 1))
    yhat_RS[[i]] <- yhat_r*s_sd

    loglam  <- seq(1, 9, 0.25)
    gcvs    <- vapply(loglam, function(ll)
      sum(smooth.basis(Xt, y_std, fdPar(basis_obj, 2, 10^ll))$gcv), numeric(1))
    lopt    <- 10^loglam[which.min(gcvs)]
    ss_fit  <- smooth.basis(Xt, y_std, fdPar(basis_obj, 2, lopt))
    yhat_s  <- as.numeric(B_i %*% ss_fit$fd$coefs)
    RSS_s   <- sum((y_std - yhat_s)^2)
    R2_SS[i]   <- 1 - (RSS_s / (N - ss_fit$df)) / (TSS / (N - 1))
    yhat_SS[[i]] <- yhat_s*s_sd

    reg_par <- seq(0.001, 1, length = 9)
    R2_l    <- vapply(reg_par, function(la) {
      fl  <- glmnet(B_i, y_std, alpha = 1, lambda = la / (2 * N), standardize = FALSE, intercept = FALSE)
      bl  <- as.numeric(coef(fl)[-1])
      yhl <- as.numeric(B_i %*% bl)
      1 - (sum((y_std - yhl)^2) / (N - max(1, sum(abs(bl) > 0)))) / (TSS / (N - 1))
    }, numeric(1))
    lop  <- which.max(R2_l)
    fl_f <- glmnet(B_i, y_std, alpha = 1, lambda = reg_par[lop] / (2 * N), standardize = FALSE, intercept = FALSE)
    bl_f <- as.numeric(coef(fl_f)[-1])
    yhat_la <- as.numeric(B_i %*% bl_f)
    R2_LASSO[i]    <- R2_l[lop]
    yhat_LASSO[[i]] <- yhat_la*s_sd
    beta_LASSO_list[[i]] <- bl_f

    BL1 <- blasso(B_i, y_std, T = 10000, thin = 1, beta = rep(-1, K), icept = TRUE, normalize = FALSE, verb = 0)
    BL2 <- blasso(B_i, y_std, T = 10000, thin = 1, beta = rep(1, K), icept = TRUE, normalize = FALSE, verb = 0)
    keep   <- seq(1, 5000, by = 50)
    df_bl  <- rbind(as.data.frame(cbind(beta0 = BL1$mu, BL1$beta))[5001:10000, ][keep, ],
                    as.data.frame(cbind(beta0 = BL2$mu, BL2$beta))[5001:10000, ][keep, ])
    b_bl   <- colMeans(df_bl)
    yhat_b <- as.numeric(cbind(1, B_i) %*% as.matrix(b_bl))
    RSS_b  <- sum((y_std - yhat_b)^2)
    R2_BL[i] <- 1 - (RSS_b / (N - max(1, sum(abs(b_bl[-1]) > 0)))) / (TSS / (N - 1))
    yhat_BL[[i]]  <- yhat_b*s_sd
    beta_BL_list[[i]] <- b_bl[-1]
  }

  list(yhats = list(VEM = yhat_VEM, RS = yhat_RS, SS = yhat_SS, LASSO = yhat_LASSO, BLASSO = yhat_BL),
       r2s   = data.frame(VEM = R2_VEM, RS = R2_RS, SS = R2_SS, LASSO = R2_LASSO, BLASSO = R2_BL),
       betas = list(LASSO = beta_LASSO_list, BL = beta_BL_list),
       zhats = zhat_list)
}



active_table <- function(fit_obj, res_obj, y_names) {
  start <- 1
  do.call(rbind, lapply(seq_along(y_names), function(i) {
    K <- fit_obj$selected_K[i]
    idx <- start:(start + K - 1)
    start <<- start + K
    data.frame(Curve = y_names[i], K = K, VEM_PIP05 = sum(fit_obj$model$prob[idx] > 0.5),
               LASSO_nonzero = sum(abs(res_obj$betas$LASSO[[i]]) > 0),
               BayesLASSO_nonzero = sum(abs(res_obj$betas$BL[[i]]) > 0), stringsAsFactors = FALSE)
  }))
}


res <- run_comparisons(y_new, y_raw, fit, Xt_5min)
act <- active_table(fit, res, names(y_raw))

save(y_raw, y_new,
  fit, res, act,
  Xt_5min,
  file = "NHANES_outputs.RData"
)
