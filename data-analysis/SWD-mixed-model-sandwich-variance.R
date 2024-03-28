########################################
# Sandwich variance estimator with linear mixed models for stepped wedge designs
# assuming an exchangeable working correlation
########################################

########################################
# Input
# y: N by 1 vector of outcomes
# X: N by p design matrix (including period-specific intercept, treatment terms, and covariates)
# id: N by 1 vector of cluster id
# beta: p by 1 mean model parameter estimates
# tau2: 1 by 1 tau2 parameter estimate for variance of the cluster-random intercept
# sigma2: 1 by 1 sigma2 parameter estimate for the individual residual variance
# n: number of clusters
# n_period: number of randomization steps
# trt_model: "constant", "duration-specific", "period-specific", or "saturated"
########################################


SV_SWD_i <- function(y, X, id, beta, n, n_period, trt_model){
  sand_var_components <- map(1:n, function(i){
    y_i <- y[id == i]
    X_i <- X[id == i, ]
    m_i <- nrow(X_i)
    psi_i <- 2 * t(X_i) %*% (y_i - X_i %*% beta)
    B_i <- -2 * t(X_i) %*% X_i
    list(psi_i, B_i)
  })
  B <- Reduce("+", map(sand_var_components, ~.[[2]])) / n
  IF <-  - solve(B) %*% Reduce(cbind, map(sand_var_components, ~.[[1]]))
  
  if (trt_model == "constant"){
    rb_se <- sd(IF[n_period+1,])*sqrt(n-1)/sqrt(n)/sqrt(n)
    return(c(beta[n_period+1], rb_se))
  } else  if (trt_model == "duration-specific"){
    # trt_index <- n_period:(2*n_period-2)
    trt_index <- n_period + 1:n_period
    rb_se <- map_dbl(trt_index, ~sd(IF[.,])) *sqrt(n-1)/sqrt(n)/sqrt(n)
    rb_se_avg <- sqrt(mean(var(t(IF[trt_index,]))))*sqrt(n-1)/sqrt(n)/sqrt(n)
    return(c(beta[trt_index], mean(beta[trt_index]), rb_se, rb_se_avg))
  } else  if (trt_model == "period-specific"){
    trt_index <- n_period + 1:n_period
    rb_se <- map_dbl(trt_index, ~sd(IF[.,])) *sqrt(n-1)/sqrt(n)/sqrt(n)
    rb_se_avg <- sqrt(mean(var(t(IF[trt_index,]))))*sqrt(n-1)/sqrt(n)/sqrt(n)
    return(c(beta[trt_index], mean(beta[trt_index]), rb_se, rb_se_avg))
  } else if (trt_model == "saturated") {
    trt_index <- n_period + 1:(n_period*(n_period+1)/2)
    rb_se <- map_dbl(trt_index, ~sd(IF[.,])) *sqrt(n-1)/sqrt(n)/sqrt(n)
    rb_se_avg <- sqrt(mean(var(t(IF[trt_index,]))))*sqrt(n-1)/sqrt(n)/sqrt(n)
    return(c(beta[trt_index], mean(beta[trt_index]), rb_se, rb_se_avg))
  }
}

SV_SWD_e <- function(y, X, id, beta, tau2, sigma2, n, n_period, trt_model){
  sand_var_components <- map(1:n, function(i){
    y_i <- y[id == i]
    X_i <- X[id == i, ]
    m_i <- nrow(X_i)
    Sigma_i <- tau2 * matrix(1, nrow = m_i, ncol = m_i) + sigma2 * diag(m_i)
    Sigma_i_inv <- solve(Sigma_i)
    psi_i <- rbind(
      2 * t(X_i) %*% Sigma_i_inv %*%  (y_i - X_i %*% beta),
      - sum(diag(Sigma_i_inv)) + sum((Sigma_i_inv %*%  (y_i - X_i %*% beta))^2),
      - sum(Sigma_i_inv) + (t(y_i - X_i %*% beta) %*% rowSums(Sigma_i_inv))^2
    )
    B11 <- -2 * t(X_i) %*% Sigma_i_inv %*% X_i
    B12 <- -2 * t(X_i) %*% Sigma_i_inv %*% Sigma_i_inv %*% (y_i - X_i %*% beta)
    B13 <- -2 * t(X_i) %*% rowSums(Sigma_i_inv) %*% colSums(Sigma_i_inv) %*% (y_i - X_i %*% beta)
    B22 <- sum(diag(Sigma_i_inv %*% Sigma_i_inv)) - 
      2 * t(y_i - X_i %*% beta) %*% Sigma_i_inv %*% Sigma_i_inv %*% Sigma_i_inv %*% (y_i - X_i %*% beta)
    B23 <- sum(Sigma_i_inv %*% Sigma_i_inv) - 
      2 * t(y_i - X_i %*% beta) %*% rowSums(Sigma_i_inv) * t(y_i - X_i %*% beta) %*% rowSums(Sigma_i_inv %*% Sigma_i_inv)
    B33 <- sum(Sigma_i_inv)^2 - 2 * sum(Sigma_i_inv) * (t(y_i - X_i %*% beta) %*% rowSums(Sigma_i_inv))^2
    B_i <- rbind(
      cbind(B11, B12, B13),
      cbind(t(B12), B22, B23),
      cbind(t(B13), t(B23), B33)
    )
    list(psi_i, B_i)
  })
  
  B <- Reduce("+", map(sand_var_components, ~.[[2]])) / n
  IF <-  - solve(B) %*% Reduce(cbind, map(sand_var_components, ~.[[1]]))
  IF_gee <- - solve(B[1:length(beta), 1:length(beta)]) %*% Reduce(cbind, map(sand_var_components, ~.[[1]][1:length(beta)]))
  
  if (trt_model == "constant"){
    rb_se <- sd(IF[n_period+1,])*sqrt(n-1)/sqrt(n)/sqrt(n)
    rb_se_gee <- sd(IF_gee[n_period+1,])*sqrt(n-1)/sqrt(n)/sqrt(n)
    return(c(beta[n_period+1], rb_se, rb_se_gee))
  } else  if (trt_model == "duration-specific"){
    # trt_index <- n_period:(2*n_period-2)
    trt_index <- n_period + 1:n_period
    rb_se <- map_dbl(trt_index, ~sd(IF[.,])) *sqrt(n-1)/sqrt(n)/sqrt(n)
    rb_se_gee <- map_dbl(trt_index, ~sd(IF_gee[.,])) *sqrt(n-1)/sqrt(n)/sqrt(n)
    rb_se_avg <- sqrt(mean(var(t(IF[trt_index,]))))*sqrt(n-1)/sqrt(n)/sqrt(n)
    rb_se_gee_avg <- sqrt(mean(var(t(IF_gee[trt_index,]))))*sqrt(n-1)/sqrt(n)/sqrt(n)
    return(c(beta[trt_index], mean(beta[trt_index]), rb_se, rb_se_avg, rb_se_gee, rb_se_gee_avg))
  } else  if (trt_model == "period-specific"){
    trt_index <- n_period + 1:n_period
    rb_se <- map_dbl(trt_index, ~sd(IF[.,])) *sqrt(n-1)/sqrt(n)/sqrt(n)
    rb_se_gee <- map_dbl(trt_index, ~sd(IF_gee[.,])) *sqrt(n-1)/sqrt(n)/sqrt(n)
    rb_se_avg <- sqrt(mean(var(t(IF[trt_index,]))))*sqrt(n-1)/sqrt(n)/sqrt(n)
    rb_se_gee_avg <- sqrt(mean(var(t(IF_gee[trt_index,]))))*sqrt(n-1)/sqrt(n)/sqrt(n)
    return(c(beta[trt_index], mean(beta[trt_index]), rb_se, rb_se_avg, rb_se_gee, rb_se_gee_avg))
  } else if (trt_model == "saturated") {
    trt_index <- n_period + 1:(n_period*(n_period+1)/2)
    rb_se <- map_dbl(trt_index, ~sd(IF[.,])) *sqrt(n-1)/sqrt(n)/sqrt(n)
    rb_se_gee <- map_dbl(trt_index, ~sd(IF_gee[.,])) *sqrt(n-1)/sqrt(n)/sqrt(n)
    rb_se_avg <- sqrt(mean(var(t(IF[trt_index,]))))*sqrt(n-1)/sqrt(n)/sqrt(n)
    rb_se_gee_avg <- sqrt(mean(var(t(IF_gee[trt_index,]))))*sqrt(n-1)/sqrt(n)/sqrt(n)
    return(c(beta[trt_index], mean(beta[trt_index]), rb_se, rb_se_avg, rb_se_gee, rb_se_gee_avg))
  }
  
}



SV_SWD_ne <- function(y, X, id, beta, tau2, sigma2, kappa2, n, n_period, period, trt_model){
  sand_var_components <- map(1:n, function(i){
    y_i <- y[id == i]
    X_i <- X[id == i, ]
    period_i <- period[id == i]
    m_i <- nrow(X_i)
    n_ij <- as.data.frame(cbind(y_i, period_i)) %>% group_by(period_i) %>% summarise(N_ij = n()) %>% .$N_ij
    bdiag_mat <- bdiag(map(n_ij, ~matrix(1, nrow = ., ncol =.)))
    Sigma_i <- tau2 * matrix(1, nrow = m_i, ncol = m_i) + sigma2 * diag(m_i) + kappa2 * bdiag_mat
    Sigma_i_inv <- solve(Sigma_i)
    psi_i <- rbind(
      2 * t(X_i) %*% Sigma_i_inv %*%  (y_i- X_i %*% beta),
      - sum(diag(Sigma_i_inv)) + sum((Sigma_i_inv %*%  (y_i- X_i %*% beta))^2),
      - sum(Sigma_i_inv) + (t(y_i- X_i %*% beta) %*% rowSums(Sigma_i_inv))^2,
      - sum(diag(Sigma_i_inv %*% bdiag_mat)) + t(y_i- X_i %*% beta) %*% Sigma_i_inv %*% bdiag_mat %*% Sigma_i_inv %*% (y_i- X_i %*% beta)
    )
    B11 <- -2 * t(X_i) %*% Sigma_i_inv %*% X_i
    B12 <- -2 * t(X_i) %*% Sigma_i_inv %*% Sigma_i_inv %*% (y_i- X_i %*% beta)
    B13 <- -2 * t(X_i) %*% rowSums(Sigma_i_inv) %*% colSums(Sigma_i_inv) %*% (y_i- X_i %*% beta)
    B14 <- -2 * t(X_i) %*% Sigma_i_inv %*% bdiag_mat %*% Sigma_i_inv %*% (y_i- X_i %*% beta)
    B22 <- sum(diag(Sigma_i_inv %*% Sigma_i_inv)) - 
      2 * t(y_i- X_i %*% beta) %*% Sigma_i_inv %*% Sigma_i_inv %*% Sigma_i_inv %*% (y_i- X_i %*% beta)
    B23 <- sum(Sigma_i_inv %*% Sigma_i_inv) - 
      2 * t(y_i- X_i %*% beta) %*% rowSums(Sigma_i_inv) * t(y_i- X_i %*% beta) %*% rowSums(Sigma_i_inv %*% Sigma_i_inv)
    B24 <- sum(diag(Sigma_i_inv %*% Sigma_i_inv %*% bdiag_mat)) - 
      2 * t(y_i- X_i %*% beta) %*% Sigma_i_inv %*% bdiag_mat %*% Sigma_i_inv %*% Sigma_i_inv %*% (y_i- X_i %*% beta)
    B33 <- sum(Sigma_i_inv)^2 - 2 * sum(Sigma_i_inv) * (t(y_i- X_i %*% beta) %*% rowSums(Sigma_i_inv))^2
    B34 <- sum(Sigma_i_inv %*% bdiag_mat %*% Sigma_i_inv) -
      2 * t(y_i- X_i %*% beta) %*% Sigma_i_inv %*% bdiag_mat %*% rowSums(Sigma_i_inv) * colSums(Sigma_i_inv) %*% (y_i- X_i %*% beta)
    B44 <- sum(diag(Sigma_i_inv %*% bdiag_mat %*% Sigma_i_inv %*% bdiag_mat)) -
      2 * t(y_i- X_i %*% beta) %*% Sigma_i_inv %*% bdiag_mat %*% Sigma_i_inv %*% bdiag_mat %*% Sigma_i_inv %*% (y_i- X_i %*% beta) 
    B_i <- rbind(
      cbind(B11, B12, B13, B14),
      cbind(t(B12), B22, B23, B24),
      cbind(t(B13), t(B23), B33, B34),
      cbind(t(B14), t(B24), t(B34), B44)
    )
    list(psi_i, B_i)
  })
  
  B <- Reduce("+", map(sand_var_components, ~.[[2]])) / n
  IF <-  - solve(B) %*% Reduce(cbind, map(sand_var_components, ~.[[1]]))
  IF_gee <- - solve(B[1:length(beta), 1:length(beta)]) %*% Reduce(cbind, map(sand_var_components, ~.[[1]][1:length(beta)]))
  
  if (trt_model == "constant"){
    rb_se <- sd(IF[n_period+1,])*sqrt(n-1)/sqrt(n)/sqrt(n)
    rb_se_gee <- sd(IF_gee[n_period+1,])*sqrt(n-1)/sqrt(n)/sqrt(n)
    return(c(beta[n_period+1], rb_se, rb_se_gee))
  } else  if (trt_model == "duration-specific"){
    # trt_index <- n_period:(2*n_period-2)
    trt_index <- n_period + 1:n_period
    rb_se <- map_dbl(trt_index, ~sd(IF[.,])) *sqrt(n-1)/sqrt(n)/sqrt(n)
    rb_se_gee <- map_dbl(trt_index, ~sd(IF_gee[.,])) *sqrt(n-1)/sqrt(n)/sqrt(n)
    rb_se_avg <- sqrt(mean(var(as.matrix(t(IF[trt_index,])))))*sqrt(n-1)/sqrt(n)/sqrt(n)
    rb_se_gee_avg <- sqrt(mean(var(as.matrix(t(IF_gee[trt_index,])))))*sqrt(n-1)/sqrt(n)/sqrt(n)
    return(c(beta[trt_index], mean(beta[trt_index]), rb_se, rb_se_avg, rb_se_gee, rb_se_gee_avg))
  } else  if (trt_model == "period-specific"){
    trt_index <- n_period + 1:n_period
    rb_se <- map_dbl(trt_index, ~sd(IF[.,])) *sqrt(n-1)/sqrt(n)/sqrt(n)
    rb_se_gee <- map_dbl(trt_index, ~sd(IF_gee[.,])) *sqrt(n-1)/sqrt(n)/sqrt(n)
    rb_se_avg <- sqrt(mean(var(as.matrix(t(IF[trt_index,])))))*sqrt(n-1)/sqrt(n)/sqrt(n)
    rb_se_gee_avg <- sqrt(mean(var(as.matrix(t(IF_gee[trt_index,])))))*sqrt(n-1)/sqrt(n)/sqrt(n)
    return(c(beta[trt_index], mean(beta[trt_index]), rb_se, rb_se_avg, rb_se_gee, rb_se_gee_avg))
  } else if (trt_model == "saturated") {
    trt_index <- n_period + 1:(n_period*(n_period+1)/2)
    rb_se <- map_dbl(trt_index, ~sd(IF[.,])) *sqrt(n-1)/sqrt(n)/sqrt(n)
    rb_se_gee <- map_dbl(trt_index, ~sd(IF_gee[.,])) *sqrt(n-1)/sqrt(n)/sqrt(n)
    rb_se_avg <- sqrt(mean(var(as.matrix(t(IF[trt_index,])))))*sqrt(n-1)/sqrt(n)/sqrt(n)
    rb_se_gee_avg <- sqrt(mean(var(as.matrix(t(IF_gee[trt_index,])))))*sqrt(n-1)/sqrt(n)/sqrt(n)
    return(c(beta[trt_index], mean(beta[trt_index]), rb_se, rb_se_avg, rb_se_gee, rb_se_gee_avg))
  }
  
}
