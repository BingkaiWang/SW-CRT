########################################
# continuous outcome
# sandwich variance
# with nested exchangeable correlation
########################################

########################################
# Input
# y: N by 1 vector of outcomes
# X: N by p design matrix (including intercept)
# id: N by 1 vector of cluster id
# beta: p by 1 mean model parameter estimates
# tau2: 1 by 1 tau2 parameter estimate
# sigma2: 1 by 1 sigma2 parameter estimate
########################################

SV_NE <- function(y, X, id, beta, tau2, sigma2, p=1, kappa2, period){
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

  B <- as.matrix(Reduce("+", map(sand_var_components, ~.[[2]])) / n)
  IF <- as.matrix(- solve(B) %*% Reduce(cbind, map(sand_var_components, ~.[[1]])))
  IF_gee <- as.matrix(- solve(B[1:length(beta), 1:length(beta)]) %*% Reduce(cbind, map(sand_var_components, ~.[[1]][1:length(beta)])))
  
  if (p==1){
    rb_se <- sd(IF[1,])*sqrt(n-1)/sqrt(n)/sqrt(n)
    rb_se_gee <- sd(IF_gee[1,])*sqrt(n-1)/sqrt(n)/sqrt(n)
    return(c(beta[1], rb_se_gee, rb_se))
  } else{
    rb_se_5 <- sqrt(sum(cov(t(IF[2:6,]))))*sqrt(n-1)/sqrt(n)/sqrt(n)/5
    rb_se <- sqrt(diag(cov(t(IF[2:6,]))))*sqrt(n-1)/sqrt(n)/sqrt(n)
    
    rb_se_gee_5 <- sqrt(sum(cov(t(IF_gee[2:6,]))))*sqrt(n-1)/sqrt(n)/sqrt(n)/5
    rb_se_gee <- sqrt(diag(cov(t(IF_gee[2:6,]))))*sqrt(n-1)/sqrt(n)/sqrt(n)
    return(c(beta[1], rb_se_gee_5, rb_se_gee, rb_se_5, rb_se))
  }
  
}


