########################################
# continuous outcome
# Bias-corrected Variance
# with simple exchangeable correlation
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

SV <- function(y, X, id, beta, tau2, sigma2, p=1){
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


