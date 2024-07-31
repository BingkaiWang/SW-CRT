#############################################################
# Simulation for SWD (with Data Generation Process B2)

# October 2022

# INPUT
# n: Total number of clusters
# t: period
# pplv: Vector of population size in each cluster
# cpl: Lower bound of the cluster-period sizes
# cpu: Upper bound of the cluster-period sizes
# beta0: Parameter of the duration-related treatment effect
# beta1: Parameter of the covariate-related treatment effect
# mu_bx: Mean vector of binary covariates
# sd_cxc: Sd vector of cluster-level random effect(s) for continuous covariates
# sd_cxe: Sd vector of random error(s) for continuous covariates
# sdb: Sd of the cluster randomized effects
# sdep: Sd of the random errors
#############################################################
library(lme4)
library(performance)
library(Matrix)
library(tidyverse)

ATE <- function(n, t, pplv, cpl, cpu, beta0, beta1, mu_bx, sd_cxc, sd_cxe, sdb, sdep){
  
  # Source program code
  source("contGEE_BCV2.R")
  source("SV.R")
  source("SV_NE.R")
  
  # true cluster-period size
  mvt <- rep(pplv, each=t)
  mvts <- c(0, cumsum(mvt))
  # cluster id
  cluster <- rep(1:n, times=(pplv*t))
  # period id
  period <- rep(rep(1:t, times=n), times=mvt)
  
  # period effects
  bas <- 0.25
  end <- 0.27
  p <- seq(bas, end, length.out=t)
  pv <- rep(rep(p, times=n), times=mvt)
  
  # covariate(s)
  xn <- length(mu_bx) + length(sd_cxc)
  xmt <- NULL
  for (k in 1:length(mu_bx)){
    xvk <- NULL
    for (i in 1:n){ # independent binary covariates within each cluster
      xvki <- rbinom(n=pplv[i], size=1, prob=mu_bx[k])
      xvk <- c(xvk, rep(xvki, times=t))
    }
    xmt <- cbind(xmt, as.vector(xvk))
  }
  for (k in 1:length(sd_cxc)){
    xvk <- NULL
    for (i in 1:n){ # correlated continuous covariates within each cluster
      xvki <- rnorm(1, 0, sd_cxc[k]) + rnorm(pplv[i], 0, sd_cxe[k])
      xvk <- c(xvk, rep(xvki, times=t))
    }
    xmt <- cbind(xmt, as.vector(xvk))
  }
  
  # cluster randomized effects
  gamma_ba <- sdb[1]^(-2)
  gamma_bb <- 1/(sdb[1]^2)
  bvt <- rgamma(n, shape=gamma_ba, rate=gamma_bb)
  bvt <- bvt - mean(bvt)
  bv <- rep(rep(bvt, each=t), times=mvt)
  # cluster-period randomized effects
  gamma_ca <- sdb[2]^(-2)
  gamma_cb <- 1/(sdb[2]^2)
  cvt <- rgamma(n*t, shape=gamma_ca, rate=gamma_cb)
  cvt <- cvt - mean(cvt)
  cv <- rep(cvt, times=mvt)
  # random errors
  epvt <- rpois(sum(pplv), lambda=sdep^2)
  epvt <- epvt - mean(epvt)
  pplvs <- c(0, cumsum(pplv))
  epv <-NULL
  for (i in 1:n){
    epv <- c(epv, rep(epvt[(pplvs[i]+1):(pplvs[i+1])], times=t))
  }
  
  # cluster intervention effect
  cievt <- rnorm(n, 0, 0.5)
  cievt <- cievt - mean(cievt)
  ciev <- rep(rep(cievt, each=t), times=mvt)
  
  # data frame
  df <- data.frame(cbind(cluster, period, xmt))
  df$period <- as.factor(df$period)
  xname <- NULL
  for (k in 1:xn){
    xname <- c(xname, paste0("X", k))
  }
  colnames(df) <- c("CLUSTER", "PERIOD", xname)

  # potential outcomes
  df$X1_mean <- ave(df$X1, df$CLUSTER)
  df$X33_mean <- ave(df$X3^3, df$CLUSTER)
  
  xs <- period/2*3*(df$X1) + df$X2 + period/t*6*(df$X3^2) + df$X4
  
  yd <- pv + bv + cv + epv + xs
  for (j in 1:(t-1)){
    betaj <- beta0[j] + beta1[j]/2*period*(df$X1-df$X1_mean) + beta1[j]*period/t*(df$X3^3-df$X33_mean)
    yj <- pv + bv + cv + epv + xs + betaj + ciev
    yd <- cbind(yd, yj)
  }
  
  # observed cluster sizes and cluster-period sizes
  mv <- NULL  # cluster-period sizes
  idv <- NULL  # id vector of random selected observed individuals
  for (i in 1:n){
    mviv <- sample(cpl:cpu, t, replace=T)
    mv <- c(mv, mviv)
    
    a <- t*(i-1)
    for (j in 1:t){
      idviv <- sample(seq(mvts[a+j]+1, mvts[a+j+1], 1), mviv[j], replace=FALSE)
      idv <- c(idv, idviv)
    }
  }
  
  # total number of observed individuals
  N <- sum(mv)
  
  # potential outcomes and covariates of observed individuals
  dfo <- df[idv, ]
  ydo <- yd[idv, ]

  # indicator of intervention arm; i = 1 treated; i = 0 non treated
  trt_rd <- sample(rep(1:(t-1), each=n/(t-1)))
  design <- matrix(0, t-1, t)
  design[upper.tri(design)] <- 1
  
  design_eff <- matrix(0, t-1, t)
  uptri <- NULL
  for (j in 1:(t-1)){
    uptri <- c(uptri, j:1)
  }
  design_eff[upper.tri(design_eff)] <- uptri
  
  trt_seq <- eff_seq <- NULL
  for (i in 1:n){
    trt_seq <- rbind(trt_seq, rep(design[trt_rd[i], ]))
    eff_seq <- rbind(eff_seq, rep(design_eff[trt_rd[i], ]))
  }
  dfo$TRT <- rep(c(t(trt_seq)), times = mv)
  eff <- rep(c(t(eff_seq)), times = mv)
  dfo$D <- eff
  
  des_ma <- dname <- NULL
  for (j in 1:t){
    jj <- j-1
    des_ma <- cbind(des_ma, as.numeric(eff==jj))
    dname <- c(dname, paste0("D", jj))
  }
  dfonames <- colnames(dfo)
  dfo <- cbind(dfo, des_ma)
  colnames(dfo) <- c(dfonames, dname)
  dfo$D <- as.factor(dfo$D)
  
  # observed values
  y <- rowSums(ydo*des_ma)
  dfo$Y <- y
  id <- dfo$CLUSTER
  period_o <- dfo$PERIOD
  phi <- var(y)

  ###################################### simple exchangeable correlation ######################################
  # fit model 1t
  fit_lme_un_t <- lmer(Y ~ 0 + TRT + PERIOD + (1 | CLUSTER), data = dfo)
  beta <- as.numeric(summary(fit_lme_un_t)$coefficients[, 1])
  X <- model.matrix(~ 0 + TRT + PERIOD, data = dfo)
  alpha_un_t <- icc(fit_lme_un_t)$ICC_adjusted
  fit_un_t <- contGEEBCV(y=y,X=X,beta=beta,alpha=alpha_un_t,phi=phi,id=id)
  
  tau2 <- summary(fit_lme_un_t)$varcor$CLUSTER[1]
  sigma2 <- summary(fit_lme_un_t)$sigma^2
  c_un_t <- SV(y, X, id, beta, tau2, sigma2)
  
  # fit model 1d
  fit_lme_un_d <- lmer(Y ~ D + PERIOD + (1 | CLUSTER), data = dfo)
  beta <- as.numeric(summary(fit_lme_un_d)$coefficients[, 1])
  X <- model.matrix(~ D + PERIOD, data = dfo)
  alpha_un_d <- icc(fit_lme_un_d)$ICC_adjusted
  fit_un_d <- contGEEBCV(y=y,X=X,beta=beta,alpha=alpha_un_d,phi=phi,id=id)
  
  tau2 <- summary(fit_lme_un_d)$varcor$CLUSTER[1]
  sigma2 <- summary(fit_lme_un_d)$sigma^2
  c_un_d <- SV(y, X, id, beta, tau2, sigma2, 5)
  
  # fit model 2t
  fit_lme_ca_t <- lmer(Y ~ 0 + TRT + PERIOD + X1 + X3 + (1 | CLUSTER), data = dfo)
  beta <- as.numeric(summary(fit_lme_ca_t)$coefficients[, 1])
  X <- model.matrix(~ 0 + TRT + PERIOD + X1 + X3, data = dfo)
  alpha_ca_t <- icc(fit_lme_ca_t)$ICC_adjusted
  fit_ca_t <- contGEEBCV(y=y,X=X,beta=beta,alpha=alpha_ca_t,phi=phi,id=id)
  
  tau2 <- summary(fit_lme_ca_t)$varcor$CLUSTER[1]
  sigma2 <- summary(fit_lme_ca_t)$sigma^2
  c_ca_t <- SV(y, X, id, beta, tau2, sigma2)
  
  # fit model 2d
  fit_lme_ca_d <- lmer(Y ~ D + PERIOD + X1 + X3 + (1 | CLUSTER), data = dfo)
  beta <- as.numeric(summary(fit_lme_ca_d)$coefficients[, 1])
  X <- model.matrix(~ D + PERIOD + X1 + X3, data = dfo)
  alpha_ca_d <- icc(fit_lme_ca_d)$ICC_adjusted
  fit_ca_d <- contGEEBCV(y=y,X=X,beta=beta,alpha=alpha_ca_d,phi=phi,id=id)
  
  tau2 <- summary(fit_lme_ca_d)$varcor$CLUSTER[1]
  sigma2 <- summary(fit_lme_ca_d)$sigma^2
  c_ca_d <- SV(y, X, id, beta, tau2, sigma2, 5)
  
  # fit model 3t
  fit_lme_cac_t <- lmer(Y ~ 0 + TRT + PERIOD + X1 + X2 + X3 + X4 + (1 | CLUSTER), data = dfo)
  beta <- as.numeric(summary(fit_lme_cac_t)$coefficients[, 1])
  X <- model.matrix(~ 0 + TRT + PERIOD + X1 + X2 + X3 + X4, data = dfo)
  alpha_cac_t <- icc(fit_lme_cac_t)$ICC_adjusted
  fit_cac_t <- contGEEBCV(y=y,X=X,beta=beta,alpha=alpha_cac_t,phi=phi,id=id)
  
  tau2 <- summary(fit_lme_cac_t)$varcor$CLUSTER[1]
  sigma2 <- summary(fit_lme_cac_t)$sigma^2
  c_cac_t <- SV(y, X, id, beta, tau2, sigma2)
  
  # fit model 3d
  fit_lme_cac_d <- lmer(Y ~ D + PERIOD + X1 + X2 + X3 + X4 + (1 | CLUSTER), data = dfo)
  beta <- as.numeric(summary(fit_lme_cac_d)$coefficients[, 1])
  X <- model.matrix(~ D + PERIOD + X1 + X2 + X3 + X4, data = dfo)
  alpha_cac_d <- icc(fit_lme_cac_d)$ICC_adjusted
  fit_cac_d <- contGEEBCV(y=y,X=X,beta=beta,alpha=alpha_cac_d,phi=phi,id=id)
  
  tau2 <- summary(fit_lme_cac_d)$varcor$CLUSTER[1]
  sigma2 <- summary(fit_lme_cac_d)$sigma^2
  c_cac_d <- SV(y, X, id, beta, tau2, sigma2, 5)

  # estimates
  estimates <- c(summary(fit_lme_un_t)$coefficients[1, 1], sum(summary(fit_lme_un_d)$coefficients[2:6, 1])/5, summary(fit_lme_un_d)$coefficients[2:6, 1],
                 summary(fit_lme_ca_t)$coefficients[1, 1], sum(summary(fit_lme_ca_d)$coefficients[2:6, 1])/5, summary(fit_lme_ca_d)$coefficients[2:6, 1],
                 summary(fit_lme_cac_t)$coefficients[1, 1], sum(summary(fit_lme_cac_d)$coefficients[2:6, 1])/5, summary(fit_lme_cac_d)$coefficients[2:6, 1])
  naive_se <- c(summary(fit_lme_un_t)$coefficients[1, 2], sqrt(sum(summary(fit_lme_un_d)$vcov[2:6, 2:6]))/5, summary(fit_lme_un_d)$coefficients[2:6, 2],
                summary(fit_lme_ca_t)$coefficients[1, 2], sqrt(sum(summary(fit_lme_ca_d)$vcov[2:6, 2:6]))/5, summary(fit_lme_ca_d)$coefficients[2:6, 2],
                summary(fit_lme_cac_t)$coefficients[1, 2], sqrt(sum(summary(fit_lme_cac_d)$vcov[2:6, 2:6]))/5, summary(fit_lme_cac_d)$coefficients[2:6, 2])
  
  naive_se2 <- c(sqrt(fit_un_t$naive[1, 1]), sqrt(sum(fit_un_d$naive[2:6, 2:6]))/5, sqrt(diag(fit_un_d$naive)[2:6]),
                 sqrt(fit_ca_t$naive[1, 1]), sqrt(sum(fit_ca_d$naive[2:6, 2:6]))/5, sqrt(diag(fit_ca_d$naive)[2:6]),
                 sqrt(fit_cac_t$naive[1, 1]), sqrt(sum(fit_cac_d$naive[2:6, 2:6]))/5, sqrt(diag(fit_cac_d$naive)[2:6]))
  robust_se <- c(sqrt(fit_un_t$robust[1, 1]), sqrt(sum(fit_un_d$robust[2:6, 2:6]))/5, sqrt(diag(fit_un_d$robust)[2:6]),
                 sqrt(fit_ca_t$robust[1, 1]), sqrt(sum(fit_ca_d$robust[2:6, 2:6]))/5, sqrt(diag(fit_ca_d$robust)[2:6]),
                 sqrt(fit_cac_t$robust[1, 1]), sqrt(sum(fit_cac_d$robust[2:6, 2:6]))/5, sqrt(diag(fit_cac_d$robust)[2:6]))
  KC_se <- c(sqrt(fit_un_t$varKC[1, 1]), sqrt(sum(fit_un_d$varKC[2:6, 2:6]))/5, sqrt(diag(fit_un_d$varKC)[2:6]),
             sqrt(fit_ca_t$varKC[1, 1]), sqrt(sum(fit_ca_d$varKC[2:6, 2:6]))/5, sqrt(diag(fit_ca_d$varKC)[2:6]),
             sqrt(fit_cac_t$varKC[1, 1]), sqrt(sum(fit_cac_d$varKC[2:6, 2:6]))/5, sqrt(diag(fit_cac_d$varKC)[2:6]))
  MD_se <- c(sqrt(fit_un_t$varMD[1, 1]), sqrt(sum(fit_un_d$varMD[2:6, 2:6]))/5, sqrt(diag(fit_un_d$varMD)[2:6]),
             sqrt(fit_ca_t$varMD[1, 1]), sqrt(sum(fit_ca_d$varMD[2:6, 2:6]))/5, sqrt(diag(fit_ca_d$varMD)[2:6]),
             sqrt(fit_cac_t$varMD[1, 1]), sqrt(sum(fit_cac_d$varMD[2:6, 2:6]))/5, sqrt(diag(fit_cac_d$varMD)[2:6]))
  
  robust_se_o <- c(c_un_t[2], c_un_d[2], c_un_d[3:7],
                   c_ca_t[2], c_ca_d[2], c_ca_d[3:7],
                   c_cac_t[2], c_cac_d[2], c_cac_d[3:7])
  robust_se_c <- c(c_un_t[3], c_un_d[8], c_un_d[9:13],
                   c_ca_t[3], c_ca_d[8], c_ca_d[9:13],
                   c_cac_t[3], c_cac_d[8], c_cac_d[9:13])
  
  alphas <- c(alpha_un_t, alpha_un_d, alpha_ca_t, alpha_ca_d, alpha_cac_t, alpha_cac_d)
  
  tau2 <- c(summary(fit_lme_un_t)$varcor$CLUSTER[1], summary(fit_lme_un_d)$varcor$CLUSTER[1], 
            summary(fit_lme_ca_t)$varcor$CLUSTER[1], summary(fit_lme_ca_d)$varcor$CLUSTER[1], 
            summary(fit_lme_cac_t)$varcor$CLUSTER[1], summary(fit_lme_cac_d)$varcor$CLUSTER[1])
  sigma2 <- c(summary(fit_lme_un_t)$sigma^2, summary(fit_lme_un_d)$sigma^2, 
              summary(fit_lme_ca_t)$sigma^2, summary(fit_lme_ca_d)$sigma^2, 
              summary(fit_lme_cac_t)$sigma^2, summary(fit_lme_cac_d)$sigma^2)
  
  ###################################### nested exchangeable correlation ######################################
  dfo <- dfo %>% mutate(CLUSTER_PERIOD = paste0(CLUSTER, PERIOD))
  # fit model 1t
  fit_lme_un_t <- lmer(Y ~ 0 + TRT + PERIOD + (1 | CLUSTER) + (1 | CLUSTER_PERIOD), dfo)
  beta <- as.numeric(summary(fit_lme_un_t)$coefficients[, 1])
  X <- model.matrix(~ 0 + TRT + PERIOD, data = dfo)
  
  tau2 <- summary(fit_lme_un_t)$varcor$CLUSTER[1]
  kappa2 <-  summary(fit_lme_un_t)$varcor$CLUSTER_PERIOD[1]
  sigma2 <- summary(fit_lme_un_t)$sigma^2
  c_un_t <- SV_NE(y, X, id, beta, tau2, sigma2, 1, kappa2, period_o)
  
  # fit model 1d
  fit_lme_un_d <- lmer(Y ~ D + PERIOD + (1 | CLUSTER) + (1 | CLUSTER_PERIOD), dfo)
  beta <- as.numeric(summary(fit_lme_un_d)$coefficients[, 1])
  X <- model.matrix(~ D + PERIOD, data = dfo)
  
  tau2 <- summary(fit_lme_un_d)$varcor$CLUSTER[1]
  kappa2 <-  summary(fit_lme_un_d)$varcor$CLUSTER_PERIOD[1]
  sigma2 <- summary(fit_lme_un_d)$sigma^2
  c_un_d <- SV_NE(y, X, id, beta, tau2, sigma2, 5, kappa2, period_o)
  
  # fit model 2t
  fit_lme_ca_t <- lmer(Y ~ 0 + TRT + PERIOD + X1 + X3 + (1 | CLUSTER) + (1 | CLUSTER_PERIOD), dfo)
  beta <- as.numeric(summary(fit_lme_ca_t)$coefficients[, 1])
  X <- model.matrix(~ 0 + TRT + PERIOD + X1 + X3, data = dfo)
  
  tau2 <- summary(fit_lme_ca_t)$varcor$CLUSTER[1]
  kappa2 <-  summary(fit_lme_ca_t)$varcor$CLUSTER_PERIOD[1]
  sigma2 <- summary(fit_lme_ca_t)$sigma^2
  c_ca_t <- SV_NE(y, X, id, beta, tau2, sigma2, 1, kappa2, period_o)
  
  # fit model 2d
  fit_lme_ca_d <- lmer(Y ~ D + PERIOD + X1 + X3 + (1 | CLUSTER) + (1 | CLUSTER_PERIOD), dfo)
  beta <- as.numeric(summary(fit_lme_ca_d)$coefficients[, 1])
  X <- model.matrix(~ D + PERIOD + X1 + X3, data = dfo)
  
  tau2 <- summary(fit_lme_ca_d)$varcor$CLUSTER[1]
  kappa2 <-  summary(fit_lme_ca_d)$varcor$CLUSTER_PERIOD[1]
  sigma2 <- summary(fit_lme_ca_d)$sigma^2
  c_ca_d <- SV_NE(y, X, id, beta, tau2, sigma2, 5, kappa2, period_o)
  
  # fit model 3t
  fit_lme_cac_t <- lmer(Y ~ 0 + TRT + PERIOD + X1 + X2 + X3 + X4 + (1 | CLUSTER) + (1 | CLUSTER_PERIOD), dfo)
  beta <- as.numeric(summary(fit_lme_cac_t)$coefficients[, 1])
  X <- model.matrix(~ 0 + TRT + PERIOD + X1 + X2 + X3 + X4, data = dfo)
   
  tau2 <- summary(fit_lme_cac_t)$varcor$CLUSTER[1]
  kappa2 <-  summary(fit_lme_cac_t)$varcor$CLUSTER_PERIOD[1]
  sigma2 <- summary(fit_lme_cac_t)$sigma^2
  c_cac_t <- SV_NE(y, X, id, beta, tau2, sigma2, 1, kappa2, period_o)
  
  # fit model 3d
  fit_lme_cac_d <- lmer(Y ~ D + PERIOD + X1 + X2 + X3 + X4 + (1 | CLUSTER) + (1 | CLUSTER_PERIOD), dfo)
  beta <- as.numeric(summary(fit_lme_cac_d)$coefficients[, 1])
  X <- model.matrix(~ D + PERIOD + X1 + X2 + X3 + X4, data = dfo)
  
  tau2 <- summary(fit_lme_cac_d)$varcor$CLUSTER[1]
  kappa2 <-  summary(fit_lme_cac_d)$varcor$CLUSTER_PERIOD[1]
  sigma2 <- summary(fit_lme_cac_d)$sigma^2
  c_cac_d <- SV_NE(y, X, id, beta, tau2, sigma2, 5, kappa2, period_o)
  
  # estimates
  estimates_ne <- c(summary(fit_lme_un_t)$coefficients[1, 1], sum(summary(fit_lme_un_d)$coefficients[2:6, 1])/5, summary(fit_lme_un_d)$coefficients[2:6, 1],
                 summary(fit_lme_ca_t)$coefficients[1, 1], sum(summary(fit_lme_ca_d)$coefficients[2:6, 1])/5, summary(fit_lme_ca_d)$coefficients[2:6, 1],
                 summary(fit_lme_cac_t)$coefficients[1, 1], sum(summary(fit_lme_cac_d)$coefficients[2:6, 1])/5, summary(fit_lme_cac_d)$coefficients[2:6, 1])
  naive_se_ne <- c(summary(fit_lme_un_t)$coefficients[1, 2], sqrt(sum(summary(fit_lme_un_d)$vcov[2:6, 2:6]))/5, summary(fit_lme_un_d)$coefficients[2:6, 2],
                summary(fit_lme_ca_t)$coefficients[1, 2], sqrt(sum(summary(fit_lme_ca_d)$vcov[2:6, 2:6]))/5, summary(fit_lme_ca_d)$coefficients[2:6, 2],
                summary(fit_lme_cac_t)$coefficients[1, 2], sqrt(sum(summary(fit_lme_cac_d)$vcov[2:6, 2:6]))/5, summary(fit_lme_cac_d)$coefficients[2:6, 2])
  
  robust_se_o_ne <- c(c_un_t[2], c_un_d[2], c_un_d[3:7],
                   c_ca_t[2], c_ca_d[2], c_ca_d[3:7],
                   c_cac_t[2], c_cac_d[2], c_cac_d[3:7])
  robust_se_c_ne <- c(c_un_t[3], c_un_d[8], c_un_d[9:13],
                   c_ca_t[3], c_ca_d[8], c_ca_d[9:13],
                   c_cac_t[3], c_cac_d[8], c_cac_d[9:13])
  
  tau2_ne <- c(summary(fit_lme_un_t)$varcor$CLUSTER[1], summary(fit_lme_un_d)$varcor$CLUSTER[1], 
            summary(fit_lme_ca_t)$varcor$CLUSTER[1], summary(fit_lme_ca_d)$varcor$CLUSTER[1], 
            summary(fit_lme_cac_t)$varcor$CLUSTER[1], summary(fit_lme_cac_d)$varcor$CLUSTER[1])
  kappa2_ne <- c(summary(fit_lme_un_t)$varcor$CLUSTER_PERIOD[1], summary(fit_lme_un_d)$varcor$CLUSTER_PERIOD[1], 
              summary(fit_lme_ca_t)$varcor$CLUSTER_PERIOD[1], summary(fit_lme_ca_d)$varcor$CLUSTER_PERIOD[1], 
              summary(fit_lme_cac_t)$varcor$CLUSTER_PERIOD[1], summary(fit_lme_cac_d)$varcor$CLUSTER_PERIOD[1])
  sigma2_ne <- c(summary(fit_lme_un_t)$sigma^2, summary(fit_lme_un_d)$sigma^2, 
              summary(fit_lme_ca_t)$sigma^2, summary(fit_lme_ca_d)$sigma^2, 
              summary(fit_lme_cac_t)$sigma^2, summary(fit_lme_cac_d)$sigma^2)
  
  # estimand
  estimand <- c(mean(yd[, -1] - yd[, 1]), mean(yd[, -1] - yd[, 1]))
  for (j in 1:(t-1)){
    estimand <- c(estimand, mean(yd[, (j+1)] - yd[, 1]))
  }

  # results [(7+3*7*8+1+1)+(7+3*7*4+6+6+6)]
  return(c(estimand, estimates, naive_se, naive_se2, robust_se, KC_se, MD_se, robust_se_o, robust_se_c, tau2, sigma2,
           estimand, estimates_ne, naive_se_ne, robust_se_o_ne, robust_se_c_ne, tau2_ne, kappa2_ne, sigma2_ne))
}

simulation <- function(n, t, pplv, cpl, cpu, beta0, beta1, mu_bx, sd_cxc, sd_cxe, sdb, sdep, nsim){
  results <- NULL
  nerror <- 0
  sj <- 0
  for (ss in 1:nsim){
    passMODEL <- 1
    
    while(passMODEL == 1){
      sj <- sj+1
      set.seed(sj)
      res_ss <- try(ATE(n=n, t=t, pplv=pplv, cpl=cpl, cpu=cpu, beta0=beta0, beta1=beta1, mu_bx=mu_bx, sd_cxc=sd_cxc, sd_cxe=sd_cxe, sdb=sdb, sdep=sdep), silent=TRUE)
      if(class(res_ss)[1]=="try-error" | is.null(res_ss)==1){passMODEL <- 1; nerror <- nerror + 1; next}
      passMODEL <- 0
    }
    
    results <- rbind(results, res_ss)
    
    # Loop control
    if(ss %% 50 == 0){print(paste("Iteration", ss, "NERROR =", nerror))}
  }
  return(results)
}

simuRES <- function(data_all, n, t){
  ###################################### simple exchangeable correlation ######################################
  data <- data_all[, 1:177]
  
  # estimand
  estimand <- colMeans(data[, 1:7])
  estimand_3 <-  c(estimand, estimand, estimand)
  
  # estimate
  estimate <- colMeans(data[, 8:28])
  
  # bias
  bias <- estimate - estimand
  relative_bias <- ((estimate - estimand)/estimand)*100
  
  # empirical standard error
  ese <- apply(data[, 8:28], 2, sd)
  
  # relative efficiency
  re_ca <- ese[1:7]^2/ese[8:14]^2
  re_cac <- ese[1:7]^2/ese[15:21]^2
  re <- c(rep(1, 7), re_ca, re_cac)
  
  ase <- ec <- NULL
  for (se in 1:7){
    datase <- data[, c(1:28, (29+(se-1)*21):(28+se*21))]
    
    asese <- ecse <- NULL
    for (i in 1:3){
      datai <- datase[, c((8+(i-1)*7):(7+i*7), (29+(i-1)*7):(28+i*7))]
      
      # averaged standard error
      asese <- c(asese, apply(datai[, 8:14], 2, mean))
      
      # empirical coverage using t value with n-(1+t+p) degrees of freedom
      df_t <- 1 + t + 2*(i-1)
      df_d <- 2*t + 2*(i-1)
      
      eci <- mean((datai[, 1]-qt(0.975, df=n-df_t)*datai[, 8]<=estimand[1])+(datai[, 1]+qt(0.975, df=n-df_t)*datai[, 8]>=estimand[1]) == 2, na.rm=TRUE)
      for (j in 2:7){
        eci <- c(eci, mean((datai[, j]-qt(0.975, df=n-df_d)*datai[, (7+j)]<=estimand[j])+(datai[, j]+qt(0.975, df=n-df_d)*datai[, (7+j)]>=estimand[j]) == 2, na.rm=TRUE))
      }
      
      ecse <- c(ecse, eci)
    }

    ase <- cbind(ase, asese)
    ec <- cbind(ec, ecse)
  }
  
  res_se <- cbind(estimand_3, estimate, bias, relative_bias, ese, re, ase, ec)
  
  ###################################### nested exchangeable correlation ######################################
  data <- data_all[, 178:286]
  
  # estimand
  estimand <- colMeans(data[, 1:7])
  estimand_3 <-  c(estimand, estimand, estimand)
  
  # estimate
  estimate <- colMeans(data[, 8:28])
  
  # bias
  bias <- estimate - estimand
  relative_bias <- ((estimate - estimand)/estimand)*100
  
  # empirical standard error
  ese <- apply(data[, 8:28], 2, sd)
  
  # relative efficiency
  re_ca <- ese[1:7]^2/ese[8:14]^2
  re_cac <- ese[1:7]^2/ese[15:21]^2
  re <- c(rep(1, 7), re_ca, re_cac)
  
  ase <- ec <- NULL
  for (se in 1:3){
    datase <- data[, c(1:28, (29+(se-1)*21):(28+se*21))]
    
    asese <- ecse <- NULL
    for (i in 1:3){
      datai <- datase[, c((8+(i-1)*7):(7+i*7), (29+(i-1)*7):(28+i*7))]
      
      # averaged standard error
      asese <- c(asese, apply(datai[, 8:14], 2, mean))
      
      # empirical coverage using t value with n-(1+t+p) degrees of freedom
      df_t <- 1 + t + 2*(i-1)
      df_d <- 2*t + 2*(i-1)
      
      eci <- mean((datai[, 1]-qt(0.975, df=n-df_t)*datai[, 8]<=estimand[1])+(datai[, 1]+qt(0.975, df=n-df_t)*datai[, 8]>=estimand[1]) == 2, na.rm=TRUE)
      for (j in 2:7){
        eci <- c(eci, mean((datai[, j]-qt(0.975, df=n-df_d)*datai[, (7+j)]<=estimand[j])+(datai[, j]+qt(0.975, df=n-df_d)*datai[, (7+j)]>=estimand[j]) == 2, na.rm=TRUE))
      }
      
      ecse <- c(ecse, eci)
    }
    
    ase <- cbind(ase, asese)
    ec <- cbind(ec, ecse)
  }
  
  res_ne <- cbind(estimand_3, estimate, bias, relative_bias, ese, re, ase, ec)
  
  # save results
  res <- cbind(res_se, res_ne)
  colnames(res) <- c("ATE", "Estimate", "ABias", "ARBias", "EmpSE", "RE", 
                     "ASE_MB", "ASE_MB2", "ASE_ROB", "ASE_KC", "ASE_MD", "ASE_ROB_O", "ASE_ROB_C", 
                     "ECP_MB", "ECP_MB2", "ECP_ROB", "ECP_KC", "ECP_MD", "ECP_ROB_O", "ECP_ROB_C",
                     "ATE_NE", "Estimate_NE", "ABias_NE", "ARBias_NE", "EmpSE_NE", "RE_NE", 
                     "ASE_MB_NE", "ASE_ROB_O_NE", "ASE_ROB_C_NE", 
                     "ECP_MB_NE", "ECP_ROB_O_NE", "ECP_ROB_C_NE")
  rownames(res) <- rep(c("TRT", "D_AVG", "D1", "D2", "D3", "D4", "D5"), times=3)
  
  return(res)
}

# ---------------------------------------------------------------------------------


# Scenario 2 (30 clusters)
n <- 30; t <- 6
pplv <- rep(1000, n)
cpl <- 5; cpu <- 50
beta0 <- c(1, 1.5, 2, 2.5, 3); beta1 <- beta0/2
mu_bx <- c(0.5, 0.8)
sd_cxc <- c(sqrt(0.1), sqrt(0.1)); sd_cxe <- c(sqrt(0.4), sqrt(0.9))
sdb <- c(sqrt(0.06), sqrt(0.04)); sdep <- sqrt(0.9)
nsim <- 1000

cres2_30 <- simulation(n, t, pplv, cpl, cpu, beta0, beta1, mu_bx, sd_cxc, sd_cxe, sdb, sdep, nsim)
save(cres2_30, file = "Data Results B/cres2_30.RData")
res2 <- simuRES(cres2_30, n, t)
res2

# Scenario 2 (100 clusters)
n <- 100; t <- 6
pplv <- rep(1000, n)
cpl <- 5; cpu <- 50
beta0 <- c(1, 1.5, 2, 2.5, 3); beta1 <- beta0/2
mu_bx <- c(0.5, 0.8)
sd_cxc <- c(sqrt(0.1), sqrt(0.1)); sd_cxe <- c(sqrt(0.4), sqrt(0.9))
sdb <- c(sqrt(0.06), sqrt(0.04)); sdep <- sqrt(0.9)
nsim <- 1000

cres2_100 <- simulation(n, t, pplv, cpl, cpu, beta0, beta1, mu_bx, sd_cxc, sd_cxe, sdb, sdep, nsim)
save(cres2_100, file = "Data Results B/cres2_100.RData")
res2 <- simuRES(cres2_100, n, t)
res2


# ---------------------------------------------------------------------------------


