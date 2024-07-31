#############################################################
# Simulation for SWD (with Data Generation Process A2)

# September 2022

# INPUT
# n: Total number of clusters
# t: period
# pplv: Vector of population size in each cluster
# cpl: Lower bound of the cluster-period sizes
# cpu: Upper bound of the cluster-period sizes
# beta0: Parameter of the constant intervention effect
# beta1: Parameter of the intervention effect related to covariates
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
  beta <- beta0 + beta1/2*period*(df$X1-df$X1_mean) + beta1*period/t*(df$X3^3-df$X33_mean)
  
  df$Y0 <- pv + bv + cv + epv + xs
  df$Y1 <- pv + bv + cv + epv + xs + beta + ciev
  
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
  
  # indicator of intervention arm; i = 1 treated; i = 0 non treated
  trt_rd <- sample(rep(1:(t-1), each=n/(t-1)))
  design <- matrix(0, t-1, t)
  design[upper.tri(design)] <- 1
  trt_seq <- NULL
  for (i in 1:n){
    trt_seq <- rbind(trt_seq, rep(design[trt_rd[i], ]))
  }
  dfo$TRT <- rep(c(t(trt_seq)), times = mv)
  # observed values
  dfo$Y <- dfo$TRT*dfo$Y1 + (1 - dfo$TRT)*dfo$Y0
  y <- dfo$Y
  id <- dfo$CLUSTER
  period_o <- dfo$PERIOD
  phi <- var(y)
  
  ###################################### simple exchangeable correlation ######################################
  # fit model 1
  fit_lme_un <- lmer(Y ~ 0 + TRT + PERIOD + (1 | CLUSTER), data = dfo)
  beta <- as.numeric(summary(fit_lme_un)$coefficients[, 1])
  X <- model.matrix(~ 0 + TRT + PERIOD, data = dfo)
  alpha_un <- icc(fit_lme_un)$ICC_adjusted
  fit_un <- contGEEBCV(y=y,X=X,beta=beta,alpha=alpha_un,phi=phi,id=id)
  
  tau2 <- summary(fit_lme_un)$varcor$CLUSTER[1]
  sigma2 <- summary(fit_lme_un)$sigma^2
  c_un <- SV(y, X, id, beta, tau2, sigma2)
  
  # fit model 2
  fit_lme_ca <- lmer(Y ~ 0 + TRT + PERIOD + X1 + X3 + (1 | CLUSTER), data = dfo)
  beta <- as.numeric(summary(fit_lme_ca)$coefficients[, 1])
  X <- model.matrix(~ 0 + TRT + PERIOD + X1 + X3, data = dfo)
  alpha_ca <- icc(fit_lme_ca)$ICC_adjusted
  fit_ca <- contGEEBCV(y=y,X=X,beta=beta,alpha=alpha_ca,phi=phi,id=id)
  
  tau2 <- summary(fit_lme_ca)$varcor$CLUSTER[1]
  sigma2 <- summary(fit_lme_ca)$sigma^2
  c_ca <- SV(y, X, id, beta, tau2, sigma2)
  
  # fit model 3
  fit_lme_cac <- lmer(Y ~ 0 + TRT + PERIOD + X1 + X2 + X3 + X4 + (1 | CLUSTER), data = dfo)
  beta <- as.numeric(summary(fit_lme_cac)$coefficients[, 1])
  X <- model.matrix(~ 0 + TRT + PERIOD + X1 + X2 + X3 + X4, data = dfo)
  alpha_cac <- icc(fit_lme_cac)$ICC_adjusted
  fit_cac <- contGEEBCV(y=y,X=X,beta=beta,alpha=alpha_cac,phi=phi,id=id)
  
  tau2 <- summary(fit_lme_cac)$varcor$CLUSTER[1]
  sigma2 <- summary(fit_lme_cac)$sigma^2
  c_cac <- SV(y, X, id, beta, tau2, sigma2)
  
  # estimates
  estimates <- c(summary(fit_lme_un)$coefficients[1, 1], summary(fit_lme_ca)$coefficients[1, 1], summary(fit_lme_cac)$coefficients[1, 1])
  naive_se <-c(summary(fit_lme_un)$coefficients[1, 2], summary(fit_lme_ca)$coefficients[1, 2], summary(fit_lme_cac)$coefficients[1, 2])
  
  naive_se2 <- sqrt(c(fit_un$naive[1, 1], fit_ca$naive[1, 1], fit_cac$naive[1, 1]))
  robust_se <- sqrt(c(fit_un$robust[1, 1], fit_ca$robust[1, 1], fit_cac$robust[1, 1]))
  KC_se <- sqrt(c(fit_un$varKC[1, 1], fit_ca$varKC[1, 1], fit_cac$varKC[1, 1]))
  MD_se <- sqrt(c(fit_un$varMD[1, 1], fit_ca$varMD[1, 1], fit_cac$varMD[1, 1]))
  
  robust_se_o <- c(c_un[2], c_ca[2], c_cac[2])
  robust_se_c <- c(c_un[3], c_ca[3], c_cac[3])
  
  alphas <- c(alpha_un, alpha_ca, alpha_cac)
  
  tau2 <- c(summary(fit_lme_un)$varcor$CLUSTER[1], summary(fit_lme_ca)$varcor$CLUSTER[1], summary(fit_lme_cac)$varcor$CLUSTER[1])
  sigma2 <- c(summary(fit_lme_un)$sigma^2, summary(fit_lme_ca)$sigma^2, summary(fit_lme_cac)$sigma^2)
  
  ###################################### nested exchangeable correlation ######################################
  dfo <- dfo %>% mutate(CLUSTER_PERIOD = paste0(CLUSTER, PERIOD))
  # fit model 1
  fit_lme_un <- lmer(Y ~ 0 + TRT + PERIOD + (1 | CLUSTER) + (1 | CLUSTER_PERIOD), dfo)
  beta <- as.numeric(summary(fit_lme_un)$coefficients[, 1])
  X <- model.matrix(~ 0 + TRT + PERIOD, data = dfo)
  
  tau2 <- summary(fit_lme_un)$varcor$CLUSTER[1]
  kappa2 <-  summary(fit_lme_un)$varcor$CLUSTER_PERIOD[1]
  sigma2 <- summary(fit_lme_un)$sigma^2
  c_un <- SV_NE(y, X, id, beta, tau2, sigma2, 1, kappa2, period_o)
  
  # fit model 2
  fit_lme_ca <- lmer(Y ~ 0 + TRT + PERIOD + X1 + X3 + (1 | CLUSTER) + (1 | CLUSTER_PERIOD), dfo)
  beta <- as.numeric(summary(fit_lme_ca)$coefficients[, 1])
  X <- model.matrix(~ 0 + TRT + PERIOD + X1 + X3, data = dfo)
  
  tau2 <- summary(fit_lme_ca)$varcor$CLUSTER[1]
  kappa2 <-  summary(fit_lme_ca)$varcor$CLUSTER_PERIOD[1]
  sigma2 <- summary(fit_lme_ca)$sigma^2
  c_ca <- SV_NE(y, X, id, beta, tau2, sigma2, 1, kappa2, period_o)
  
  # fit model 3
  fit_lme_cac <- lmer(Y ~ 0 + TRT + PERIOD + X1 + X2 + X3 + X4 + (1 | CLUSTER) + (1 | CLUSTER_PERIOD), dfo)
  beta <- as.numeric(summary(fit_lme_cac)$coefficients[, 1])
  X <- model.matrix(~ 0 + TRT + PERIOD + X1 + X2 + X3 + X4, data = dfo)
  
  tau2 <- summary(fit_lme_cac)$varcor$CLUSTER[1]
  kappa2 <-  summary(fit_lme_cac)$varcor$CLUSTER_PERIOD[1]
  sigma2 <- summary(fit_lme_cac)$sigma^2
  c_cac <- SV_NE(y, X, id, beta, tau2, sigma2, 1, kappa2, period_o)
  
  # estimates
  estimates_ne <- c(summary(fit_lme_un)$coefficients[1, 1], summary(fit_lme_ca)$coefficients[1, 1], summary(fit_lme_cac)$coefficients[1, 1])
  naive_se_ne <-c(summary(fit_lme_un)$coefficients[1, 2], summary(fit_lme_ca)$coefficients[1, 2], summary(fit_lme_cac)$coefficients[1, 2])
  
  robust_se_o_ne <- c(c_un[2], c_ca[2], c_cac[2])
  robust_se_c_ne <- c(c_un[3], c_ca[3], c_cac[3])
  
  tau2_ne <- c(summary(fit_lme_un)$varcor$CLUSTER[1], summary(fit_lme_ca)$varcor$CLUSTER[1], summary(fit_lme_cac)$varcor$CLUSTER[1])
  kappa2_ne <- c(summary(fit_lme_un)$varcor$CLUSTER_PERIOD[1], summary(fit_lme_ca)$varcor$CLUSTER_PERIOD[1], summary(fit_lme_cac)$varcor$CLUSTER_PERIOD[1])
  sigma2_ne <- c(summary(fit_lme_un)$sigma^2, summary(fit_lme_ca)$sigma^2, summary(fit_lme_cac)$sigma^2)
  
  # estimand
  yd <- df$Y1 - df$Y0
  estimand <- sum(yd)/sum(mvt)
  ydo <- dfo$Y1 - dfo$Y0
  estimando <- sum(ydo)/sum(mv)
  
  # results [(2+3*8+3+1+1+1)+(2+3*4+3+3+3)]
  return(c(estimand, estimando, estimates, naive_se, naive_se2, robust_se, KC_se, MD_se, robust_se_o, robust_se_c, alphas, phi, tau2, sigma2,
           estimand, estimando, estimates_ne, naive_se_ne, robust_se_o_ne, robust_se_c_ne, tau2_ne, kappa2_ne, sigma2_ne))
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
  data <- data_all[, 1:32]
  
  # estimand
  estimand <- mean(data[, 1])
  estimando <- mean(data[, 2])
  
  # estimate
  estimate <- colMeans(data[, 3:5])
  
  # bias
  bias <- estimate - estimand
  relative_bias <- ((estimate - estimand)/estimand)*100
  
  # empirical standard error
  ese <- apply(data[, 3:5], 2, sd)
  
  # relative efficiency
  re_ca <- ese[1]^2/ese[2]^2
  re_cac <- ese[1]^2/ese[3]^2
  
  # averaged icc & variance
  icc <- colMeans(data[, 27:29])
  phi <- mean(data[, 30])
  
  ase <- ec <- NULL
  for (i in 1:7){
    datai <- data[, c(3:5, (6+(i-1)*3):(5+i*3))]
    
    # averaged standard error
    ase_i <- apply(datai[, 4:6], 2, mean)
    
    # empirical coverage using t value with n-(1+t+p) degrees of freedom
    df_un <- 1 + t
    ec_un <- mean((datai[, 1]-qt(0.975, df=n-df_un)*datai[, 4]<=estimand)+(datai[, 1]+qt(0.975, df=n-df_un)*datai[, 4]>=estimand) == 2, na.rm=TRUE)
    
    df_ca <- 1 + t + 2
    ec_ca <- mean((datai[, 2]-qt(0.975, df=n-df_ca)*datai[, 5]<=estimand)+(datai[, 2]+qt(0.975, df=n-df_ca)*datai[, 5]>=estimand) == 2, na.rm=TRUE)
    
    df_cac <- 1 + t + 4
    ec_cac <- mean((datai[, 3]-qt(0.975, df=n-df_cac)*datai[, 6]<=estimand)+(datai[, 3]+qt(0.975, df=n-df_cac)*datai[, 6]>=estimand) == 2, na.rm=TRUE)
    
    eci <- c(ec_un, ec_ca, ec_cac)
    
    ase <- cbind(ase, ase_i)
    ec <- cbind(ec, eci)
  }
  
  res_se <- cbind(estimand, estimando, estimate, bias, relative_bias, ese, c(1, re_ca, re_cac), ase, ec, icc, phi)
  
  ###################################### nested exchangeable correlation ######################################
  data <- data_all[, 33:55]
  
  # estimand
  estimand <- mean(data[, 1])
  estimando <- mean(data[, 2])
  
  # estimate
  estimate <- colMeans(data[, 3:5])
  
  # bias
  bias <- estimate - estimand
  relative_bias <- ((estimate - estimand)/estimand)*100
  
  # empirical standard error
  ese <- apply(data[, 3:5], 2, sd)
  
  # relative efficiency
  re_ca <- ese[1]^2/ese[2]^2
  re_cac <- ese[1]^2/ese[3]^2
  
  ase <- ec <- NULL
  for (i in 1:3){
    datai <- data[, c(3:5, (6+(i-1)*3):(5+i*3))]
    
    # averaged standard error
    ase_i <- apply(datai[, 4:6], 2, mean)
    
    # empirical coverage using t value with n-(1+t+p) degrees of freedom
    df_un <- 1 + t
    ec_un <- mean((datai[, 1]-qt(0.975, df=n-df_un)*datai[, 4]<=estimand)+(datai[, 1]+qt(0.975, df=n-df_un)*datai[, 4]>=estimand) == 2, na.rm=TRUE)
    
    df_ca <- 1 + t + 2
    ec_ca <- mean((datai[, 2]-qt(0.975, df=n-df_ca)*datai[, 5]<=estimand)+(datai[, 2]+qt(0.975, df=n-df_ca)*datai[, 5]>=estimand) == 2, na.rm=TRUE)
    
    df_cac <- 1 + t + 4
    ec_cac <- mean((datai[, 3]-qt(0.975, df=n-df_cac)*datai[, 6]<=estimand)+(datai[, 3]+qt(0.975, df=n-df_cac)*datai[, 6]>=estimand) == 2, na.rm=TRUE)
    
    eci <- c(ec_un, ec_ca, ec_cac)
    
    ase <- cbind(ase, ase_i)
    ec <- cbind(ec, eci)
  }
  
  res_ne <- cbind(estimand, estimando, estimate, bias, relative_bias, ese, c(1, re_ca, re_cac), ase, ec)
  
  # save results
  res <- cbind(res_se, res_ne)
  colnames(res) <- c("ATE", "ATEo", "Estimate", "ABias", "ARBias", "EmpSE", "RE", 
                     "ASE_MB", "ASE_MB2", "ASE_ROB", "ASE_KC", "ASE_MD", "ASE_ROB_O", "ASE_ROB_C", 
                     "ECP_MB", "ECP_MB2", "ECP_ROB", "ECP_KC", "ECP_MD", "ECP_ROB_O", "ECP_ROB_C", 
                     "ICC", "VAR(Y)",
                     "ATE_NE", "ATEo_NE", "Estimate_NE", "ABias_NE", "ARBias_NE", "EmpSE_NE", "RE_NE", 
                     "ASE_MB_NE", "ASE_ROB_O_NE", "ASE_ROB_C_NE",
                     "ECP_MB_NE", "ECP_ROB_O_NE", "ECP_ROB_C_NE")
  rownames(res) <- c("mixed-model unadjusted", "mixed-model adjusted (incomplete)", "mixed-model adjusted (complete)")
  
  return(res)
}

# ---------------------------------------------------------------------------------


# Scenario 2 (30 clusters)
n <- 30; t <- 6
pplv <- rep(1000, n)
cpl <- 5; cpu <- 50
beta0 <- 2; beta1 <- 1
mu_bx <- c(0.5, 0.8)
sd_cxc <- c(sqrt(0.1), sqrt(0.1)); sd_cxe <- c(sqrt(0.4), sqrt(0.9))
sdb <- c(sqrt(0.06), sqrt(0.04)); sdep <- sqrt(0.9)
nsim <- 1000

cres2_30 <- simulation(n, t, pplv, cpl, cpu, beta0, beta1, mu_bx, sd_cxc, sd_cxe, sdb, sdep, nsim)
save(cres2_30, file = "Data Results A/cres2_30.RData")
res2 <- simuRES(cres2_30, n, t)
res2

# Scenario 2 (100 clusters)
n <- 100; t <- 6
pplv <- rep(1000, n)
cpl <- 5; cpu <- 50
beta0 <- 2; beta1 <- 1
mu_bx <- c(0.5, 0.8)
sd_cxc <- c(sqrt(0.1), sqrt(0.1)); sd_cxe <- c(sqrt(0.4), sqrt(0.9))
sdb <- c(sqrt(0.06), sqrt(0.04)); sdep <- sqrt(0.9)
nsim <- 1000

cres2_100 <- simulation(n, t, pplv, cpl, cpu, beta0, beta1, mu_bx, sd_cxc, sd_cxe, sdb, sdep, nsim)
save(cres2_100, file = "Data Results A/cres2_100.RData")
res2 <- simuRES(cres2_100, n, t)
res2


# ---------------------------------------------------------------------------------


