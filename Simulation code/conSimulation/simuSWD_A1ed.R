#############################################################
# Simulation for SWD (with Data Generation Process A1)

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
library(glmmTMB)

ATE <- function(n, t, pplv, cpl, cpu, beta0, beta1, mu_bx, sd_cxc, sd_cxe, sdb, sdep){
  
  # Source program code
  source("contGEE_BCV2_ED.R")
  
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
  bv <- rep(rep(rnorm(n, 0, sdb), each=t), times=mvt)
  # random errors
  epv <-NULL
  for (i in 1:n){
    epv <- c(epv, rep(rnorm(pplv[i], 0, sdep), times=t))
  }
  
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
  beta <- beta0 + beta1/2*(df$X1-df$X1_mean) + beta1*(df$X3^3-df$X33_mean)
  
  df$Y0 <- pv + bv + epv + xs
  df$Y1 <- pv + bv + epv + xs + beta
  
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
  
  ###################################### exponential decay correlation ######################################
  
  # fit model 1
  fit_lme_un <- glmmTMB(Y ~  0 + TRT + PERIOD + ar1(PERIOD + 0 | CLUSTER), REML=T, 
                        data = dfo, family = gaussian)
  V_AR_un <- VarCorr(fit_lme_un)
  wpicc_ar_un <- round(V_AR_un$cond$CLUSTER[1,1]/sum(V_AR_un$cond$CLUSTER[1,1],(attr(V_AR_un$cond,"sc"))^2),8)
  ar1_un <- round(attr(V_AR_un$cond$CLUSTER,"correlation")[1,2],8)

  beta <- as.numeric(summary(fit_lme_un)$coefficients$cond[, 1])
  X <- model.matrix(~ 0 + TRT + PERIOD, data = dfo)
  res_un <- summary(fit_lme_un)$vcov
  fit_un <- contGEEBCV(y=y,X=X,beta=beta,tau=wpicc_ar_un,rho=ar1_un,phi=phi,id=id,t=t,mv=mv)
  
  # fit model 2
  fit_lme_ca <- glmmTMB(Y ~  0 + TRT + PERIOD + X1 + X3 + ar1(PERIOD + 0 | CLUSTER), REML=T, 
                        data = dfo, family = gaussian)
  V_AR_ca <- VarCorr(fit_lme_ca)
  wpicc_ar_ca <- round(V_AR_ca$cond$CLUSTER[1,1]/sum(V_AR_ca$cond$CLUSTER[1,1],(attr(V_AR_ca$cond,"sc"))^2),8)
  ar1_ca <- round(attr(V_AR_ca$cond$CLUSTER,"correlation")[1,2],8)
  
  beta <- as.numeric(summary(fit_lme_ca)$coefficients$cond[, 1])
  X <- model.matrix(~ 0 + TRT + PERIOD + X1 + X3, data = dfo)
  res_ca <- summary(fit_lme_ca)$vcov
  fit_ca <- contGEEBCV(y=y,X=X,beta=beta,tau=wpicc_ar_ca,rho=ar1_ca,phi=phi,id=id,t=t,mv=mv)
  
  # fit model 3
  fit_lme_cac <- glmmTMB(Y ~  0 + TRT + PERIOD + X1 + X2 + X3 + X4 + ar1(PERIOD + 0 | CLUSTER), REML=T, 
                         data = dfo, family = gaussian)
  V_AR_cac <- VarCorr(fit_lme_cac)
  wpicc_ar_cac <- round(V_AR_cac$cond$CLUSTER[1,1]/sum(V_AR_cac$cond$CLUSTER[1,1],(attr(V_AR_cac$cond,"sc"))^2),8)
  ar1_cac <- round(attr(V_AR_cac$cond$CLUSTER,"correlation")[1,2],8)
  
  beta <- as.numeric(summary(fit_lme_cac)$coefficients$cond[, 1])
  X <- model.matrix(~ 0 + TRT + PERIOD + X1 + X2 + X3 + X4, data = dfo)
  res_cac <- summary(fit_lme_cac)$vcov
  fit_cac <- contGEEBCV(y=y,X=X,beta=beta,tau=wpicc_ar_cac,rho=ar1_cac,phi=phi,id=id,t=t,mv=mv)
  
  # estimates
  estimates <- c(summary(fit_lme_un)$coefficients$cond[1, 1], summary(fit_lme_ca)$coefficients$cond[1, 1], summary(fit_lme_cac)$coefficients$cond[1, 1])
  naive_se <-c(summary(fit_lme_un)$coefficients$cond[1, 2], summary(fit_lme_ca)$coefficients$cond[1, 2], summary(fit_lme_cac)$coefficients$cond[1, 2])
  
  naive_se2 <- sqrt(c(fit_un$naive[1, 1], fit_ca$naive[1, 1], fit_cac$naive[1, 1]))
  robust_se <- sqrt(c(fit_un$robust[1, 1], fit_ca$robust[1, 1], fit_cac$robust[1, 1]))
  KC_se <- sqrt(c(fit_un$varKC[1, 1], fit_ca$varKC[1, 1], fit_cac$varKC[1, 1]))
  MD_se <- sqrt(c(fit_un$varMD[1, 1], fit_ca$varMD[1, 1], fit_cac$varMD[1, 1]))
  
  wpicc_ars <- c(wpicc_ar_un, wpicc_ar_ca, wpicc_ar_cac)
  ar1s <- c(ar1_un, ar1_ca, ar1_cac)
  sigma2 <- c(summary(fit_lme_un)$sigma^2, summary(fit_lme_ca)$sigma^2, summary(fit_lme_cac)$sigma^2)
  
  # estimand
  yd <- df$Y1 - df$Y0
  estimand <- sum(yd)/sum(mvt)
  ydo <- dfo$Y1 - dfo$Y0
  estimando <- sum(ydo)/sum(mv)
  
  # results [(2+3*8+1+3)=30]
  return(c(estimand, estimando, estimates, naive_se, naive_se2, robust_se, KC_se, MD_se, wpicc_ars, ar1s, phi, sigma2))
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
  ###################################### exponential decay correlation ######################################
  data <- data_all
  
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
  tau <- colMeans(data[, 21:23])
  rho <- colMeans(data[, 24:26])
  phi <- mean(data[, 27])
  
  ase <- ec <- NULL
  for (i in 1:5){
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
  
  res_ed <- cbind(estimand, estimando, estimate, bias, relative_bias, ese, c(1, re_ca, re_cac), ase, ec, tau, rho, phi)
  
  # save results
  res <- res_ed
  colnames(res) <- c("ATE", "ATEo", "Estimate", "ABias", "ARBias", "EmpSE", "RE", 
                     "ASE_MB", "ASE_MB2", "ASE_ROB", "ASE_KC", "ASE_MD", 
                     "ECP_MB", "ECP_MB2", "ECP_ROB", "ECP_KC", "ECP_MD", 
                     "WP-ICC", "DC", "VAR(Y)")
  rownames(res) <- c("mixed-model unadjusted", "mixed-model adjusted (incomplete)", "mixed-model adjusted (complete)")
  
  return(res)
}

# ---------------------------------------------------------------------------------


# Scenario 1 (30 clusters)
n <- 30; t <- 6
pplv <- rep(1000, n)
cpl <- 5; cpu <- 50
beta0 <- 2; beta1 <- 1
mu_bx <- c(0.5, 0.8)
sd_cxc <- c(sqrt(0.1), sqrt(0.1)); sd_cxe <- c(sqrt(0.4), sqrt(0.9))
sdb <- sqrt(0.1); sdep <- sqrt(0.9)
nsim <- 1000

cres1_30 <- simulation(n, t, pplv, cpl, cpu, beta0, beta1, mu_bx, sd_cxc, sd_cxe, sdb, sdep, nsim)
save(cres1_30, file = "Data Results A/cres1ed_30.RData")
res1 <- simuRES(cres1_30, n, t)
res1

# Scenario 1 (100 clusters)
n <- 100; t <- 6
pplv <- rep(1000, n)
cpl <- 5; cpu <- 50
beta0 <- 2; beta1 <- 1
mu_bx <- c(0.5, 0.8)
sd_cxc <- c(sqrt(0.1), sqrt(0.1)); sd_cxe <- c(sqrt(0.4), sqrt(0.9))
sdb <- sqrt(0.1); sdep <- sqrt(0.9)
nsim <- 1000

cres1_100 <- simulation(n, t, pplv, cpl, cpu, beta0, beta1, mu_bx, sd_cxc, sd_cxe, sdb, sdep, nsim)
save(cres1_100, file = "Data Results A/cres1ed_100.RData")
res1 <- simuRES(cres1_100, n, t)
res1


# ---------------------------------------------------------------------------------


