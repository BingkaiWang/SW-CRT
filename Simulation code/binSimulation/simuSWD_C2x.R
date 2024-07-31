#############################################################
# Simulation for SWD (with Data Generation Process C2)
# with binary outcomes
# with covariate adjustment

# March 2023

# INPUT
# n: Total number of clusters
# t: Period
# pplv: Vector of population size in each cluster
# cpl: Lower bound of the cluster-period sizes
# cpu: Upper bound of the cluster-period sizes
# beta0: Parameter of the constant intervention effect
# beta1: Parameter of the intervention effect related to covariates
# mu_bx: Mean vector of binary covariates
# sd_cxc: Sd vector of cluster-level random effect(s) for continuous covariates
# sd_cxe: Sd vector of random error(s) for continuous covariates
# sdb: Sd of the cluster randomized effects
# period_effect: Intercept for each period
#############################################################
library(tidyverse)
library(lme4)
library(geepack)

source("sim_saturated_2x.R")
source("sim_estimand_2.R")

simulation <- function(n, t, pplv, cpl, cpu, beta0, beta1, mu_bx, sd_cxc, sd_cxe, sdb, period_effect, dg_model, e_scale, nsim){
  results <- NULL
  nerror <- 0
  sj <- 0
  for (ss in 1:nsim){
    passMODEL <- 1
    
    while(passMODEL == 1){
      sj <- sj+1
      set.seed(sj)
      res_ss <- try(sim_saturated(n=n, t=t, pplv=pplv, cpl=cpl, cpu=cpu, beta0=beta0, beta1=beta1, mu_bx=mu_bx, 
                                  sd_cxc=sd_cxc, sd_cxe=sd_cxe, sdb=sdb, period_effect=period_effect,
                                  dg_model=dg_model, e_scale=e_scale), silent=TRUE)
      if(class(res_ss)[1]=="try-error" | is.null(res_ss)==1){passMODEL <- 1; nerror <- nerror + 1; next}
      passMODEL <- 0
    }
    
    results <- rbind(results, res_ss)
    
    # Loop control
    if(ss %% 50 == 0){print(paste("Iteration", ss, "NERROR =", nerror))}
  }
  return(results)
}

simulation_emd <- function(n, t, pplv, cpl, cpu, beta0, beta1, mu_bx, sd_cxc, sd_cxe, sdb, period_effect, dg_model, e_scale, nsim){
  results <- NULL
  nerror <- 0
  sj <- 0
  for (ss in 1:nsim){
    passMODEL <- 1
    
    while(passMODEL == 1){
      sj <- sj+1
      set.seed(sj)
      res_ss <- try(sim_estimand(n=n, t=t, pplv=pplv, cpl=cpl, cpu=cpu, beta0=beta0, beta1=beta1, mu_bx=mu_bx, 
                                 sd_cxc=sd_cxc, sd_cxe=sd_cxe, sdb=sdb, period_effect=period_effect,
                                 dg_model=dg_model, e_scale=e_scale), silent=TRUE)
      if(class(res_ss)[1]=="try-error" | is.null(res_ss)==1){passMODEL <- 1; nerror <- nerror + 1; next}
      passMODEL <- 0
    }
    
    results <- rbind(results, res_ss)
    
    # Loop control
    if(ss %% 50 == 0){print(paste("Iteration", ss, "NERROR =", nerror))}
  }
  
  ee_emd <- colMeans(results)
  
  if (e_scale == "rd"){
    emd <- c(ee_emd["period1-duration1"]-ee_emd["period1-duration0"],
             ee_emd["period2-duration1"]-ee_emd["period2-duration0"],
             ee_emd["period2-duration2"]-ee_emd["period2-duration0"])
  } else if (e_scale == "rr"){
    emd <- c(ee_emd["period1-duration1"]/ee_emd["period1-duration0"],
             ee_emd["period2-duration1"]/ee_emd["period2-duration0"],
             ee_emd["period2-duration2"]/ee_emd["period2-duration0"])
  } else if (e_scale == "or"){
    or_emd_fun <- function(emd1, emd0){
      (emd1/(1-emd1))/(emd0/(1-emd0))
    }
    emd <- c(or_emd_fun(ee_emd["period1-duration1"], ee_emd["period1-duration0"]),
             or_emd_fun(ee_emd["period2-duration1"], ee_emd["period2-duration0"]),
             or_emd_fun(ee_emd["period2-duration2"], ee_emd["period2-duration0"]))
  }
  
  return(emd)
}

# ---------------------------------------------------------------------------------


# Scenario 2

############
# rd
############

t <- 3
ppl <- 5000
cpl <- 50; cpu <- 100
beta0 <- 2; beta1 <- 1
mu_bx <- 0.1
sd_cxc <- sqrt(0.1); sd_cxe <- sqrt(0.4)
sdb <- c(sqrt(0.06), sqrt(0.04))
period_effect <- seq(0.1, 0.12, length.out=t)
dg_model <- "logit"
e_scale <- "rd"
nsim <- 1000

n <- 1000
pplv <- rep(ppl, n)
bres2_rd_emd <- simulation_emd(n, t, pplv, cpl, cpu, beta0, beta1, mu_bx, sd_cxc, sd_cxe, sdb, period_effect, dg_model, e_scale, nsim)

nsim <- 1000
bres2_rd_emd <- matrix(bres2_rd_emd, 1, 3)[rep(1, nsim), ]

# 30 clusters
n <- 30
pplv <- rep(ppl, n)
bres2_rd_30 <- simulation(n, t, pplv, cpl, cpu, beta0, beta1, mu_bx, sd_cxc, sd_cxe, sdb, period_effect, dg_model, e_scale, nsim)
bres2_rd_30 <- cbind(bres2_rd_emd, bres2_rd_30)
save(bres2_rd_30, file = "Data Results B/bres2_rd_30.RData")

# 100 clusters
n <- 100
pplv <- rep(ppl, n)
bres2_rd_100 <- simulation(n, t, pplv, cpl, cpu, beta0, beta1, mu_bx, sd_cxc, sd_cxe, sdb, period_effect, dg_model, e_scale, nsim)
bres2_rd_100 <- cbind(bres2_rd_emd, bres2_rd_100)
save(bres2_rd_100, file = "Data Results B/bres2_rd_100.RData")


############
# rr
############

t <- 3
ppl <- 5000
cpl <- 50; cpu <- 100
beta0 <- 2; beta1 <- 1
mu_bx <- 0.1
sd_cxc <- sqrt(0.1); sd_cxe <- sqrt(0.4)
sdb <- c(sqrt(0.06), sqrt(0.04))
period_effect <- seq(0.1, 0.12, length.out=t)
dg_model <- "logit"
e_scale <- "rr"
nsim <- 1000

n <- 1000
pplv <- rep(ppl, n)
bres2_rr_emd <- simulation_emd(n, t, pplv, cpl, cpu, beta0, beta1, mu_bx, sd_cxc, sd_cxe, sdb, period_effect, dg_model, e_scale, nsim)

nsim <- 1000
bres2_rr_emd <- matrix(bres2_rr_emd, 1, 3)[rep(1, nsim), ]

# 30 clusters
n <- 30
pplv <- rep(ppl, n)
bres2_rr_30 <- simulation(n, t, pplv, cpl, cpu, beta0, beta1, mu_bx, sd_cxc, sd_cxe, sdb, period_effect, dg_model, e_scale, nsim)
bres2_rr_30 <- cbind(bres2_rr_emd, bres2_rr_30)
save(bres2_rr_30, file = "Data Results B/bres2_rr_30.RData")

# 100 clusters
n <- 100
pplv <- rep(ppl, n)
bres2_rr_100 <- simulation(n, t, pplv, cpl, cpu, beta0, beta1, mu_bx, sd_cxc, sd_cxe, sdb, period_effect, dg_model, e_scale, nsim)
bres2_rr_100 <- cbind(bres2_rr_emd, bres2_rr_100)
save(bres2_rr_100, file = "Data Results B/bres2_rr_100.RData")


############
# or
############

t <- 3
ppl <- 5000
cpl <- 50; cpu <- 100
beta0 <- 0.9; beta1 <- 0.5
mu_bx <- 0.1
sd_cxc <- sqrt(0.1); sd_cxe <- sqrt(0.4)
sdb <- c(sqrt(0.06), sqrt(0.04))
period_effect <- seq(0.1, 0.12, length.out=t)
dg_model <- "logit"
e_scale <- "or"
nsim <- 1000

n <- 1000
pplv <- rep(ppl, n)
bres2_or_emd <- simulation_emd(n, t, pplv, cpl, cpu, beta0, beta1, mu_bx, sd_cxc, sd_cxe, sdb, period_effect, dg_model, e_scale, nsim)

nsim <- 1000
bres2_or_emd <- matrix(bres2_or_emd, 1, 3)[rep(1, nsim), ]

# 30 clusters
n <- 30
pplv <- rep(ppl, n)
bres2_or_30 <- simulation(n, t, pplv, cpl, cpu, beta0, beta1, mu_bx, sd_cxc, sd_cxe, sdb, period_effect, dg_model, e_scale, nsim)
bres2_or_30 <- cbind(bres2_or_emd, bres2_or_30)
save(bres2_or_30, file = "Data Results B/bres2_or_30.RData")

# 100 clusters
n <- 100
pplv <- rep(ppl, n)
bres2_or_100 <- simulation(n, t, pplv, cpl, cpu, beta0, beta1, mu_bx, sd_cxc, sd_cxe, sdb, period_effect, dg_model, e_scale, nsim)
bres2_or_100 <- cbind(bres2_or_emd, bres2_or_100)
save(bres2_or_100, file = "Data Results B/bres2_or_100.RData")

# 300 clusters
n <- 300
pplv <- rep(ppl, n)
bres2_or_300 <- simulation(n, t, pplv, cpl, cpu, beta0, beta1, mu_bx, sd_cxc, sd_cxe, sdb, period_effect, dg_model, e_scale, nsim)
bres2_or_300 <- cbind(bres2_or_emd, bres2_or_300)
save(bres2_or_300, file = "Data Results B/bres2_or_300.RData")

# 1000 clusters
n <- 1000
pplv <- rep(ppl, n)
bres2_or_1000 <- simulation(n, t, pplv, cpl, cpu, beta0, beta1, mu_bx, sd_cxc, sd_cxe, sdb, period_effect, dg_model, e_scale, nsim)
bres2_or_1000 <- cbind(bres2_or_emd, bres2_or_1000)
save(bres2_or_1000, file = "Data Results B/bres2_or_1000.RData")


