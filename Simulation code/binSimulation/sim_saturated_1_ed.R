#########################################
# binary outcome
# Analysis with data generation process 1
# (larger variance of beta_i for cluster-level intervention)
# without covariate adjustment
#########################################

sim_saturated <- function(n, t, pplv, cpl, cpu, beta0, beta1, mu_bx, sd_cxc, sd_cxe, sdb, period_effect, dg_model, e_scale){
  rm(list=ls())
  
  source("contGEE_BCV2_ED.R")
  
  # expit function ------
  expit <- function(x){y <- exp(x)/(1+exp(x));return(y)}
  
  # generating complete data with all potential outcomes -------
  complete_data <- map(1:n, function(i){
    df_i <- data.frame(cluster_id = rep(paste0(i), pplv[i]),
                       trt_seq = rep(sample(1:t, 1), pplv[i]),
                       X1 = rbinom(pplv[i], size = 1, prob=mu_bx),
                       X2 = rnorm(1, 0, sd_cxc) + rnorm(pplv[i], 0, sd_cxe)
    )
    alpha_i <- rnorm(1, 0, sdb[1])
    theta_ij <- rnorm(3, 0, sdb[2])
    beta_i <- rnorm(1, 0, 0.5)
    Yd <- mud <- vector("list", t+1)
    for(j in 1:(t+1)){
      Yd[[j]] <- mud[[j]] <- matrix(NA, nrow = pplv[i], ncol = t)
      colnames(Yd[[j]]) <- paste0("Y", j-1, "-period", 1:t)
      colnames(mud[[j]]) <- paste0("mu", j-1, "-period", 1:t)
      for(jj in 1:t){
        mud[[j]][,jj] <- period_effect[jj] + alpha_i + theta_ij[jj] + 
          df_i$X1 + jj/t*df_i$X2^2 + 
          ifelse(j==1, 0, beta0*(0.8+0.2*(j-1)) + beta1/2*(df_i$X1-mean(df_i$X1)) + beta1*0.2*(j-1)*(df_i$X2^2-mean(df_i$X2^2)) + beta_i)
        if (dg_model == "logit"){
          # logit model
          mud[[j]][,jj] <- expit(mud[[j]][,jj])
        } else if (dg_model == "probit"){
          # probit model
          mud[[j]][,jj] <- pnorm(mud[[j]][,jj])
        }
        Yd[[j]][,jj] <- rbinom(pplv[i], 1, mud[[j]][,jj])
      }
    }
    cbind(df_i, do.call("cbind",Yd), do.call("cbind",mud))
  })
  
  # generating all observed_data (the last period is dropped) --------
  observed_data <- map(complete_data, function(df_i){
    n_i <- nrow(df_i)
    dfo_i <- map_dfr(1:(t-1), function(j){
      nj <-sample(cpl:cpu, size = 1)
      trt_duration <- ifelse(df_i$trt_seq[1] <= j, j-df_i$trt_seq[1]+1, 0)
      Y <- df_i[,paste0("Y",trt_duration,"-period", j)]
      df_i %>%
        mutate(period = as.character(j)) %>%
        mutate(trt_duration = as.character(trt_duration)) %>%
        mutate(Y = Y) %>%
        dplyr::select(cluster_id:X2, period,trt_duration,Y) %>%
        .[sample(1:n_i, size = nj, replace = F),] 
    })
    rownames(dfo_i) <- NULL
    dfo_i
  })
  dfo <- map_dfr(observed_data, ~(.)) %>%
    mutate(trt_saturated = ifelse(trt_duration=="0", trt_duration, paste0("period",period,"duration",trt_duration))) %>%# create treatment indicators
    mutate(cluster_id = as.numeric(cluster_id))
  
  dfo_cp <- dfo %>%
    group_by(cluster_id, as.numeric(period)) %>%
    summarise(mv = n(), .groups = 'drop')
  mv <- dfo_cp$mv
  
  ###############################################################################
  # fit model 1 (exponential decay correlation, saturated treatment effect)------
  ###############################################################################
  
  fit_lme_un <- glmmTMB(Y ~  0 + period + trt_saturated + ar1(period + 0 | cluster_id), REML=T, 
                        data = dfo, family = gaussian)
  V_AR_un <- VarCorr(fit_lme_un)
  wpicc_ar_un <- round(V_AR_un$cond$cluster_id[1,1]/sum(V_AR_un$cond$cluster_id[1,1],(attr(V_AR_un$cond,"sc"))^2),8)
  ar1_un <- round(attr(V_AR_un$cond$cluster_id,"correlation")[1,2],8)
  
  beta <- as.numeric(summary(fit_lme_un)$coefficients$cond[, 1])
  X <- model.matrix(~ 0 + period + trt_saturated, data = dfo)
  phi <- summary(fit_lme_un)$sigma^2
  res_un <- summary(fit_lme_un)$vcov
  fit_un <- contGEEBCV(y=dfo$Y,X=X,beta=beta,tau=wpicc_ar_un,rho=ar1_un,phi=phi,id=dfo$cluster_id,t=2,mv=mv)
  rob <- matrix(0, nrow = t-1+t*(t-1)/2 + 2, ncol = t-1+t*(t-1)/2 + 2)
  rob[1:(t-1+t*(t-1)/2), 1:(t-1+t*(t-1)/2)] <- fit_un$robust
  
  beta_intercept <- beta[1:(t-1)]
  beta_trt <- beta[t:(t-1+t*(t-1)/2)]
  beta_X <- c(0, 0)
  
  covariates_names <- c("X1", "X2")
  mean_X <- t(map_dfr(observed_data, ~colMeans(.[,covariates_names])))
  
  groupmean_estimates <- rep(NA, t-1+t*(t-1)/2)
  for(j in 1:(t-1)){
    groupmean_estimates[j*(j-1)/2+j] <- beta_intercept[j] + mean(t(beta_X) %*% mean_X)
    names(groupmean_estimates)[j*(j-1)/2+j] <- paste0("period",j,"-duration",0)
    for(d in 1:j){
      groupmean_estimates[j*(j-1)/2+j+d] <- beta_intercept[j] + beta_trt[j*(j-1)/2+d] + mean(t(beta_X) %*% mean_X)
      names(groupmean_estimates)[j*(j-1)/2+j+d] <- paste0("period",j,"-duration",d)
    }
  }
  rob_groupmean <- matrix(NA, nrow = t-1+t*(t-1)/2, ncol = t-1+t*(t-1)/2)
  rownames(rob_groupmean) <- rep(NA, t-1+t*(t-1)/2)
  for(j in 1:(t-1)){
    rob_groupmean[j*(j-1)/2+j, j*(j-1)/2+j] <- t(c(1, rowMeans(mean_X))) %*% rob[c(j, t-1+t*(t-1)/2+1, t-1+t*(t-1)/2+2), c(j, t-1+t*(t-1)/2+1, t-1+t*(t-1)/2+2)] %*% c(1, rowMeans(mean_X))
    rownames(rob_groupmean)[j*(j-1)/2+j] <- paste0("period",j,"-duration",0)
    for(d in 1:j){
      rob_groupmean[j*(j-1)/2+j+d, j*(j-1)/2+j+d] <- t(c(1, 1, rowMeans(mean_X))) %*% rob[c(j, t-1+j-1+d, t-1+t*(t-1)/2+1, t-1+t*(t-1)/2+2), c(j, t-1+j-1+d, t-1+t*(t-1)/2+1, t-1+t*(t-1)/2+2)] %*% c(1, 1, rowMeans(mean_X))
      rob_groupmean[j*(j-1)/2+j, j*(j-1)/2+j+d] <- t(c(1, 0, rowMeans(mean_X))) %*% rob[c(j, t-1+j-1+d, t-1+t*(t-1)/2+1, t-1+t*(t-1)/2+2), c(j, t-1+j-1+d, t-1+t*(t-1)/2+1, t-1+t*(t-1)/2+2)] %*% c(1, 1, rowMeans(mean_X))
      rownames(rob_groupmean)[j*(j-1)/2+j+d] <- paste0("period",j,"-duration",d)
    }
  }
  rob_groupmean[lower.tri(rob_groupmean)] = t(rob_groupmean)[lower.tri(rob_groupmean)]
  
  if (e_scale == "rd"){
    # Example: computing risk difference: E[Y_ij(d)] - E[Y_ij(0)] for period j = 1,2 ------
    rd_est <- c(groupmean_estimates["period1-duration1"]-groupmean_estimates["period1-duration0"],
                groupmean_estimates["period2-duration1"]-groupmean_estimates["period2-duration0"],
                groupmean_estimates["period2-duration2"]-groupmean_estimates["period2-duration0"])
    rd_se_fun <- function(est, variance){
      sqrt(c(1, -1) %*% variance %*% c(1, -1))
    }
    rd_se <- c(rd_se_fun(groupmean_estimates[c(2,1)], rob_groupmean[c(2,1),c(2,1)]),
               rd_se_fun(groupmean_estimates[c(4,3)], rob_groupmean[c(4,3),c(4,3)]),
               rd_se_fun(groupmean_estimates[c(5,3)], rob_groupmean[c(5,3),c(5,3)]))
    
    # res ------
    res <- as.numeric(c(rd_est, rd_se))
  } else if (e_scale == "rr"){
    # Example: computing risk ratio: E[Y_ij(d)]/E[Y_ij(0)] for period j = 1,2 ------
    rr_est <- c(groupmean_estimates["period1-duration1"]/groupmean_estimates["period1-duration0"],
                groupmean_estimates["period2-duration1"]/groupmean_estimates["period2-duration0"],
                groupmean_estimates["period2-duration2"]/groupmean_estimates["period2-duration0"])
    rr_se_fun <- function(est, variance){
      sqrt(c(1/est[2], -est[1]/est[2]^2) %*% variance %*% c(1/est[2], -est[1]/est[2]^2))
    }
    rr_se <- c(rr_se_fun(groupmean_estimates[c(2,1)], rob_groupmean[c(2,1),c(2,1)]),
               rr_se_fun(groupmean_estimates[c(4,3)], rob_groupmean[c(4,3),c(4,3)]),
               rr_se_fun(groupmean_estimates[c(5,3)], rob_groupmean[c(5,3),c(5,3)]))
    
    # res ------
    res <- as.numeric(c(rr_est, rr_se))
  } else if (e_scale == "or"){
    # Example: computing odds ratio: E[Y_ij(d)] to E[Y_ij(0)] for period j = 1,2 ------
    or_est_fun <- function(est1, est0){
      (est1/(1-est1))/(est0/(1-est0))
    }
    or_est <- c(or_est_fun(groupmean_estimates["period1-duration1"], groupmean_estimates["period1-duration0"]),
                or_est_fun(groupmean_estimates["period2-duration1"], groupmean_estimates["period2-duration0"]),
                or_est_fun(groupmean_estimates["period2-duration2"], groupmean_estimates["period2-duration0"]))
    or_se_fun <- function(est, variance){
      sqrt(c((1-est[2])/(est[2]*(1-est[1])^2), -est[1]/((1-est[1])*est[2]^2)) %*% variance %*% c((1-est[2])/(est[2]*(1-est[1])^2), -est[1]/((1-est[1])*est[2]^2)))
    }
    or_se <- c(or_se_fun(groupmean_estimates[c(2,1)], rob_groupmean[c(2,1),c(2,1)]),
               or_se_fun(groupmean_estimates[c(4,3)], rob_groupmean[c(4,3),c(4,3)]),
               or_se_fun(groupmean_estimates[c(5,3)], rob_groupmean[c(5,3),c(5,3)]))
    
    # res ------
    res <- as.numeric(c(or_est, or_se))
  }
  
  ###############################################################################
  # fit model 2 (exponential decay correlation, saturated treatment effect)------
  ###############################################################################
  
  fit_lme_ca <- glmmTMB(Y ~  0 + period + trt_saturated + X1 + X2 + ar1(period + 0 | cluster_id), REML=T, 
                        data = dfo, family = gaussian)
  V_AR_ca <- VarCorr(fit_lme_ca)
  wpicc_ar_ca <- round(V_AR_ca$cond$cluster_id[1,1]/sum(V_AR_ca$cond$cluster_id[1,1],(attr(V_AR_ca$cond,"sc"))^2),8)
  ar1_ca <- round(attr(V_AR_ca$cond$cluster_id,"correlation")[1,2],8)
  
  beta <- as.numeric(summary(fit_lme_ca)$coefficients$cond[, 1])
  X <- model.matrix(~ 0 + period + trt_saturated + X1 + X2, data = dfo)
  phi <- summary(fit_lme_ca)$sigma^2
  res_ca <- summary(fit_lme_ca)$vcov
  fit_ca <- contGEEBCV(y=dfo$Y,X=X,beta=beta,tau=wpicc_ar_ca,rho=ar1_ca,phi=phi,id=dfo$cluster_id,t=2,mv=mv)
  rob <- fit_ca$robust
  
  beta_intercept <- beta[1:(t-1)]
  beta_trt <- beta[t:(t-1+t*(t-1)/2)]
  beta_X <- beta[(t-1+t*(t-1)/2) + 1:length(covariates_names)]
  
  covariates_names <- c("X1", "X2")
  mean_X <- t(map_dfr(observed_data, ~colMeans(.[,covariates_names])))
  
  groupmean_estimates <- rep(NA, t-1+t*(t-1)/2)
  for(j in 1:(t-1)){
    groupmean_estimates[j*(j-1)/2+j] <- beta_intercept[j] + mean(t(beta_X) %*% mean_X)
    names(groupmean_estimates)[j*(j-1)/2+j] <- paste0("period",j,"-duration",0)
    for(d in 1:j){
      groupmean_estimates[j*(j-1)/2+j+d] <- beta_intercept[j] + beta_trt[j*(j-1)/2+d] + mean(t(beta_X) %*% mean_X)
      names(groupmean_estimates)[j*(j-1)/2+j+d] <- paste0("period",j,"-duration",d)
    }
  }
  rob_groupmean <- matrix(NA, nrow = t-1+t*(t-1)/2, ncol = t-1+t*(t-1)/2)
  rownames(rob_groupmean) <- rep(NA, t-1+t*(t-1)/2)
  for(j in 1:(t-1)){
    rob_groupmean[j*(j-1)/2+j, j*(j-1)/2+j] <- t(c(1, rowMeans(mean_X))) %*% rob[c(j, t-1+t*(t-1)/2+1, t-1+t*(t-1)/2+2), c(j, t-1+t*(t-1)/2+1, t-1+t*(t-1)/2+2)] %*% c(1, rowMeans(mean_X))
    rownames(rob_groupmean)[j*(j-1)/2+j] <- paste0("period",j,"-duration",0)
    for(d in 1:j){
      rob_groupmean[j*(j-1)/2+j+d, j*(j-1)/2+j+d] <- t(c(1, 1, rowMeans(mean_X))) %*% rob[c(j, t-1+j-1+d, t-1+t*(t-1)/2+1, t-1+t*(t-1)/2+2), c(j, t-1+j-1+d, t-1+t*(t-1)/2+1, t-1+t*(t-1)/2+2)] %*% c(1, 1, rowMeans(mean_X))
      rob_groupmean[j*(j-1)/2+j, j*(j-1)/2+j+d] <- t(c(1, 0, rowMeans(mean_X))) %*% rob[c(j, t-1+j-1+d, t-1+t*(t-1)/2+1, t-1+t*(t-1)/2+2), c(j, t-1+j-1+d, t-1+t*(t-1)/2+1, t-1+t*(t-1)/2+2)] %*% c(1, 1, rowMeans(mean_X))
      rownames(rob_groupmean)[j*(j-1)/2+j+d] <- paste0("period",j,"-duration",d)
    }
  }
  rob_groupmean[lower.tri(rob_groupmean)] = t(rob_groupmean)[lower.tri(rob_groupmean)]
  
  if (e_scale == "rd"){
    # Example: computing risk difference: E[Y_ij(d)] - E[Y_ij(0)] for period j = 1,2 ------
    rd_est <- c(groupmean_estimates["period1-duration1"]-groupmean_estimates["period1-duration0"],
                groupmean_estimates["period2-duration1"]-groupmean_estimates["period2-duration0"],
                groupmean_estimates["period2-duration2"]-groupmean_estimates["period2-duration0"])
    rd_se_fun <- function(est, variance){
      sqrt(c(1, -1) %*% variance %*% c(1, -1))
    }
    rd_se <- c(rd_se_fun(groupmean_estimates[c(2,1)], rob_groupmean[c(2,1),c(2,1)]),
               rd_se_fun(groupmean_estimates[c(4,3)], rob_groupmean[c(4,3),c(4,3)]),
               rd_se_fun(groupmean_estimates[c(5,3)], rob_groupmean[c(5,3),c(5,3)]))
    
    # res ------
    res <- as.numeric(c(res, rd_est, rd_se))
  } else if (e_scale == "rr"){
    # Example: computing risk ratio: E[Y_ij(d)]/E[Y_ij(0)] for period j = 1,2 ------
    rr_est <- c(groupmean_estimates["period1-duration1"]/groupmean_estimates["period1-duration0"],
                groupmean_estimates["period2-duration1"]/groupmean_estimates["period2-duration0"],
                groupmean_estimates["period2-duration2"]/groupmean_estimates["period2-duration0"])
    rr_se_fun <- function(est, variance){
      sqrt(c(1/est[2], -est[1]/est[2]^2) %*% variance %*% c(1/est[2], -est[1]/est[2]^2))
    }
    rr_se <- c(rr_se_fun(groupmean_estimates[c(2,1)], rob_groupmean[c(2,1),c(2,1)]),
               rr_se_fun(groupmean_estimates[c(4,3)], rob_groupmean[c(4,3),c(4,3)]),
               rr_se_fun(groupmean_estimates[c(5,3)], rob_groupmean[c(5,3),c(5,3)]))
    
    # res ------
    res <- as.numeric(c(res, rr_est, rr_se))
  } else if (e_scale == "or"){
    # Example: computing odds ratio: E[Y_ij(d)] to E[Y_ij(0)] for period j = 1,2 ------
    or_est_fun <- function(est1, est0){
      (est1/(1-est1))/(est0/(1-est0))
    }
    or_est <- c(or_est_fun(groupmean_estimates["period1-duration1"], groupmean_estimates["period1-duration0"]),
                or_est_fun(groupmean_estimates["period2-duration1"], groupmean_estimates["period2-duration0"]),
                or_est_fun(groupmean_estimates["period2-duration2"], groupmean_estimates["period2-duration0"]))
    or_se_fun <- function(est, variance){
      sqrt(c((1-est[2])/(est[2]*(1-est[1])^2), -est[1]/((1-est[1])*est[2]^2)) %*% variance %*% c((1-est[2])/(est[2]*(1-est[1])^2), -est[1]/((1-est[1])*est[2]^2)))
    }
    or_se <- c(or_se_fun(groupmean_estimates[c(2,1)], rob_groupmean[c(2,1),c(2,1)]),
               or_se_fun(groupmean_estimates[c(4,3)], rob_groupmean[c(4,3),c(4,3)]),
               or_se_fun(groupmean_estimates[c(5,3)], rob_groupmean[c(5,3),c(5,3)]))
    
    ###############################################################################
    # fit model 3 (exponential decay correlation, saturated treatment effect)------
    ###############################################################################
    
    fit_glmm_un <- glmmTMB(Y ~  0 + period + trt_saturated + ar1(period + 0 | cluster_id), REML=T, 
                           data = dfo, family = binomial(link = logit))
    beta_glmm_un <- as.numeric(summary(fit_glmm_un)$coefficients$cond[, 1])[3:5]
    naive_se_glmm_un <- as.numeric(summary(fit_glmm_un)$coefficients$cond[, 2])[3:5]
    
    or_se_fun <- function(est, se){
      exp(est) * se
    }
    
    or_est_glmm_un <- exp(beta_glmm_un)
    or_naive_se_glmm_un <- c(or_se_fun(beta_glmm_un[1], naive_se_glmm_un[1]),
                             or_se_fun(beta_glmm_un[2], naive_se_glmm_un[2]),
                             or_se_fun(beta_glmm_un[3], naive_se_glmm_un[3]))
    
    ###############################################################################
    # fit model 4 (exponential decay correlation, saturated treatment effect)------
    ###############################################################################
    
    fit_glmm_ca <- glmmTMB(Y ~  0 + period + trt_saturated + X1 + X2 + ar1(period + 0 | cluster_id), REML=T, 
                           data = dfo, family = binomial(link = logit))
    beta_glmm_ca <- as.numeric(summary(fit_glmm_ca)$coefficients$cond[, 1])[3:5]
    naive_se_glmm_ca <- as.numeric(summary(fit_glmm_ca)$coefficients$cond[, 2])[3:5]
    
    or_est_glmm_ca <- exp(beta_glmm_ca)
    or_naive_se_glmm_ca <- c(or_se_fun(beta_glmm_ca[1], naive_se_glmm_ca[1]),
                             or_se_fun(beta_glmm_ca[2], naive_se_glmm_ca[2]),
                             or_se_fun(beta_glmm_ca[3], naive_se_glmm_ca[3]))
    
    # res ------
    res <- as.numeric(c(res, or_est, or_se, or_est_glmm_un, or_naive_se_glmm_un, or_est_glmm_ca, or_naive_se_glmm_ca))
  }
  
  return(res)
}


