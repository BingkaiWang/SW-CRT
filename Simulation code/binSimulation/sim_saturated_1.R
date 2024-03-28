#########################################
# binary outcome
# Analysis with data generation process 1
# (larger variance of beta_i for cluster-level intervention)
# without covariate adjustment
#########################################

sim_saturated <- function(n, t, pplv, cpl, cpu, beta0, beta1, mu_bx, sd_cxc, sd_cxe, sdb, period_effect, dg_model, e_scale){
  rm(list=ls())
  
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
        select(cluster_id:X2, period,trt_duration,Y) %>%
        .[sample(1:n_i, size = nj, replace = F),] 
    })
    rownames(dfo_i) <- NULL
    dfo_i
  })
  dfo <- map_dfr(observed_data, ~(.)) %>%
    mutate(trt_saturated = ifelse(trt_duration=="0", trt_duration, paste0("period",period,"duration",trt_duration)))# create treatment indicators
  
  ###############################################################################
  # fit model (nested exchangeable correlation, saturated treatment effect)------
  ###############################################################################
  
  fit_lme_ca <- lmer(Y ~ 0 +  period + trt_saturated + (1 | cluster_id) + (1 | cluster_period), 
                     data = dfo %>% mutate(cluster_period = paste0(cluster_id, period)))
  beta <- as.numeric(summary(fit_lme_ca)$coefficients[, 1])
  tau2 <- summary(fit_lme_ca)$varcor$cluster_id[1]
  kappa2 <-  summary(fit_lme_ca)$varcor$cluster_period[1]
  sigma2 <- summary(fit_lme_ca)$sigma^2
  QQ <-     model.matrix(~0 +  period + trt_saturated , data = dfo)
  Qo <- map(1:n, function(i) QQ[dfo$cluster_id==as.character(i),] )
  
  
  # influence function for model parameters -------
  sand_var_components <- map(1:n, function(i){
    dfo_i <- observed_data[[i]]
    m_i <- nrow(dfo_i)
    n_ij <- dfo_i %>% group_by(period) %>% summarise(N_ij = n()) %>% .$N_ij
    bdiag_mat <- bdiag(map(n_ij, ~matrix(1, nrow = ., ncol =.)))
    Sigma_i <- tau2 * matrix(1, nrow = m_i, ncol = m_i) + sigma2 * diag(m_i) + kappa2 * bdiag_mat
    Sigma_i_inv <- solve(Sigma_i)
    psi_i <- rbind(
      2 * t(Qo[[i]]) %*% Sigma_i_inv %*%  (dfo_i$Y - Qo[[i]] %*% beta),
      - sum(diag(Sigma_i_inv)) + sum((Sigma_i_inv %*%  (dfo_i$Y - Qo[[i]] %*% beta))^2),
      - sum(Sigma_i_inv) + (t(dfo_i$Y - Qo[[i]] %*% beta) %*% rowSums(Sigma_i_inv))^2,
      - sum(diag(Sigma_i_inv %*% bdiag_mat)) + t(dfo_i$Y - Qo[[i]] %*% beta) %*% Sigma_i_inv %*% bdiag_mat %*% Sigma_i_inv %*% (dfo_i$Y - Qo[[i]] %*% beta)
    )
    B11 <- -2 * t(Qo[[i]]) %*% Sigma_i_inv %*% Qo[[i]]
    B12 <- -2 * t(Qo[[i]]) %*% Sigma_i_inv %*% Sigma_i_inv %*% (dfo_i$Y - Qo[[i]] %*% beta)
    B13 <- -2 * t(Qo[[i]]) %*% rowSums(Sigma_i_inv) %*% colSums(Sigma_i_inv) %*% (dfo_i$Y - Qo[[i]] %*% beta)
    B14 <- -2 * t(Qo[[i]]) %*% Sigma_i_inv %*% bdiag_mat %*% Sigma_i_inv %*% (dfo_i$Y - Qo[[i]] %*% beta)
    B22 <- sum(diag(Sigma_i_inv %*% Sigma_i_inv)) - 
      2 * t(dfo_i$Y - Qo[[i]] %*% beta) %*% Sigma_i_inv %*% Sigma_i_inv %*% Sigma_i_inv %*% (dfo_i$Y - Qo[[i]] %*% beta)
    B23 <- sum(Sigma_i_inv %*% Sigma_i_inv) - 
      2 * t(dfo_i$Y - Qo[[i]] %*% beta) %*% rowSums(Sigma_i_inv) * t(dfo_i$Y - Qo[[i]] %*% beta) %*% rowSums(Sigma_i_inv %*% Sigma_i_inv)
    B24 <- sum(diag(Sigma_i_inv %*% Sigma_i_inv %*% bdiag_mat)) - 
      2 * t(dfo_i$Y - Qo[[i]] %*% beta) %*% Sigma_i_inv %*% bdiag_mat %*% Sigma_i_inv %*% Sigma_i_inv %*% (dfo_i$Y - Qo[[i]] %*% beta)
    B33 <- sum(Sigma_i_inv)^2 - 2 * sum(Sigma_i_inv) * (t(dfo_i$Y - Qo[[i]] %*% beta) %*% rowSums(Sigma_i_inv))^2
    B34 <- sum(Sigma_i_inv %*% bdiag_mat %*% Sigma_i_inv) -
      2 * t(dfo_i$Y - Qo[[i]] %*% beta) %*% Sigma_i_inv %*% bdiag_mat %*% rowSums(Sigma_i_inv) * colSums(Sigma_i_inv) %*% (dfo_i$Y - Qo[[i]] %*% beta)
    B44 <- sum(diag(Sigma_i_inv %*% bdiag_mat %*% Sigma_i_inv %*% bdiag_mat)) -
      2 * t(dfo_i$Y - Qo[[i]] %*% beta) %*% Sigma_i_inv %*% bdiag_mat %*% Sigma_i_inv %*% bdiag_mat %*% Sigma_i_inv %*% (dfo_i$Y - Qo[[i]] %*% beta) 
    B_i <- rbind(
      cbind(B11, B12, B13, B14),
      cbind(t(B12), B22, B23, B24),
      cbind(t(B13), t(B23), B33, B34),
      cbind(t(B14), t(B24), t(B34), B44)
    )
    list(psi_i, B_i)
  })
  B <- Reduce("+", map(sand_var_components, ~.[[2]])) / n
  IF <-  - solve(B) %*% Reduce(cbind, map(sand_var_components, ~.[[1]])) # influence functions of model-parameters
  
  # computing point estimates and IF of group means (variable: groupmean_estimates and IF_groupmean) -------
  IF_intercept <- IF[1:(t-1),]
  IF_trt <- IF[t:(t-1+t*(t-1)/2),]
  covariates_names <- c("X1", "X2")
  IF_betaX <- IF[(t-1+t*(t-1)/2) + 1:length(covariates_names),]
  mean_X <- t(map_dfr(observed_data, ~colMeans(.[,covariates_names])))
  
  beta_intercept <- beta[1:(t-1)]
  beta_trt <- beta[t:(t-1+t*(t-1)/2)]
  beta_X <- c(0, 0)
  
  groupmean_estimates <- rep(NA, t-1+t*(t-1)/2)
  for(j in 1:(t-1)){
    groupmean_estimates[j*(j-1)/2+j] <- beta_intercept[j] + mean(t(beta_X) %*% mean_X)
    names(groupmean_estimates)[j*(j-1)/2+j] <- paste0("period",j,"-duration",0)
    for(d in 1:j){
      groupmean_estimates[j*(j-1)/2+j+d] <- beta_intercept[j] + beta_trt[j*(j-1)/2+d] + mean(t(beta_X) %*% mean_X)
      names(groupmean_estimates)[j*(j-1)/2+j+d] <- paste0("period",j,"-duration",d)
    }
  }
  IF_groupmean <- matrix(NA, nrow = t-1+t*(t-1)/2, ncol = n)
  rownames(IF_groupmean) <- rep(NA, t-1+t*(t-1)/2)
  for(j in 1:(t-1)){
    IF_groupmean[j*(j-1)/2+j,] <- IF_intercept[j,] + as.vector(rowMeans(mean_X) %*% IF_betaX) +
      beta_intercept[j] + t(beta_X) %*% mean_X - groupmean_estimates[j*(j-1)/2+j]
    rownames(IF_groupmean)[j*(j-1)/2+j] <- paste0("period",j,"-duration",0)
    for(d in 1:j){
      IF_groupmean[j*(j-1)/2+j+d,] <- IF_intercept[j,] + IF_trt[j*(j-1)/2+d,] + as.vector(rowMeans(mean_X) %*% IF_betaX) +
        beta_intercept[j] + beta_trt[j*(j-1)/2+d] + t(beta_X) %*% mean_X - groupmean_estimates[j*(j-1)/2+j+d]
      rownames(IF_groupmean)[j*(j-1)/2+j+d] <- paste0("period",j,"-duration",d)
    }
  }
  
  if (e_scale == "rd"){
    # Example: computing risk difference: E[Y_ij(d)] - E[Y_ij(0)] for period j = 1,2 ------
    rd_est <- c(groupmean_estimates["period1-duration1"]-groupmean_estimates["period1-duration0"],
                groupmean_estimates["period2-duration1"]-groupmean_estimates["period2-duration0"],
                groupmean_estimates["period2-duration2"]-groupmean_estimates["period2-duration0"])
    rd_se_fun <- function(est, joint_IF){
      sqrt(c(1, -1) %*% var(joint_IF) %*% c(1, -1))
    }
    rd_se <- c(rd_se_fun(groupmean_estimates[c(2,1)], t(IF_groupmean[c(2,1),])),
               rd_se_fun(groupmean_estimates[c(4,3)], t(IF_groupmean[c(4,3),])),
               rd_se_fun(groupmean_estimates[c(5,3)], t(IF_groupmean[c(5,3),])))/sqrt(n)
    
    # res ------
    res <- as.numeric(c(rd_est, rd_se))
  } else if (e_scale == "rr"){
    # Example: computing risk ratio: E[Y_ij(d)]/E[Y_ij(0)] for period j = 1,2 ------
    rr_est <- c(groupmean_estimates["period1-duration1"]/groupmean_estimates["period1-duration0"],
                groupmean_estimates["period2-duration1"]/groupmean_estimates["period2-duration0"],
                groupmean_estimates["period2-duration2"]/groupmean_estimates["period2-duration0"])
    rr_se_fun <- function(est, joint_IF){
      sqrt(c(1/est[2], -est[1]/est[2]^2) %*% var(joint_IF) %*% c(1/est[2], -est[1]/est[2]^2))
    }
    rr_se <- c(rr_se_fun(groupmean_estimates[c(2,1)], t(IF_groupmean[c(2,1),])),
               rr_se_fun(groupmean_estimates[c(4,3)], t(IF_groupmean[c(4,3),])),
               rr_se_fun(groupmean_estimates[c(5,3)], t(IF_groupmean[c(5,3),])))/sqrt(n)
    
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
    or_se_fun <- function(est, joint_IF){
      sqrt(c((1-est[2])/(est[2]*(1-est[1])^2), -est[1]/((1-est[1])*est[2]^2)) %*% var(joint_IF) %*% c((1-est[2])/(est[2]*(1-est[1])^2), -est[1]/((1-est[1])*est[2]^2)))
    }
    or_se <- c(or_se_fun(groupmean_estimates[c(2,1)], t(IF_groupmean[c(2,1),])),
               or_se_fun(groupmean_estimates[c(4,3)], t(IF_groupmean[c(4,3),])),
               or_se_fun(groupmean_estimates[c(5,3)], t(IF_groupmean[c(5,3),])))/sqrt(n)
    
    # res ------
    res <- as.numeric(c(or_est, or_se))
  }
  
  ###############################################################################
  # fit model (gee with independence correlation structure)------
  ###############################################################################
  
  fit_gee <- glm(Y~0+period + trt_saturated, data = dfo, family = "binomial")
  QQ <- model.matrix(~0 +  period + trt_saturated, data = dfo)
  Qo <- map(1:n, function(i) QQ[dfo$cluster_id==as.character(i),] )
  
  # influence function for group means -------
  mu_levels <- data.frame(period = c("1", "1", "2", "2", "2"),
                          trt_saturated = c("0", "period1duration1", "0", "period2duration1", "period2duration2"))
  mu_hat <- map_dbl(1:nrow(mu_levels), function(iter){
    mean(predict(fit_gee, type = "response", newdata = mutate(dfo, period = mu_levels$period[iter], trt_saturated = mu_levels$trt_saturated[iter])))
  })
  names(mu_hat) <- c("period1-duration0", "period1-duration1", "period2-duration0", "period2-duration1", "period2-duration2")
  
  sand_var_components <- map(1:n, function(i){
    dfo_i <- filter(dfo, cluster_id == as.character(i))
    m_i <- nrow(dfo_i)
    mu_hat_i <- map(1:nrow(mu_levels), function(iter){
      predict(fit_gee, type = "response", newdata = mutate(dfo_i, period = mu_levels$period[iter], trt_saturated = mu_levels$trt_saturated[iter]))
    })
    psi_1 <- m_i * mu_hat - map_dbl(mu_hat_i, sum)
    pred_i <- predict(fit_gee,  newdata = dfo_i,type = "response")
    psi_2 <- 2 * t(Qo[[i]]) %*% (dfo_i$Y - pred_i)
    B11 <- diag(m_i, nrow = length(psi_1), ncol = length(psi_1))
    B22 <- -2 * t(Qo[[i]]) %*% (((pred_i * (1-pred_i)) %*% t(rep(1, ncol(Qo[[i]])))) * Qo[[i]])
    B21 <- matrix(0, ncol = ncol(B11), nrow = nrow(B22))
    B12 <- map_dfr(1:nrow(mu_levels), function(iter){
      Qm_i <- Qo[[i]]
      Qm_i[,1:5] <- 0
      Qm_i[,paste0("period", mu_levels$period[iter])] <- 1
      if( mu_levels$trt_saturated[iter] != "0"){
        Qm_i[,paste0("trt_saturated", mu_levels$trt_saturated[iter])] <- 1
      }
      - colSums(((mu_hat_i[[iter]] * (1 - mu_hat_i[[iter]])) %*% t(rep(1, ncol(Qm_i)))) * Qm_i)
    }) %>% as.matrix
    B_i <- rbind(
      cbind(B11, B12),
      cbind(B21, B22)
    )
    list(c(psi_1, psi_2), B_i)
  })
  B <- Reduce("+", map(sand_var_components, ~.[[2]])) / n
  IF <-  - solve(B) %*% Reduce(cbind, map(sand_var_components, ~.[[1]])) # influence functions of all model-parameters
  IF_groupmean <- IF[1:5,]
  groupmean_estimates <- mu_hat
  
  if (e_scale == "rd"){
    # Example: computing risk difference: E[Y_ij(d)] - E[Y_ij(0)] for period j = 1,2 ------
    rd_est <- c(groupmean_estimates["period1-duration1"]-groupmean_estimates["period1-duration0"],
                groupmean_estimates["period2-duration1"]-groupmean_estimates["period2-duration0"],
                groupmean_estimates["period2-duration2"]-groupmean_estimates["period2-duration0"])
    rd_se_fun <- function(est, joint_IF){
      sqrt(c(1, -1) %*% var(joint_IF) %*% c(1, -1))
    }
    rd_se <- c(rd_se_fun(groupmean_estimates[c(2,1)], t(IF_groupmean[c(2,1),])),
               rd_se_fun(groupmean_estimates[c(4,3)], t(IF_groupmean[c(4,3),])),
               rd_se_fun(groupmean_estimates[c(5,3)], t(IF_groupmean[c(5,3),])))/sqrt(n)
    
    # res ------
    res <- as.numeric(c(res, rd_est, rd_se))
  } else if (e_scale == "rr"){
    # Example: computing risk ratio: E[Y_ij(d)]/E[Y_ij(0)] for period j = 1,2 ------
    rr_est <- c(groupmean_estimates["period1-duration1"]/groupmean_estimates["period1-duration0"],
                groupmean_estimates["period2-duration1"]/groupmean_estimates["period2-duration0"],
                groupmean_estimates["period2-duration2"]/groupmean_estimates["period2-duration0"])
    rr_se_fun <- function(est, joint_IF){
      sqrt(c(1/est[2], -est[1]/est[2]^2) %*% var(joint_IF) %*% c(1/est[2], -est[1]/est[2]^2))
    }
    rr_se <- c(rr_se_fun(groupmean_estimates[c(2,1)], t(IF_groupmean[c(2,1),])),
               rr_se_fun(groupmean_estimates[c(4,3)], t(IF_groupmean[c(4,3),])),
               rr_se_fun(groupmean_estimates[c(5,3)], t(IF_groupmean[c(5,3),])))/sqrt(n)
    
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
    or_se_fun <- function(est, joint_IF){
      sqrt(c((1-est[2])/(est[2]*(1-est[1])^2), -est[1]/((1-est[1])*est[2]^2)) %*% var(joint_IF) %*% c((1-est[2])/(est[2]*(1-est[1])^2), -est[1]/((1-est[1])*est[2]^2)))
    }
    or_se <- c(or_se_fun(groupmean_estimates[c(2,1)], t(IF_groupmean[c(2,1),])),
               or_se_fun(groupmean_estimates[c(4,3)], t(IF_groupmean[c(4,3),])),
               or_se_fun(groupmean_estimates[c(5,3)], t(IF_groupmean[c(5,3),])))/sqrt(n)
    
    ###############################################################################
    # fit model using glmer (nested exchangeable correlation, saturated treatment effect)------
    ###############################################################################
    
    fit_glm_ca <- glmer(Y ~ 0 +  period + trt_saturated + (1 | cluster_id) + (1 | cluster_period), family = binomial, 
                        data = dfo %>% mutate(cluster_period = paste0(cluster_id, period)))
    beta_glm <- as.numeric(summary(fit_glm_ca)$coefficients[, 1])[3:5]
    se_glm <- as.numeric(summary(fit_glm_ca)$coefficients[, 2])[3:5]
    
    or_est_glm <- exp(beta_glm)
    or_se_glm_fun <- function(est, se){
      exp(est) * se
    }
    or_se_glm <- c(or_se_glm_fun(beta_glm[1], se_glm[1]),
                   or_se_glm_fun(beta_glm[2], se_glm[2]),
                   or_se_glm_fun(beta_glm[3], se_glm[3]))
    
    ###############################################################################
    # fit model using gee------
    ###############################################################################
    
    fit_geeglm <- geeglm(Y ~ 0 + period + trt_saturated, data = dfo, id = cluster_id, family = "binomial", corstr = "exchangeable")
    beta_gee <- as.numeric(summary(fit_geeglm)$coefficients[, 1])[3:5]
    se_gee <- as.numeric(summary(fit_geeglm)$coefficients[, 2])[3:5]
    
    or_est_gee <- exp(beta_gee)
    or_se_gee_fun <- function(est, se){
      exp(est) * se
    }
    or_se_gee <- c(or_se_gee_fun(beta_gee[1], se_gee[1]),
                   or_se_gee_fun(beta_gee[2], se_gee[2]),
                   or_se_gee_fun(beta_gee[3], se_gee[3]))
    
    # res ------
    res <- as.numeric(c(res, or_est, or_se, or_est_glm, or_se_glm, or_est_gee, or_se_gee))
  }
  
  return(res)
}


