#########################################
# binary outcome
# Estimand with data generation process 2
#########################################

sim_estimand <- function(n, t, pplv, cpl, cpu, beta0, beta1, mu_bx, sd_cxc, sd_cxe, sdb, period_effect, dg_model, e_scale){
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
    beta_i <- rnorm(1, 0, 0.2)
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
  
  # estimand ------
  ee_true <- map(complete_data, function(df_i){
    ee <- matrix(NA, 1, 5)
    i <- 0
    for (pi in 1:(t-1)){
      for (di in 0:pi){
        i <- i + 1
        ee[1, i] <- mean(df_i[,paste0("mu",di,"-period", pi)])
      }
    }
    ee
  })
  ee_emd <- colMeans(do.call(rbind, ee_true), na.rm = TRUE)
  
  names(ee_emd) <- c("period1-duration0", "period1-duration1",
                     "period2-duration0", "period2-duration1", "period2-duration2")
  
  return(ee_emd)
}


