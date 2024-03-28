rm(list=ls())
set.seed(123)
library(tidyverse)
library(lme4)
setwd("~/Dropbox (Penn)/Mixed model SW-CRT Data Examples/Application/Data 2/")
# source("../../Simulation code/conSimulation/SV.R")
# source("../../Simulation code/conSimulation/SV_NE.R")
# source("../../Simulation code/conSimulation/contGEE_BCV2.R")

source("SWD-mixed-model-sandwich-variance.R")

data2 <- readxl::read_xlsx("2. SMI_data_and_dictionary.xlsx")
data2 <- data2 %>%
  dplyr::select(cluster_id = phc_code, 
                period = PHASE, 
                trt_onset = block,
                trt_indi = TRT, 
                Y_continuous = SBP_cal, 
                Y_binary = primary_event, 
                BP_under_140_base, gender, Age, SBP_cal_base, PH_CVD_base, PCT_GOVT_80) %>%
  # mutate(Y_continuous = Y_continuous - SBP_cal_base) %>%
  filter(period > 1) %>%
  mutate(period = as.character(period - 1), 
         trt_indi = as.numeric(trt_indi == "Intervention"), 
         Y_binary = Y_binary == "Yes") %>%
  mutate(trt_duration = as.character(ifelse(as.numeric(period) >= trt_onset, as.numeric(period) - trt_onset + 1, 0))) %>%
  mutate(trt_period = ifelse(trt_indi, paste0("period", period, trt_indi), "0")) %>%
  mutate(trt_saturated = ifelse(trt_duration=="0", "0", paste0("period",period,"duration",trt_duration))) %>%
  mutate(cluster_period = paste0(cluster_id, period))
data2 <- filter(data2, period!=3)

n <- 18
n_period <- 2

## Estimate constant treatment effect --------------------
# unadjusted, independence correlation
fit.constant <- lm(Y_continuous ~ 0 + period + trt_indi, data = data2)
beta <- as.numeric(summary(fit.constant)$coefficients[, 1])
X <- model.matrix(~ 0 + period + trt_indi, data = data2)
unadj_i <- SV_SWD_i(y=data2$Y_continuous, X, id=data2$cluster_id, beta, n, n_period, trt_model = "constant")
est_unadj_i <- unadj_i[1]
sd_unadj_i <- unadj_i[2]

# unadjusted, exchangeable correlation
fit.constant <- lmer(Y_continuous ~ 0 + period + trt_indi + (1 | cluster_id), data = data2)
beta <- as.numeric(summary(fit.constant)$coefficients[, 1])
X <- model.matrix(~ 0 + period + trt_indi, data = data2)
tau2 <- summary(fit.constant)$varcor$cluster_id[1]
sigma2 <- summary(fit.constant)$sigma^2
unadj_e <- SV_SWD_e(y=data2$Y_continuous, X, id=data2$cluster_id, beta, tau2, sigma2, n, n_period, trt_model = "constant")
est_unadj_e <- unadj_e[1]
sd_unadj_e <- unadj_e[2]

# unadjusted, nested exchangeable
fit.constant.ne <- lmer(Y_continuous ~ 0 + period + trt_indi + (1 | cluster_id) +  (1 | cluster_period), data = data2)
beta <- as.numeric(summary(fit.constant.ne)$coefficients[, 1])
X <- model.matrix(~ 0 + period + trt_indi, data = data2)
tau2 <- summary(fit.constant.ne)$varcor$cluster_id[1]
sigma2 <- summary(fit.constant.ne)$sigma^2
kappa2 <- summary(fit.constant.ne)$varcor$cluster_period[1]
unadj_ne <- SV_SWD_ne(y=data2$Y_continuous, X, id=data2$cluster_id, beta, tau2, sigma2, kappa2, n, n_period, period = data2$period, trt_model = "constant")
est_unadj_ne <- unadj_ne[1]
sd_unadj_ne <- unadj_ne[2]

# adjusted, independence correlation
fit.constant <- lm(Y_continuous ~ 0 + period + trt_indi + BP_under_140_base + gender + Age + SBP_cal_base + PH_CVD_base + PCT_GOVT_80, data = data2)
beta <- as.numeric(summary(fit.constant)$coefficients[, 1])
X <- model.matrix(~ 0 + period + trt_indi + BP_under_140_base + gender + Age + SBP_cal_base + PH_CVD_base + PCT_GOVT_80, data = data2)
unadj_i <- SV_SWD_i(y=data2$Y_continuous, X, id=data2$cluster_id, beta, n, n_period, trt_model = "constant")
est_adj_i <- unadj_i[1]
sd_adj_i <- unadj_i[2]

# adjusted, exchangeable correlation
fit.constant <- lmer(Y_continuous ~ 0 + period + trt_indi + BP_under_140_base + gender + Age + SBP_cal_base + PH_CVD_base + PCT_GOVT_80 + (1 | cluster_id), data = data2)
beta <- as.numeric(summary(fit.constant)$coefficients[, 1])
X <- model.matrix(~ 0 + period + trt_indi + BP_under_140_base + gender + Age + SBP_cal_base + PH_CVD_base + PCT_GOVT_80, data = data2)
tau2 <- summary(fit.constant)$varcor$cluster_id[1]
sigma2 <- summary(fit.constant)$sigma^2
adj_e <- SV_SWD_e(y=data2$Y_continuous, X, id=data2$cluster_id, beta, tau2, sigma2, n, n_period, trt_model = "constant")
est_adj_e <- adj_e[1]
sd_adj_e <- adj_e[2]

# adjusted, nested exchangeable
fit.constant.ne <- lmer(Y_continuous ~ 0 + period + trt_indi + BP_under_140_base + gender + Age + SBP_cal_base + (1 | cluster_id) +  (1 | cluster_period), data = data2)
beta <- as.numeric(summary(fit.constant.ne)$coefficients[, 1])
X <- model.matrix(~ 0 + period + trt_indi  + BP_under_140_base + gender + Age + SBP_cal_base, data = data2)
tau2 <- summary(fit.constant.ne)$varcor$cluster_id[1]
sigma2 <- summary(fit.constant.ne)$sigma^2
kappa2 <- summary(fit.constant.ne)$varcor$cluster_period[1]
adj_ne <- SV_SWD_ne(y=data2$Y_continuous, X, id=data2$cluster_id, beta, tau2, sigma2, kappa2, 
                    n, n_period, period = data2$period, trt_model = "constant")
est_adj_ne <- adj_ne[1]
sd_adj_ne <- adj_ne[2]

results_constant_te <- matrix(data = c(est_unadj_i, est_unadj_e, est_unadj_ne, est_adj_i, est_adj_e, est_adj_ne,
                                       sd_unadj_i, sd_unadj_e, sd_unadj_ne, sd_adj_i, sd_adj_e, sd_adj_ne), nrow = 6, ncol = 2,
                              dimnames = list(c("unadj_i", "unadj_e", "unadj_ne", "adj_i","adj_e", "adj_ne"), c("Est", "SD")))
round(results_constant_te,2)

## Estimate duration-specific treatment effect------------
## unadjusted, independence correlation
fit.ds <- lm(Y_continuous ~ 0 + period + trt_duration, data = data2)
beta <- as.numeric(summary(fit.ds)$coefficients[, 1])
X <- model.matrix(~ 0 + period + trt_duration, data = data2)
unadj_i <- SV_SWD_i(y=data2$Y_continuous, X, id=data2$cluster_id, beta, n, n_period, trt_model = "duration-specific")

## unadjusted, exchangeable correlation
fit.ds <- lmer(Y_continuous ~ 0 + period + trt_duration + (1 | cluster_id), data = data2)
beta <- as.numeric(summary(fit.ds)$coefficients[, 1])
X <- model.matrix(~ 0 + period + trt_duration, data = data2)
tau2 <- summary(fit.ds)$varcor$cluster_id[1]
sigma2 <- summary(fit.ds)$sigma^2
unadj_e <- SV_SWD_e(y=data2$Y_continuous, X, id=data2$cluster_id, beta, tau2, sigma2, n, n_period, trt_model = "duration-specific")

## unadjusted, nested exchangeable correlation
fit.ds.ne <- lmer(Y_continuous ~ 0 + period + trt_duration + (1 | cluster_id) +  (1 | cluster_period), data = data2)
beta <- as.numeric(summary(fit.ds.ne)$coefficients[, 1])
X <- model.matrix(~ 0 + period + trt_duration, data = data2)
tau2 <- summary(fit.ds.ne)$varcor$cluster_id[1]
sigma2 <- summary(fit.ds.ne)$sigma^2
kappa2 <- summary(fit.ds.ne)$varcor$cluster_period[1]
unadj_ne <- SV_SWD_ne(y=data2$Y_continuous, X, id=data2$cluster_id, beta, tau2, sigma2, kappa2, 
                      n, n_period, period = data2$period, trt_model = "duration-specific")

## adjusted, independence correlation
fit.ds <- lm(Y_continuous ~ 0 + period + trt_duration + BP_under_140_base + gender + Age + SBP_cal_base, data = data2)
beta <- as.numeric(summary(fit.ds)$coefficients[, 1])
X <- model.matrix(~ 0 + period + trt_duration + BP_under_140_base + gender + Age + SBP_cal_base, data = data2)
adj_i <- SV_SWD_i(y=data2$Y_continuous, X, id=data2$cluster_id, beta, n, n_period, trt_model = "duration-specific")

## adjusted, exchangeable correlation
fit.ds <- lmer(Y_continuous ~ 0 + period + trt_duration + BP_under_140_base + gender + Age + SBP_cal_base + (1 | cluster_id), data = data2)
beta <- as.numeric(summary(fit.ds)$coefficients[, 1])
X <- model.matrix(~ 0 + period + trt_duration + BP_under_140_base + gender + Age + SBP_cal_base, data = data2)
tau2 <- summary(fit.ds)$varcor$cluster_id[1]
sigma2 <- summary(fit.ds)$sigma^2
adj_e <- SV_SWD_e(y=data2$Y_continuous, X, id=data2$cluster_id, beta, tau2, sigma2, n, n_period, trt_model = "duration-specific")

## adjusted, nested exchangeable correlation
fit.ds.ne <- lmer(Y_continuous ~ 0 + period + trt_duration + BP_under_140_base + gender + Age + SBP_cal_base + (1 | cluster_id) +  (1 | cluster_period), data = data2)
beta <- as.numeric(summary(fit.ds.ne)$coefficients[, 1])
X <- model.matrix(~ 0 + period + trt_duration + BP_under_140_base + gender + Age + SBP_cal_base, data = data2)
tau2 <- summary(fit.ds.ne)$varcor$cluster_id[1]
sigma2 <- summary(fit.ds.ne)$sigma^2
kappa2 <- summary(fit.ds.ne)$varcor$cluster_period[1]
adj_ne <- SV_SWD_ne(y=data2$Y_continuous, X, id=data2$cluster_id, beta, tau2, sigma2, kappa2, 
                      n, n_period, period = data2$period, trt_model = "duration-specific")


index_set <- c(1,4,2,5,3,6)

results_ds_te <- rbind(unadj_i[index_set], unadj_e[index_set], unadj_ne[index_set], adj_i[index_set], adj_e[index_set], adj_ne[index_set])
rownames(results_ds_te) <- c("unadj_i", "unadj_e", "unadj_ne", "adj_i", "adj_e", "adj_ne")
colnames(results_ds_te) <- c("Est1", "SD1", "Est2", "SD2","Est-avg", "SD-avg")
round(results_ds_te,2)


## Estimate period-specific treatment effect --------------------
# unadjusted, independence correlation
fit.ps <- lm(Y_continuous ~ 0 + period + trt_period, data = data2)
beta <- as.numeric(summary(fit.ps)$coefficients[, 1])
X <- model.matrix(~ 0 + period + trt_period, data = data2)
unadj_i <- SV_SWD_i(y=data2$Y_continuous, X, id=data2$cluster_id, beta, n, n_period, trt_model = "period-specific")

# unadjusted, exchangeable correlation
fit.ps <- lmer(Y_continuous ~ 0 + period + trt_period + (1 | cluster_id), data = data2)
beta <- as.numeric(summary(fit.ps)$coefficients[, 1])
X <- model.matrix(~ 0 + period + trt_period, data = data2)
tau2 <- summary(fit.ps)$varcor$cluster_id[1]
sigma2 <- summary(fit.ps)$sigma^2
unadj_e <- SV_SWD_e(y=data2$Y_continuous, X, id=data2$cluster_id, beta, tau2, sigma2, n, n_period, trt_model = "period-specific")

# unadjusted, nested exchangeable
fit.ps.ne <- lmer(Y_continuous ~ 0 + period + trt_period + (1 | cluster_id) +  (1 | cluster_period), data = data2)
beta <- as.numeric(summary(fit.ps.ne)$coefficients[, 1])
X <- model.matrix(~ 0 + period + trt_period, data = data2)
tau2 <- summary(fit.ps.ne)$varcor$cluster_id[1]
sigma2 <- summary(fit.ps.ne)$sigma^2
kappa2 <- summary(fit.ps.ne)$varcor$cluster_period[1]
unadj_ne <- SV_SWD_ne(y=data2$Y_continuous, X, id=data2$cluster_id, beta, tau2, sigma2, kappa2, 
                      n, n_period, period = data2$period, trt_model = "period-specific")

# adjusted, independence correlation
fit.ps <- lm(Y_continuous ~ 0 + period + trt_period + BP_under_140_base + gender + Age + SBP_cal_base, data = data2)
beta <- as.numeric(summary(fit.ps)$coefficients[, 1])
X <- model.matrix(~ 0 + period + trt_period + BP_under_140_base + gender + Age + SBP_cal_base, data = data2)
adj_i <- SV_SWD_i(y=data2$Y_continuous, X, id=data2$cluster_id, beta, n, n_period, trt_model = "period-specific")

# adjusted, exchangeable correlation
fit.ps <- lmer(Y_continuous ~ 0 + period + trt_period + BP_under_140_base + gender + Age + SBP_cal_base + (1 | cluster_id), data = data2)
beta <- as.numeric(summary(fit.ps)$coefficients[, 1])
X <- model.matrix(~ 0 + period + trt_period + BP_under_140_base + gender + Age + SBP_cal_base, data = data2)
tau2 <- summary(fit.ps)$varcor$cluster_id[1]
sigma2 <- summary(fit.ps)$sigma^2
adj_e <- SV_SWD_e(y=data2$Y_continuous, X, id=data2$cluster_id, beta, tau2, sigma2, n, n_period, trt_model = "period-specific")

# adjusted, nested exchangeable
fit.ps.ne <- lmer(Y_continuous ~ 0 + period + trt_period + BP_under_140_base + gender + Age + SBP_cal_base + (1 | cluster_id) +  (1 | cluster_period), data = data2)
beta <- as.numeric(summary(fit.ps.ne)$coefficients[, 1])
X <- model.matrix(~ 0 + period + trt_period + BP_under_140_base + gender + Age + SBP_cal_base, data = data2)
tau2 <- summary(fit.ps.ne)$varcor$cluster_id[1]
sigma2 <- summary(fit.ps.ne)$sigma^2
kappa2 <- summary(fit.ps.ne)$varcor$cluster_period[1]
adj_ne <- SV_SWD_ne(y=data2$Y_continuous, X, id=data2$cluster_id, beta, tau2, sigma2, kappa2, 
                    n, n_period, period = data2$period, trt_model = "period-specific")

index_set <- c(1, 4, 2, 5, 3, 6)

results_ps_te <- rbind(unadj_i[index_set], unadj_e[index_set], unadj_ne[index_set], adj_i[index_set], adj_e[index_set], adj_ne[index_set])
rownames(results_ps_te) <- c("unadj_i", "unadj_e", "unadj_ne", "adj_i", "adj_e", "adj_ne")
colnames(results_ps_te) <- c("Est1", "SD1", "Est2", "SD2", "Est-avg", "SD-avg")
round(results_ps_te,2)


## Estimate saturated treatment effect --------------------
# unadjusted, independence correlation
fit.saturated <- lm(Y_continuous ~ 0 + period + trt_saturated, data = data2)
beta <- as.numeric(summary(fit.saturated)$coefficients[, 1])
X <- model.matrix(~ 0 + period + trt_saturated, data = data2)
unadj_i <- SV_SWD_i(y=data2$Y_continuous, X, id=data2$cluster_id, beta, n, n_period, trt_model = "saturated")

# unadjusted, exchangeable correlation
fit.saturated <- lmer(Y_continuous ~ 0 + period + trt_saturated + (1 | cluster_id), data = data2)
beta <- as.numeric(summary(fit.saturated)$coefficients[, 1])
X <- model.matrix(~ 0 + period + trt_saturated, data = data2)
tau2 <- summary(fit.saturated)$varcor$cluster_id[1]
sigma2 <- summary(fit.saturated)$sigma^2
unadj_e <- SV_SWD_e(y=data2$Y_continuous, X, id=data2$cluster_id, beta, tau2, sigma2, n, n_period, trt_model = "saturated")

# unadjusted, nested exchangeable
fit.saturated.ne <- lmer(Y_continuous ~ 0 + period + trt_saturated + (1 | cluster_id) +  (1 | cluster_period), data = data2)
beta <- as.numeric(summary(fit.saturated.ne)$coefficients[, 1])
X <- model.matrix(~ 0 + period + trt_saturated, data = data2)
tau2 <- summary(fit.saturated.ne)$varcor$cluster_id[1]
sigma2 <- summary(fit.saturated.ne)$sigma^2
kappa2 <- summary(fit.saturated.ne)$varcor$cluster_period[1]
unadj_ne <- SV_SWD_ne(y=data2$Y_continuous, X, id=data2$cluster_id, beta, tau2, sigma2, kappa2, 
                      n, n_period, period = data2$period, trt_model = "saturated")

# adjusted, independence correlation
fit.saturated <- lm(Y_continuous ~ 0 + period + trt_saturated + BP_under_140_base + gender + Age + SBP_cal_base, data = data2)
beta <- as.numeric(summary(fit.saturated)$coefficients[, 1])
X <- model.matrix(~ 0 + period + trt_saturated + BP_under_140_base + gender + Age + SBP_cal_base, data = data2)
adj_i <- SV_SWD_i(y=data2$Y_continuous, X, id=data2$cluster_id, beta, n, n_period, trt_model = "saturated")

# adjusted, exchangeable correlation
fit.saturated <- lmer(Y_continuous ~ 0 + period + trt_saturated + BP_under_140_base + gender + Age + SBP_cal_base + (1 | cluster_id), data = data2)
beta <- as.numeric(summary(fit.saturated)$coefficients[, 1])
X <- model.matrix(~ 0 + period + trt_saturated + BP_under_140_base + gender + Age + SBP_cal_base, data = data2)
tau2 <- summary(fit.saturated)$varcor$cluster_id[1]
sigma2 <- summary(fit.saturated)$sigma^2
adj_e <- SV_SWD_e(y=data2$Y_continuous, X, id=data2$cluster_id, beta, tau2, sigma2, n, n_period, trt_model = "saturated")

# adjusted, nested exchangeable
fit.saturated.ne <- lmer(Y_continuous ~ 0 + period + trt_saturated + BP_under_140_base + gender + Age + SBP_cal_base + (1 | cluster_id) +  (1 | cluster_period), data = data2)
beta <- as.numeric(summary(fit.saturated.ne)$coefficients[, 1])
X <- model.matrix(~ 0 + period + trt_saturated + BP_under_140_base + gender + Age + SBP_cal_base, data = data2)
tau2 <- summary(fit.saturated.ne)$varcor$cluster_id[1]
sigma2 <- summary(fit.saturated.ne)$sigma^2
kappa2 <- summary(fit.saturated.ne)$varcor$cluster_period[1]
adj_ne <- SV_SWD_ne(y=data2$Y_continuous, X, id=data2$cluster_id, beta, tau2, sigma2, kappa2, 
                    n, n_period, period = data2$period, trt_model = "saturated")

index_set <- c(1,5,2,6,3,7,4,8)

results_saturated_te <- rbind(unadj_i[index_set], unadj_e[index_set], unadj_ne[index_set], adj_i[index_set], adj_e[index_set], adj_ne[index_set])
rownames(results_saturated_te) <- c("unadj_i", "unadj_e", "unadj_ne", "adj_i", "adj_e", "adj_ne")
colnames(results_saturated_te) <- c("Est11", "SD11", "Est21", "SD21", "Est22", "SD22", "Est-avg", "SD-avg")
round(results_saturated_te,2)


########## summary results ---------
summary_idx <- 4:6
results <- rbind(results_constant_te[summary_idx,],
                 results_ds_te[summary_idx,c(1,2)],
                 results_ds_te[summary_idx,c(3,4)],
                 results_ds_te[summary_idx,c(5,6)],
                 results_ps_te[summary_idx,c(1,2)],
                 results_ps_te[summary_idx,c(3,4)],
                 results_ps_te[summary_idx,c(5,6)],
                 results_saturated_te[summary_idx,c(1,2)],
                 results_saturated_te[summary_idx,c(3,4)],
                 results_saturated_te[summary_idx,c(5,6)],
                 results_saturated_te[summary_idx,c(7,8)])
results <- mutate(as.data.frame(results), 
                  ci = paste0("(", round(results[,1] -  qt(0.975, 18) *  results[,2],2), ", ", round(results[,2] + qt(0.975, 18) *  results[,2],2), ")"))
xtable::xtable(results)

results <- cbind(results, 
                 ci.lower = results[,1] -  qt(0.975, 18) *  results[,2],
                 ci.upper = results[,2] + qt(0.975, 18) *  results[,2])

results_constant_te <- cbind(results_constant_te, 
                             ci.lower = results_constant_te[,1] - qt(0.975, 18) * results_constant_te[,2],
                             ci.upper = results_constant_te[,1] + qt(0.975, 18) * results_constant_te[,2])
xtable::xtable(results_constant_te)

results_ds_te_d1 <- cbind(results_ds_te[,c(1,2)],
                       ci.lower = results_ds_te[,1] - qt(0.975, 18) * results_ds_te[,2],
                       ci.upper = results_ds_te[,1] + qt(0.975, 18) * results_ds_te[,2])
results_ds_te_d2 <- cbind(results_ds_te[,c(3,4)],
                          ci.lower = results_ds_te[,3] - qt(0.975, 18) * results_ds_te[,4],
                          ci.upper = results_ds_te[,3] + qt(0.975, 18) * results_ds_te[,4])
results_ds_te_davg <- cbind(results_ds_te[,c(5,6)],
                          ci.lower = results_ds_te[,5] - qt(0.975, 18) * results_ds_te[,6],
                          ci.upper = results_ds_te[,5] + qt(0.975, 18) * results_ds_te[,6])
xtable::xtable(results_ds_te_d1)
xtable::xtable(results_ds_te_d2)
xtable::xtable(results_ds_te_davg)


results_ps_te_p1 <- cbind(results_ps_te[,c(1,2)],
                          ci.lower = results_ps_te[,1] - qt(0.975, 18) * results_ps_te[,2],
                          ci.upper = results_ps_te[,1] + qt(0.975, 18) * results_ps_te[,2])
results_ps_te_p2 <- cbind(results_ps_te[,c(3,4)],
                          ci.lower = results_ps_te[,3] - qt(0.975, 18) * results_ps_te[,4],
                          ci.upper = results_ps_te[,3] + qt(0.975, 18) * results_ps_te[,4])
results_ps_te_pavg <- cbind(results_ps_te[,c(5,6)],
                          ci.lower = results_ps_te[,5] - qt(0.975, 18) * results_ps_te[,6],
                          ci.upper = results_ps_te[,5] + qt(0.975, 18) * results_ps_te[,6])
xtable::xtable(results_ps_te_p1)
xtable::xtable(results_ps_te_p2)
xtable::xtable(results_ps_te_pavg)

results_saturated_te_s11 <- cbind(results_saturated_te[,c(1,2)],
                          ci.lower = results_saturated_te[,1] - qt(0.975, 18) * results_saturated_te[,2],
                          ci.upper = results_saturated_te[,1] + qt(0.975, 18) * results_saturated_te[,2])
results_saturated_te_s21 <- cbind(results_saturated_te[,c(3,4)],
                          ci.lower = results_saturated_te[,3] - qt(0.975, 18) * results_saturated_te[,4],
                          ci.upper = results_saturated_te[,3] + qt(0.975, 18) * results_saturated_te[,4])
results_saturated_te_s22 <- cbind(results_saturated_te[,c(5,6)],
                          ci.lower = results_saturated_te[,5] - qt(0.975, 18) * results_saturated_te[,6],
                          ci.upper = results_saturated_te[,5] + qt(0.975, 18) * results_saturated_te[,6])
results_saturated_te_savg <- cbind(results_saturated_te[,c(7,8)],
                            ci.lower = results_saturated_te[,7] - qt(0.975, 18) * results_saturated_te[,8],
                            ci.upper = results_saturated_te[,7] + qt(0.975, 18) * results_saturated_te[,8])
xtable::xtable(results_saturated_te_s11)
xtable::xtable(results_saturated_te_s21)
xtable::xtable(results_saturated_te_s22)
xtable::xtable(results_saturated_te_savg)

# likelihood ratio test -----------
lrt <- function(model1.fit, model2.fit, df_diff){
  teststat <- as.numeric(-2 * (logLik(model1.fit) - logLik(model2.fit)))
  pchisq(teststat, df_diff, lower.tail = F)
}

fit.constant <- lm(Y_continuous ~ 0 + period + trt_indi, data = data2)
fit.ds <- lm(Y_continuous ~ 0 + period + trt_duration, data = data2)
fit.ps <- lm(Y_continuous ~ 0 + period + trt_period, data = data2)
fit.saturated <- lm(Y_continuous ~ 0 + period + trt_saturated, data = data2)
lrt(fit.constant, fit.ds, 1)
lrt(fit.constant, fit.ps, 1)
lrt(fit.ds, fit.saturated, 1)
lrt(fit.ps, fit.saturated, 1)

fit.constant <- lmer(Y_continuous ~ 0 + period + trt_indi + (1 | cluster_id), data = data2)
fit.ds <- lmer(Y_continuous ~ 0 + period + trt_duration + (1 | cluster_id), data = data2)
fit.ps <- lmer(Y_continuous ~ 0 + period + trt_period + (1 | cluster_id), data = data2)
fit.saturated <- lmer(Y_continuous ~ 0 + period + trt_saturated + (1 | cluster_id), data = data2)
lrt(fit.constant, fit.ds, 1)
lrt(fit.constant, fit.ps, 1)
lrt(fit.ds, fit.saturated, 1)
lrt(fit.ps, fit.saturated, 1)


fit.constant.ne <- lmer(Y_continuous ~ 0 + period + trt_indi + (1 | cluster_id) +  (1 | cluster_period), data = data2)
fit.ds.ne <- lmer(Y_continuous ~ 0 + period + trt_duration + (1 | cluster_id) +  (1 | cluster_period), data = data2)
fit.ps.ne <- lmer(Y_continuous ~ 0 + period + trt_period + (1 | cluster_id) +  (1 | cluster_period), data = data2)
fit.saturated.ne <- lmer(Y_continuous ~ 0 + period + trt_saturated + (1 | cluster_id) +  (1 | cluster_period), data = data2)
lrt(fit.constant.ne, fit.ds.ne, 1)
lrt(fit.constant.ne, fit.ps.ne, 1)
lrt(fit.ds.ne, fit.saturated.ne, 1)
lrt(fit.ps.ne, fit.saturated.ne, 1)
lrt(fit.constant.ne, fit.saturated.ne, 2)

