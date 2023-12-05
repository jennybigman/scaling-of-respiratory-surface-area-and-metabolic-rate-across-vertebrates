# Bigman et al. metabolic rate and respiratory surface area scaling Science Advances

# 04 - 'Combined models'

library(rstan)
library(beepr)
library(loo)

# set options
rstan_options(auto_write = TRUE)

options(mc.cores = parallel::detectCores())

d_mat <- diag(1, 109, 109)

## model names correspond to Table S1 ##

## C1 ##

dat <- list(N=nrow(Overlap.MR.RSA),
            RSA=Overlap.MR.RSA[,'LogRSAcm2'],
            MRSA=Overlap.MR.RSA[,'LogCenteredMeanMassRSA'],
            MR=Overlap.MR.RSA[, 'LogMeanWholeOrganismMRWatts'],
            MMR=Overlap.MR.RSA[, 'LogCenteredMeanMassMR'],
            Temp = Overlap.MR.RSA[, 'Scaled.Centered.MR.InverseTemp'],
            d_mat = d_mat,
            vcov_mat = vcov_mat)

C_MOD1 <- stan(file = "Combined_Mod1.stan",
               data = dat,
               iter=5000,
               warmup = 1000,
               thin = 10,
               chains=4,
               control=list("adapt_delta"=0.81),
               pars=c("aMR",
                      "aRSA",
                      "bMassMR",
                      "bMassRSA",
                      "bResid",
                      "bTemp",
                      "lambda_MR",
                      "lambda_RSA",
                      "sigma_MR",
                      "sigma_RSA",
                      "resid",
                      "log_lik",
                      "sigma_total_MR",
                      "sigma_total_RSA"))


## C2 ##

dat <- list(N=nrow(Overlap.MR.RSA),
            RSA=Overlap.MR.RSA[,'LogRSAcm2'],
            MRSA=Overlap.MR.RSA[,'LogCenteredMeanMassRSA'],
            MR=Overlap.MR.RSA[, 'LogMeanWholeOrganismMRWatts'],
            MMR=Overlap.MR.RSA[, 'LogCenteredMeanMassMR'],
            Temp = Overlap.MR.RSA[, 'Scaled.Centered.MR.InverseTemp'],
            d_mat = d_mat,
            vcov_mat = vcov_mat)

C_MOD2 <- stan(file = "Combined_Mod2.stan",
               data = dat,
               iter= 5000,
               warmup = 1000,
               thin = 10,
               chains=4,
               control=list("adapt_delta"=0.81),
               pars=c("aMR",
                     "aRSA",
                     "bMassMR",
                     "bMassRSA",
                     "bResid",
                     "bTemp",
                     "bIntMassResid",
                     "lambda_MR",
                     "lambda_RSA",
                     "sigma_MR",
                     "sigma_RSA",
                     "resid",
                     "log_lik",
                     "sigma_total_MR",
                     "sigma_total_RSA"))


## C3 ##

dat <- list(N=nrow(Overlap.MR.RSA),
            RSA=Overlap.MR.RSA[,'LogRSAcm2'],
            MRSA=Overlap.MR.RSA[,'LogCenteredMeanMassRSA'],
            MR=Overlap.MR.RSA[, 'LogMeanWholeOrganismMRWatts'],
            MMR=Overlap.MR.RSA[, 'LogCenteredMeanMassMR'],
            Temp = Overlap.MR.RSA[, 'Scaled.Centered.MR.InverseTemp'],
            ThermoStrat = Overlap.MR.RSA[, 'ThermoStrat'],
            d_mat = d_mat,
            vcov_mat = vcov_mat)

C_MOD3 <- stan(file = "Combined_Mod3.stan",
               data = dat,
               iter= 5000,
               warmup = 1000,
               thin = 10,
               chains=4,
               control=list("adapt_delta"=0.81),
               pars=c("aMR",
                      "aRSA",
                      "bMassMR",
                      "bMassRSA",
                      "bResid",
                      "bTemp",
                      "bTherm",
                      "lambda_MR",
                      "lambda_RSA",
                      "sigma_MR",
                      "sigma_RSA",
                      "resid",
                      "log_lik",
                      "sigma_total_MR",
                      "sigma_total_RSA"))


## C4 ##

dat <- list(N=nrow(Overlap.MR.RSA),
            RSA=Overlap.MR.RSA[,'LogRSAcm2'],
            MRSA=Overlap.MR.RSA[,'LogCenteredMeanMassRSA'],
            MR=Overlap.MR.RSA[, 'LogMeanWholeOrganismMRWatts'],
            MMR=Overlap.MR.RSA[, 'LogCenteredMeanMassMR'],
            Temp = Overlap.MR.RSA[, 'Scaled.Centered.MR.InverseTemp'],
            ThermoStrat = Overlap.MR.RSA[, "ThermoStrat"],
            d_mat = d_mat,
            vcov_mat = vcov_mat)

C_MOD4 <- stan(file = "Combined_Mod4.stan",
               data = dat,
               iter= 5000,
               warmup = 1000,
               thin = 10,
               chains=4,
               control=list("adapt_delta"=0.81),
               pars=c("aMR",
                      "aRSA",
                      "bMassMR",
                      "bMassRSA",
                      "bResid",
                      "bTemp",
                      "bTherm",
                      "bIntMassResid",
                      "lambda_MR",
                      "lambda_RSA",
                      "resid",
                      "sigma_MR",
                      "sigma_RSA",
                      "log_lik",
                      "sigma_total_MR",
                      "sigma_total_RSA"))


## C5 ##

dat <- list(N=nrow(Overlap.MR.RSA),
            RSA=Overlap.MR.RSA[,'LogRSAcm2'],
            MRSA=Overlap.MR.RSA[,'LogCenteredMeanMassRSA'],
            MR=Overlap.MR.RSA[, 'LogMeanWholeOrganismMRWatts'],
            MMR=Overlap.MR.RSA[, 'LogCenteredMeanMassMR'],
            Temp = Overlap.MR.RSA[, 'Scaled.Centered.MR.InverseTemp'],
            ThermoStrat = Overlap.MR.RSA[, "ThermoStrat"],
            d_mat = d_mat,
            vcov_mat = vcov_mat)

C_MOD5 <- stan(file = "Combined_Mod5.stan",
               data = dat,
               iter= 5000,
               warmup = 1000,
               chains=4,
               thin = 10,
               control=list("adapt_delta"=0.81),
               pars=c("aMR",
                       "aRSA",
                       "bMassMR",
                       "bMassRSA",
                       "bResid",
                       "resid_RSA",
                       "resid_MR",
                       "mu_MR",
                       "bTemp",
                       "bTherm",
                       "bIntMassTS",
                       "lambda_MR",
                       "lambda_RSA",
                       "sigma_MR",
                       "sigma_RSA",
                       "log_lik",
                       "sigma_total_MR",
                       "sigma_total_RSA"))


## C6 ##

dat <- list(N=nrow(Overlap.MR.RSA),
            RSA=Overlap.MR.RSA[,'LogRSAcm2'],
            MRSA=Overlap.MR.RSA[,'LogCenteredMeanMassRSA'],
            MR=Overlap.MR.RSA[, 'LogMeanWholeOrganismMRWatts'],
            MMR=Overlap.MR.RSA[, 'LogCenteredMeanMassMR'],
            Temp = Overlap.MR.RSA[, 'Scaled.Centered.MR.InverseTemp'],
            ThermoStrat = Overlap.MR.RSA$ThermoStrat,
            d_mat = d_mat,
            vcov_mat = vcov_mat)

C_MOD6 <- stan(file = "Combined_Mod6.stan",
               data = dat,
               iter= 5000,
               warmup = 1000,
               thin = 10,
               chains=4,
               control=list("adapt_delta"=0.81),
               pars=c("aMR",
                      "aRSA",
                      "bMassMR",
                      "bMassRSA",
                      "bResid",
                      "resid",
                      "bTemp",
                      "bTherm",
                      "bIntMassTS",
                      "bIntMassResid",
                      "lambda_MR",
                      "lambda_RSA",
                      "sigma_MR",
                      "sigma_RSA",
                      "log_lik",
                      "sigma_total_MR",
                      "sigma_total_RSA"))


## model comparison ####

## C1 ##

posterior_MR_RSA1 <- as.data.frame(C_MOD1)

residual_rsa_dist <- as.matrix(posterior_MR_RSA1[grep("resid", colnames(posterior_MR_RSA1))])
residualRSA <- colMeans(residual_rsa_dist)
Overlap.MR.RSA$residualRSA <- residualRSA

y_MR_RSA1 <- Overlap.MR.RSA$LogMeanWholeOrganismMRWatts
N_MR_RSA1 <- length(y_MR_RSA1)
S_MR_RSA1 <- nrow(posterior_MR_RSA1)
loglik_MR_RSA1 <- yloo_MR_RSA1 <- sdloo_MR_RSA1 <- matrix(nrow = S_MR_RSA1,
                                                          ncol = N_MR_RSA1)

st_MR_RSA1 <- as.matrix(posterior_MR_RSA1[grep("sigma_total_MR",
                                               colnames(posterior_MR_RSA1))])
str(st_MR_RSA1)
dim(st_MR_RSA1)

sigma_array_MR_RSA1 <- array(st_MR_RSA1, dim = c(nrow(st_MR_RSA1), 109, 109))
str(sigma_array_MR_RSA1)
dim(sigma_array_MR_RSA1)

s_MR_RSA1 <- sigma_array_MR_RSA1[2, ,]
s2_MR_RSA1 <- st_MR_RSA1[2, ]
s_MR_RSA1 == s2_MR_RSA1

for (s in 1:S_MR_RSA1) {
  p_MR_RSA1 <- posterior_MR_RSA1[s, ]
  eta_MR_RSA1 <- p_MR_RSA1$aMR + p_MR_RSA1$bMassMR * Overlap.MR.RSA$LogCenteredMeanMassMR +
                 p_MR_RSA1$bResid * Overlap.MR.RSA$residualRSA +
                 p_MR_RSA1$bTemp * Overlap.MR.RSA$Scaled.Centered.MR.InverseTemp
  Cinv_MR_RSA1 <- solve(sigma_array_MR_RSA1[s, , ])
  g_MR_RSA1 <- Cinv_MR_RSA1 %*% (y_MR_RSA1 - eta_MR_RSA1)
  cbar_MR_RSA1 <- diag(Cinv_MR_RSA1)
  yloo_MR_RSA1[s, ] <- y_MR_RSA1 - g_MR_RSA1 / cbar_MR_RSA1
  sdloo_MR_RSA1[s, ] <- sqrt(1 / cbar_MR_RSA1)
  loglik_MR_RSA1[s, ] <- dnorm(y_MR_RSA1, yloo_MR_RSA1[s, ], sdloo_MR_RSA1[s, ], log = TRUE)
}

log_ratios_MR_RSA1 <- -loglik_MR_RSA1
r_eff_MR_RSA1 <- relative_eff(exp(loglik_MR_RSA1), chain_id = rep(1:4, each = 800))
psis_result_MR_RSA1 <- psis(log_ratios_MR_RSA1, r_eff = r_eff_MR_RSA1)

plot(psis_result_MR_RSA1, label_points = TRUE)

(psis_loo_MR_RSA1 <- loo(loglik_MR_RSA1))


## C2 ##

posterior_MR_RSA2 <- as.data.frame(C_MOD2)

residual_rsa_dist <- as.matrix(posterior_MR_RSA2[grep("resid", colnames(posterior_MR_RSA2))])
residualRSA <- colMeans(residual_rsa_dist)
Overlap.MR.RSA$residualRSA <- residualRSA

y_MR_RSA2 <- Overlap.MR.RSA$LogMeanWholeOrganismMRWatts
N_MR_RSA2 <- length(y_MR_RSA2)
S_MR_RSA2 <- nrow(posterior_MR_RSA2)
loglik_MR_RSA2 <- yloo_MR_RSA2 <- sdloo_MR_RSA2 <- matrix(nrow = S_MR_RSA2,
                                                          ncol = N_MR_RSA2)

st_MR_RSA2 <- as.matrix(posterior_MR_RSA2[grep("sigma_total_MR",
                                               colnames(posterior_MR_RSA2))])
str(st_MR_RSA2)
dim(st_MR_RSA2)

sigma_array_MR_RSA2 <- array(st_MR_RSA2, dim = c(nrow(st_MR_RSA2), 109, 109))
str(sigma_array_MR_RSA2)
dim(sigma_array_MR_RSA2)

s_MR_RSA2 <- sigma_array_MR_RSA2[2, ,]
s2_MR_RSA2 <- st_MR_RSA2[2, ]
s_MR_RSA2 == s2_MR_RSA2

for (s in 1:s_MR_RSA2) {
  p_MR_RSA2 <- posterior_MR_RSA2[s, ]
  eta_MR_RSA2 <- p_MR_RSA2$aMR + p_MR_RSA2$bMassMR * Overlap.MR.RSA$LogCenteredMeanMassMR +
                 p_MR_RSA2$bResid * Overlap.MR.RSA$residualRSA +
                 p_MR_RSA2$bIntMassResid * Overlap.MR.RSA$LogCenteredMeanMassMR * Overlap.MR.RSA$residualRSA +
                 p_MR_RSA2$bTemp * Overlap.MR.RSA$Scaled.Centered.MR.InverseTemp
  Cinv_MR_RSA2 <- solve(sigma_array_MR_RSA2[s, , ])
  g_MR_RSA2 <- Cinv_MR_RSA2 %*% (y_MR_RSA2 - eta_MR_RSA2)
  cbar_MR_RSA2 <- diag(Cinv_MR_RSA2)
  yloo_MR_RSA2[s, ] <- y_MR_RSA2 - g_MR_RSA2 / cbar_MR_RSA2
  sdloo_MR_RSA2[s, ] <- sqrt(1 / cbar_MR_RSA2)
  loglik_MR_RSA2[s, ] <- dnorm(y_MR_RSA2, yloo_MR_RSA2[s, ], sdloo_MR_RSA2[s, ], log = TRUE)
}

log_ratios_MR_RSA2 <- -loglik_MR_RSA2
r_eff_MR_RSA2 <- relative_eff(exp(loglik_MR_RSA2), chain_id = rep(1:4, each = 800))
psis_result_MR_RSA2 <- psis(log_ratios_MR_RSA2, r_eff = r_eff_MR_RSA2)

plot(psis_result_MR_RSA2, label_points = TRUE)

(psis_loo_MR_RSA2 <- loo(loglik_MR_RSA2))


## C3 ##

posterior_MR_RSA3 <- as.data.frame(C_MOD3)

residual_rsa_dist <- as.matrix(posterior_MR_RSA3[grep("resid", colnames(posterior_MR_RSA3))])
residualRSA <- colMeans(residual_rsa_dist)
Overlap.MR.RSA$residualRSA <- residualRSA

y_MR_RSA3 <- Overlap.MR.RSA$LogMeanWholeOrganismMRWatts
N_MR_RSA3 <- length(y_MR_RSA3)
S_MR_RSA3 <- nrow(posterior_MR_RSA3)
loglik_MR_RSA3 <- yloo_MR_RSA3 <- sdloo_MR_RSA3 <- matrix(nrow = S_MR_RSA3,
                                                          ncol = N_MR_RSA3)

st_MR_RSA3 <- as.matrix(posterior_MR_RSA3[grep("sigma_total_MR",
                                               colnames(posterior_MR_RSA3))])
str(st_MR_RSA3)
dim(st_MR_RSA3)

sigma_array_MR_RSA3 <- array(st_MR_RSA3, dim = c(nrow(st_MR_RSA3), 109, 109))
str(sigma_array_MR_RSA3)
dim(sigma_array_MR_RSA3)

s_MR_RSA3 <- sigma_array_MR_RSA3[2, ,]
s2_MR_RSA3 <- st_MR_RSA3[2, ]
s_MR_RSA3 == s2_MR_RSA3

for (s in 1:S_MR_RSA3) {
  p_MR_RSA3 <- posterior_MR_RSA3[s, ]
  eta_MR_RSA3 <- p_MR_RSA3$aMR + p_MR_RSA3$bMassMR * Overlap.MR.RSA$LogCenteredMeanMassMR +
                 p_MR_RSA3$bResid * Overlap.MR.RSA$residualRSA +
                 p_MR_RSA3$bTemp * Overlap.MR.RSA$Scaled.Centered.MR.InverseTemp +
                 p_MR_RSA3$bTherm * Overlap.MR.RSA$ThermoStrat
  Cinv_MR_RSA3 <- solve(sigma_array_MR_RSA3[s, , ])
  g_MR_RSA3 <- Cinv_MR_RSA3 %*% (y_MR_RSA3 - eta_MR_RSA3)
  cbar_MR_RSA3 <- diag(Cinv_MR_RSA3)
  yloo_MR_RSA3[s, ] <- y_MR_RSA3 - g_MR_RSA3 / cbar_MR_RSA3
  sdloo_MR_RSA3[s, ] <- sqrt(1 / cbar_MR_RSA3)
  loglik_MR_RSA3[s, ] <- dnorm(y_MR_RSA3, yloo_MR_RSA3[s, ], sdloo_MR_RSA3[s, ], log = TRUE)
}

log_ratios_MR_RSA3 <- -loglik_MR_RSA3
r_eff_MR_RSA3 <- relative_eff(exp(loglik_MR_RSA3), chain_id = rep(1:4, each = 800))
psis_result_MR_RSA3 <- psis(log_ratios_MR_RSA3, r_eff = r_eff_MR_RSA3)

plot(psis_result_MR_RSA3, label_points = TRUE)

(psis_loo_MR_RSA3 <- loo(loglik_MR_RSA3))


## C4 ##

posterior_MR_RSA4 <- as.data.frame(C_MOD4)

residual_rsa_dist <- as.matrix(posterior_MR_RSA4[grep("resid", colnames(posterior_MR_RSA4))])
residualRSA <- colMeans(residual_rsa_dist)
Overlap.MR.RSA$residualRSA <- residualRSA

y_MR_RSA4 <- Overlap.MR.RSA$LogMeanWholeOrganismMRWatts
N_MR_RSA4 <- length(y_MR_RSA4)
S_MR_RSA4 <- nrow(posterior_MR_RSA4)
loglik_MR_RSA4 <- yloo_MR_RSA4 <- sdloo_MR_RSA4 <- matrix(nrow = S_MR_RSA4,
                                                          ncol = N_MR_RSA4)

st_MR_RSA4 <- as.matrix(posterior_MR_RSA4[grep("sigma_total_MR",
                                               colnames(posterior_MR_RSA4))])
str(st_MR_RSA4)
dim(st_MR_RSA4)

sigma_array_MR_RSA4 <- array(st_MR_RSA4, dim = c(nrow(st_MR_RSA4), 109, 109))
str(sigma_array_MR_RSA4)
dim(sigma_array_MR_RSA4)

s_MR_RSA4 <- sigma_array_MR_RSA4[2, ,]
s2_MR_RSA4 <- st_MR_RSA4[2, ]
s_MR_RSA4 == s2_MR_RSA4

for (s in 1:S_MR_RSA4) {
  p_MR_RSA4 <- posterior_MR_RSA4[s, ]
  eta_MR_RSA4 <- p_MR_RSA4$aMR + p_MR_RSA4$bMassMR * Overlap.MR.RSA$LogCenteredMeanMassMR +
                 p_MR_RSA4$bResid * Overlap.MR.RSA$residualRSA +
                 p_MR_RSA4$bIntMassResid * Overlap.MR.RSA$LogCenteredMeanMassMR * Overlap.MR.RSA$residualRSA +
                 p_MR_RSA4$bTemp * Overlap.MR.RSA$Scaled.Centered.MR.InverseTemp +
                 p_MR_RSA4$bTherm * Overlap.MR.RSA$ThermoStrat
  Cinv_MR_RSA4 <- solve(sigma_array_MR_RSA4[s, , ])
  g_MR_RSA4 <- Cinv_MR_RSA4 %*% (y_MR_RSA4 - eta_MR_RSA4)
  cbar_MR_RSA4 <- diag(Cinv_MR_RSA4)
  yloo_MR_RSA4[s, ] <- y_MR_RSA4 - g_MR_RSA4 / cbar_MR_RSA4
  sdloo_MR_RSA4[s, ] <- sqrt(1 / cbar_MR_RSA4)
  loglik_MR_RSA4[s, ] <- dnorm(y_MR_RSA4, yloo_MR_RSA4[s, ], sdloo_MR_RSA4[s, ], log = TRUE)
}

log_ratios_MR_RSA4 <- -loglik_MR_RSA4
r_eff_MR_RSA4 <- relative_eff(exp(loglik_MR_RSA4), chain_id = rep(1:4, each = 800))
psis_result_MR_RSA4 <- psis(log_ratios_MR_RSA4, r_eff = r_eff_MR_RSA4)

plot(psis_result_MR_RSA4, label_points = TRUE)

(psis_loo_MR_RSA4 <- loo(loglik_MR_RSA4))


## C5 ##

posterior_MR_RSA5 <- as.data.frame(C_MOD5)

residual_rsa_dist <- as.matrix(posterior_MR_RSA5[grep("resid", colnames(posterior_MR_RSA5))])
residualRSA <- colMeans(residual_rsa_dist)
Overlap.MR.RSA$residualRSA <- residualRSA

y_MR_RSA5 <- Overlap.MR.RSA$LogMeanWholeOrganismMRWatts
N_MR_RSA5 <- length(y_MR_RSA5)
S_MR_RSA5 <- nrow(posterior_MR_RSA5)
loglik_MR_RSA5 <- yloo_MR_RSA5 <- sdloo_MR_RSA5 <- matrix(nrow = S_MR_RSA5, ncol = N_MR_RSA5)

st_MR_RSA5 <- as.matrix(posterior_MR_RSA5[grep("sigma_total_MR", colnames(posterior_MR_RSA5))])
str(st_MR_RSA5)
dim(st_MR_RSA5)

sigma_array_MR_RSA5 <- array(st_MR_RSA5, dim = c(nrow(st_MR_RSA5), 109, 109))
str(sigma_array_MR_RSA5)
dim(sigma_array_MR_RSA5)

s1_MR_RSA5 <- sigma_array_MR_RSA5[2, ,]
s2_MR_RSA5 <- st_MR_RSA5[2, ]
s1_MR_RSA5 == s2_MR_RSA5

for (s in 1:S_MR_RSA5) {
  p_MR_RSA5 <- posterior_MR_RSA5[s, ]
  eta_MR_RSA5 <- p_MR_RSA5$aMR +
                 p_MR_RSA5$bResid * Overlap.MR.RSA$residualRSA +
                 p_MR_RSA5$bMassMR * Overlap.MR.RSA$LogCenteredMeanMassMR +
                 p_MR_RSA5$bTherm * Overlap.MR.RSA$ThermoStrat  +
                 p_MR_RSA5$bIntMassTS * Overlap.MR.RSA$LogCenteredMeanMassMR * Overlap.MR.RSA$ThermoStrat +
                 p_MR_RSA5$bTemp * Overlap.MR.RSA$Scaled.Centered.MR.InverseTemp
  Cinv_MR_RSA5 <- solve(sigma_array_MR_RSA5[s, , ])
  g_MR_RSA5 <- Cinv_MR_RSA5 %*% (y_MR_RSA5 - eta_MR_RSA5)
  cbar_MR_RSA5 <- diag(Cinv_MR_RSA5)
  yloo_MR_RSA5[s, ] <- y_MR_RSA5 - g_MR_RSA5 / cbar_MR_RSA5
  sdloo_MR_RSA5[s, ] <- sqrt(1 / cbar_MR_RSA5)
  loglik_MR_RSA5[s, ] <- dnorm(y_MR_RSA5, yloo_MR_RSA5[s, ], sdloo_MR_RSA5[s, ], log = TRUE)
}

log_ratios_MR_RSA5 <- -loglik_MR_RSA5
psis_result_MR_RSA5 <- psis(log_ratios_MR_RSA5)

plot(psis_result_MR_RSA5, label_points = TRUE)

(psis_loo_MR_RSA5 <- loo(loglik_MR_RSA5))


## C6 ##

posterior_MR_RSA6 <- as.data.frame(C_MOD6)

residual_rsa_dist <- as.matrix(posterior_MR_RSA6[grep("resid", colnames(posterior_MR_RSA6))])
residualRSA <- colMeans(residual_rsa_dist)
Overlap.MR.RSA$residualRSA <- residualRSA

y_MR_RSA6 <- Overlap.MR.RSA$LogMeanWholeOrganismMRWatts
N_MR_RSA6 <- length(y_MR_RSA6)
S_MR_RSA6 <- nrow(posterior_MR_RSA6)
loglik_MR_RSA6 <- yloo_MR_RSA6 <- sdloo_MR_RSA6 <- matrix(nrow = S_MR_RSA6,
                                                          ncol = N_MR_RSA6)

st_MR_RSA6 <- as.matrix(posterior_MR_RSA6[grep("sigma_total_MR",
                                               colnames(posterior_MR_RSA6))])
str(st_MR_RSA6)
dim(st_MR_RSA6)

sigma_array_MR_RSA6 <- array(st_MR_RSA6, dim = c(nrow(st_MR_RSA6), 109, 109))
str(sigma_array_MR_RSA6)
dim(sigma_array_MR_RSA6)

s_MR_RSA6 <- sigma_array_MR_RSA6[2, ,]
s2_MR_RSA6 <- st_MR_RSA6[2, ]
s_MR_RSA6 == s2_MR_RSA6

for (s in 1:S_MR_RSA6) {
  p_MR_RSA6 <- posterior_MR_RSA6[s, ]
  eta_MR_RSA6 <- p_MR_RSA6$aMR +
                 p_MR_RSA6$bResid * Overlap.MR.RSA$residualRSA +
                 p_MR_RSA6$bMassMR * Overlap.MR.RSA$LogCenteredMeanMassMR +
                 p_MR_RSA6$bTherm * Overlap.MR.RSA$ThermoStrat +
                 p_MR_RSA6$bTemp * Overlap.MR.RSA$Scaled.Centered.MR.InverseTemp +
                 p_MR_RSA6$bIntMassResid * Overlap.MR.RSA$LogCenteredMeanMassMR * Overlap.MR.RSA$residualRSA +
                 p_MR_RSA6$bIntMassTS * Overlap.MR.RSA$LogCenteredMeanMassMR * Overlap.MR.RSA$ThermoStrat
  Cinv_MR_RSA6 <- solve(sigma_array_MR_RSA6[s, , ])
  g_MR_RSA6 <- Cinv_MR_RSA6 %*% (y_MR_RSA6 - eta_MR_RSA6)
  cbar_MR_RSA6 <- diag(Cinv_MR_RSA6)
  yloo_MR_RSA6[s, ] <- y_MR_RSA6 - g_MR_RSA6 / cbar_MR_RSA6
  sdloo_MR_RSA6[s, ] <- sqrt(1 / cbar_MR_RSA6)
  loglik_MR_RSA6[s, ] <- dnorm(y_MR_RSA6, yloo_MR_RSA6[s, ], sdloo_MR_RSA6[s, ], log = TRUE)
}

log_ratios_MR_RSA6 <- -loglik_MR_RSA6
r_eff_MR_RSA6 <- relative_eff(exp(loglik_MR_RSA6), chain_id = rep(1:4, each = 800))
psis_result_MR_RSA6 <- psis(log_ratios_MR_RSA6, r_eff = r_eff_MR_RSA6)

plot(psis_result_MR_RSA6, label_points = TRUE)

(psis_loo_MR_RSA6 <- loo(loglik_MR_RSA6))


## compare models

Combined_model_comparison <- loo_compare(psis_loo_MR_RSA1, psis_loo_MR_RSA2,
                                         psis_loo_MR_RSA3, psis_loo_MR_RSA4,
                                         psis_loo_MR_RSA5, psis_loo_MR_RSA6)


# weights
loo_list_combined <- list(psis_loo_MR_RSA1, psis_loo_MR_RSA2,
                          psis_loo_MR_RSA3, psis_loo_MR_RSA4,
                          psis_loo_MR_RSA5, psis_loo_MR_RSA6)

combined_mod_wt <- loo_model_weights(loo_list_combined)


#### models with respiratory organ instead of thermoregulatory strategy ####

## C3_LG ##

dat <- list(N=nrow(Overlap.MR.RSA),
            RSA=Overlap.MR.RSA[,'LogRSAcm2'],
            MRSA=Overlap.MR.RSA[,'LogCenteredMeanMassRSA'],
            MR=Overlap.MR.RSA[, 'LogMeanWholeOrganismMRWatts'],
            MMR=Overlap.MR.RSA[, 'LogCenteredMeanMassMR'],
            Temp = Overlap.MR.RSA[, 'Scaled.Centered.MR.InverseTemp'],
            ThermoStrat = Overlap.MR.RSA[, 'LungsGills'],
            d_mat = d_mat,
            vcov_mat = vcov_mat)

C_MOD3_LG <- stan(file = "Combined_Mod3.stan",
               data = dat,
               iter=10000,
               warmup = 2000,
               chains=4,
               control=list("adapt_delta"=0.81),
               pars=c("aMR",
                      "aRSA",
                      "bMassMR",
                      "bMassRSA",
                      "bResid",
                      "bTemp",
                      "bTherm",
                      "lambda_MR",
                      "lambda_RSA",
                      "sigma_MR",
                      "sigma_RSA",
                      "resid",
                      "log_lik",
                      "sigma_total_MR",
                      "sigma_total_RSA"))

## C4_LG ##

dat <- list(N=nrow(Overlap.MR.RSA),
            RSA=Overlap.MR.RSA[,'LogRSAcm2'],
            MRSA=Overlap.MR.RSA[,'LogCenteredMeanMassRSA'],
            MR=Overlap.MR.RSA[, 'LogMeanWholeOrganismMRWatts'],
            MMR=Overlap.MR.RSA[, 'LogCenteredMeanMassMR'],
            Temp = Overlap.MR.RSA[, 'Scaled.Centered.MR.InverseTemp'],
            ThermoStrat = Overlap.MR.RSA[, "LungsGills"],
            d_mat = d_mat,
            vcov_mat = vcov_mat)

C_MOD4_LG <- stan(file = "Combined_Mod4.stan",
               data = dat,
               iter=10000,
               warmup = 2000,
               chains=4,
               control=list("adapt_delta"=0.81),
               pars=c("aMR",
                      "aRSA",
                      "bMassMR",
                      "bMassRSA",
                      "bResid",
                      "bTemp",
                      "bTherm",
                      "bIntMassResid",
                      "lambda_MR",
                      "lambda_RSA",
                      "resid",
                      "sigma_MR",
                      "sigma_RSA",
                      "log_lik",
                      "sigma_total_MR",
                      "sigma_total_RSA"))


## C5_LG ##

dat <- list(N=nrow(Overlap.MR.RSA),
            RSA=Overlap.MR.RSA[,'LogRSAcm2'],
            MRSA=Overlap.MR.RSA[,'LogCenteredMeanMassRSA'],
            MR=Overlap.MR.RSA[, 'LogMeanWholeOrganismMRWatts'],
            MMR=Overlap.MR.RSA[, 'LogCenteredMeanMassMR'],
            Temp = Overlap.MR.RSA[, 'Scaled.Centered.MR.InverseTemp'],
            ThermoStrat = Overlap.MR.RSA[, "LungsGills"],
            d_mat = d_mat,
            vcov_mat = vcov_mat)

C_MOD5_LG <- stan(file = "Combined_Mod5.stan",
               data = dat,
               iter=10000,
               warmup = 2000,
               chains=4,
               control=list("adapt_delta"=0.81),
               pars=c("aMR",
                       "aRSA",
                       "bMassMR",
                       "bMassRSA",
                       "bResid",
                       "resid",
                       "bTemp",
                       "bTherm",
                       "bIntMassTS",
                       "lambda_MR",
                       "lambda_RSA",
                       "sigma_MR",
                       "sigma_RSA",
                       "log_lik",
                       "sigma_total_MR",
                       "sigma_total_RSA"))

## C6_LG ##

dat <- list(N=nrow(Overlap.MR.RSA),
            RSA=Overlap.MR.RSA[,'LogRSAcm2'],
            MRSA=Overlap.MR.RSA[,'LogCenteredMeanMassRSA'],
            MR=Overlap.MR.RSA[, 'LogMeanWholeOrganismMRWatts'],
            MMR=Overlap.MR.RSA[, 'LogCenteredMeanMassMR'],
            Temp = Overlap.MR.RSA[, 'Scaled.Centered.MR.InverseTemp'],
            ThermoStrat = Overlap.MR.RSA$LungsGills,
            d_mat = d_mat,
            vcov_mat = vcov_mat)

C_MOD6_LG <- stan(file = "Combined_Mod6.stan",
               data = dat,
               iter=10000,
               warmup = 2000,
               chains=4,
               control=list("adapt_delta"=0.81),
               pars=c("aMR",
                      "aRSA",
                      "bMassMR",
                      "bMassRSA",
                      "bResid",
                      "resid",
                      "bTemp",
                      "bTherm",
                      "bIntMassTS",
                      "bIntMassResid",
                      "lambda_MR",
                      "lambda_RSA",
                      "sigma_MR",
                      "sigma_RSA",
                      "log_lik",
                      "sigma_total_MR",
                      "sigma_total_RSA"))

## model comparisons for models with respiratory organ instead of thermoregulatory strategy ##

## C3_LG ##

posterior_MR_RSA3 <- as.data.frame(C_MOD3_LG)

residual_rsa_dist <- as.matrix(posterior_MR_RSA3[grep("resid", colnames(posterior_MR_RSA3))])
residualRSA <- colMeans(residual_rsa_dist)
Overlap.MR.RSA$residualRSA <- residualRSA

y_MR_RSA3 <- Overlap.MR.RSA$LogMeanWholeOrganismMRWatts
N_MR_RSA3 <- length(y_MR_RSA3)
S_MR_RSA3 <- nrow(posterior_MR_RSA3)
loglik_MR_RSA3 <- yloo_MR_RSA3 <- sdloo_MR_RSA3 <- matrix(nrow = S_MR_RSA3, 
                                                          ncol = N_MR_RSA3)

st_MR_RSA3 <- as.matrix(posterior_MR_RSA3[grep("sigma_total_MR", 
                                               colnames(posterior_MR_RSA3))])
str(st_MR_RSA3)
dim(st_MR_RSA3)

sigma_array_MR_RSA3 <- array(st_MR_RSA3, dim = c(nrow(st_MR_RSA3), 109, 109))
str(sigma_array_MR_RSA3)
dim(sigma_array_MR_RSA3)

s_MR_RSA3 <- sigma_array_MR_RSA3[2, ,]
s2_MR_RSA3 <- st_MR_RSA3[2, ]
s_MR_RSA3 == s2_MR_RSA3


for (s in 1:S_MR_RSA3) {
  p_MR_RSA3 <- posterior_MR_RSA3[s, ] 
  eta_MR_RSA3 <- p_MR_RSA3$aMR + p_MR_RSA3$bMassMR * Overlap.MR.RSA$LogCenteredMeanMassMR + 
                 p_MR_RSA3$bResid * Overlap.MR.RSA$residualRSA +
                 p_MR_RSA3$bTemp * Overlap.MR.RSA$Scaled.Centered.MR.InverseTemp +
                 p_MR_RSA3$bTherm * Overlap.MR.RSA$ThermoStrat
  Cinv_MR_RSA3 <- solve(sigma_array_MR_RSA3[s, , ]) 
  g_MR_RSA3 <- Cinv_MR_RSA3 %*% (y_MR_RSA3 - eta_MR_RSA3) 
  cbar_MR_RSA3 <- diag(Cinv_MR_RSA3) 
  yloo_MR_RSA3[s, ] <- y_MR_RSA3 - g_MR_RSA3 / cbar_MR_RSA3
  sdloo_MR_RSA3[s, ] <- sqrt(1 / cbar_MR_RSA3) 
  loglik_MR_RSA3[s, ] <- dnorm(y_MR_RSA3, yloo_MR_RSA3[s, ], sdloo_MR_RSA3[s, ], log = TRUE)
}

log_ratios_MR_RSA3 <- -loglik_MR_RSA3
r_eff_MR_RSA3 <- relative_eff(exp(loglik_MR_RSA3), chain_id = rep(1:4, each = 800)) 
psis_result_MR_RSA3 <- psis(log_ratios_MR_RSA3, r_eff = r_eff_MR_RSA3)

plot(psis_result_MR_RSA3, label_points = TRUE)

(psis_loo_MR_RSA3_LG <- loo(loglik_MR_RSA3))


## C4_LG ##

posterior_MR_RSA4 <- as.data.frame(C_MOD4_LG)

residual_rsa_dist <- as.matrix(posterior_MR_RSA4[grep("resid", colnames(posterior_MR_RSA4))])
residualRSA <- colMeans(residual_rsa_dist)
Overlap.MR.RSA$residualRSA <- residualRSA

y_MR_RSA4 <- Overlap.MR.RSA$LogMeanWholeOrganismMRWatts
N_MR_RSA4 <- length(y_MR_RSA4)
S_MR_RSA4 <- nrow(posterior_MR_RSA4)
loglik_MR_RSA4 <- yloo_MR_RSA4 <- sdloo_MR_RSA4 <- matrix(nrow = S_MR_RSA4, 
                                                          ncol = N_MR_RSA4)

st_MR_RSA4 <- as.matrix(posterior_MR_RSA4[grep("sigma_total_MR", 
                                               colnames(posterior_MR_RSA4))])
str(st_MR_RSA4)
dim(st_MR_RSA4)

sigma_array_MR_RSA4 <- array(st_MR_RSA4, dim = c(nrow(st_MR_RSA4), 109, 109))
str(sigma_array_MR_RSA4)
dim(sigma_array_MR_RSA4)

s_MR_RSA4 <- sigma_array_MR_RSA4[2, ,]
s2_MR_RSA4 <- st_MR_RSA4[2, ]
s_MR_RSA4 == s2_MR_RSA4

for (s in 1:S_MR_RSA4) {
  p_MR_RSA4 <- posterior_MR_RSA4[s, ] 
  eta_MR_RSA4 <- p_MR_RSA4$aMR + p_MR_RSA4$bMassMR * Overlap.MR.RSA$LogCenteredMeanMassMR + 
                 p_MR_RSA4$bResid * Overlap.MR.RSA$residualRSA +
                 p_MR_RSA4$bIntMassResid * Overlap.MR.RSA$LogCenteredMeanMassMR * Overlap.MR.RSA$residualRSA +
                 p_MR_RSA4$bTemp * Overlap.MR.RSA$Scaled.Centered.MR.InverseTemp +
                 p_MR_RSA4$bTherm * Overlap.MR.RSA$ThermoStrat
  Cinv_MR_RSA4 <- solve(sigma_array_MR_RSA4[s, , ]) 
  g_MR_RSA4 <- Cinv_MR_RSA4 %*% (y_MR_RSA4 - eta_MR_RSA4) 
  cbar_MR_RSA4 <- diag(Cinv_MR_RSA4) 
  yloo_MR_RSA4[s, ] <- y_MR_RSA4 - g_MR_RSA4 / cbar_MR_RSA4
  sdloo_MR_RSA4[s, ] <- sqrt(1 / cbar_MR_RSA4) 
  loglik_MR_RSA4[s, ] <- dnorm(y_MR_RSA4, yloo_MR_RSA4[s, ], sdloo_MR_RSA4[s, ], log = TRUE)
}

log_ratios_MR_RSA4 <- -loglik_MR_RSA4
r_eff_MR_RSA4 <- relative_eff(exp(loglik_MR_RSA4), chain_id = rep(1:4, each = 800)) 
psis_result_MR_RSA4 <- psis(log_ratios_MR_RSA4, r_eff = r_eff_MR_RSA4)

plot(psis_result_MR_RSA4, label_points = TRUE)

(psis_loo_MR_RSA4_LG <- loo(loglik_MR_RSA4))


## C5_LG ##

posterior_MR_RSA5 <- as.data.frame(C_MOD5_LG)

residual_rsa_dist <- as.matrix(posterior_MR_RSA5[grep("resid", colnames(posterior_MR_RSA5))])
residualRSA <- colMeans(residual_rsa_dist)
Overlap.MR.RSA$residualRSA <- residualRSA

y_MR_RSA5 <- Overlap.MR.RSA$LogMeanWholeOrganismMRWatts
N_MR_RSA5 <- length(y_MR_RSA5)
S_MR_RSA5 <- nrow(posterior_MR_RSA5)
loglik_MR_RSA5 <- yloo_MR_RSA5 <- sdloo_MR_RSA5 <- matrix(nrow = S_MR_RSA5, 
                                                          ncol = N_MR_RSA5)

st_MR_RSA5 <- as.matrix(posterior_MR_RSA5[grep("sigma_total_MR", 
                                               colnames(posterior_MR_RSA5))])
str(st_MR_RSA5)
dim(st_MR_RSA5)

sigma_array_MR_RSA5 <- array(st_MR_RSA5, dim = c(nrow(st_MR_RSA5), 109, 109))
str(sigma_array_MR_RSA5)
dim(sigma_array_MR_RSA5)

s1_MR_RSA5 <- sigma_array_MR_RSA5[2, ,]
s2_MR_RSA5 <- st_MR_RSA5[2, ]
s1_MR_RSA5 == s2_MR_RSA5

for (s in 1:S_MR_RSA5) {
  p_MR_RSA5 <- posterior_MR_RSA5[s, ] 
  eta_MR_RSA5 <- p_MR_RSA5$aMR + 
                 p_MR_RSA5$bResid * Overlap.MR.RSA$residualRSA + 
                 p_MR_RSA5$bMassMR * Overlap.MR.RSA$LogCenteredMeanMassMR + 
                 p_MR_RSA5$bTherm * Overlap.MR.RSA$ThermoStrat  + 
                 p_MR_RSA5$bIntMassTS * Overlap.MR.RSA$LogCenteredMeanMassMR * Overlap.MR.RSA$ThermoStrat + 
                 p_MR_RSA5$bTemp * Overlap.MR.RSA$Scaled.Centered.MR.InverseTemp
  Cinv_MR_RSA5 <- solve(sigma_array_MR_RSA5[s, , ]) 
  g_MR_RSA5 <- Cinv_MR_RSA5 %*% (y_MR_RSA5 - eta_MR_RSA5) 
  cbar_MR_RSA5 <- diag(Cinv_MR_RSA5) 
  yloo_MR_RSA5[s, ] <- y_MR_RSA5 - g_MR_RSA5 / cbar_MR_RSA5 
  sdloo_MR_RSA5[s, ] <- sqrt(1 / cbar_MR_RSA5) 
  loglik_MR_RSA5[s, ] <- dnorm(y_MR_RSA5, yloo_MR_RSA5[s, ], sdloo_MR_RSA5[s, ], log = TRUE)
}

log_ratios_MR_RSA5 <- -loglik_MR_RSA5
r_eff_MR_RSA5 <- relative_eff(exp(loglik_MR_RSA5), chain_id = rep(1:4, each = 800)) 
psis_result_MR_RSA5 <- psis(log_ratios_MR_RSA5, r_eff = r_eff_MR_RSA5)

plot(psis_result_MR_RSA5, label_points = TRUE)

(psis_loo_MR_RSA5_LG <- loo(loglik_MR_RSA5))


## C6_LG ##

posterior_MR_RSA6 <- as.data.frame(C_MOD6_LG)

residual_rsa_dist <- as.matrix(posterior_MR_RSA6[grep("resid", colnames(posterior_MR_RSA6))])
residualRSA <- colMeans(residual_rsa_dist)
Overlap.MR.RSA$residualRSA <- residualRSA

y_MR_RSA6 <- Overlap.MR.RSA$LogMeanWholeOrganismMRWatts
N_MR_RSA6 <- length(y_MR_RSA6)
S_MR_RSA6 <- nrow(posterior_MR_RSA6)
loglik_MR_RSA6 <- yloo_MR_RSA6 <- sdloo_MR_RSA6 <- matrix(nrow = S_MR_RSA6, 
                                                          ncol = N_MR_RSA6)

st_MR_RSA6 <- as.matrix(posterior_MR_RSA6[grep("sigma_total_MR", 
                                               colnames(posterior_MR_RSA6))])
str(st_MR_RSA6)
dim(st_MR_RSA6)

sigma_array_MR_RSA6 <- array(st_MR_RSA6, dim = c(nrow(st_MR_RSA6), 109, 109))
str(sigma_array_MR_RSA6)
dim(sigma_array_MR_RSA6)

s_MR_RSA6 <- sigma_array_MR_RSA6[2, ,]
s2_MR_RSA6 <- st_MR_RSA6[2, ]
s_MR_RSA6 == s2_MR_RSA6

for (s in 1:S_MR_RSA6) {
  p_MR_RSA6 <- posterior_MR_RSA6[s, ] 
  eta_MR_RSA6 <- p_MR_RSA6$aMR +
                 p_MR_RSA6$bResid * Overlap.MR.RSA$residualRSA +
                 p_MR_RSA6$bMassMR * Overlap.MR.RSA$LogCenteredMeanMassMR + 
                 p_MR_RSA6$bTherm * Overlap.MR.RSA$ThermoStrat +
                 p_MR_RSA6$bTemp * Overlap.MR.RSA$Scaled.Centered.MR.InverseTemp +
                 p_MR_RSA6$bIntMassResid * Overlap.MR.RSA$LogCenteredMeanMassMR * Overlap.MR.RSA$residualRSA +
                 p_MR_RSA6$bIntMassTS * Overlap.MR.RSA$LogCenteredMeanMassMR * Overlap.MR.RSA$ThermoStrat
  Cinv_MR_RSA6 <- solve(sigma_array_MR_RSA6[s, , ]) 
  g_MR_RSA6 <- Cinv_MR_RSA6 %*% (y_MR_RSA6 - eta_MR_RSA6) 
  cbar_MR_RSA6 <- diag(Cinv_MR_RSA6) 
  yloo_MR_RSA6[s, ] <- y_MR_RSA6 - g_MR_RSA6 / cbar_MR_RSA6
  sdloo_MR_RSA6[s, ] <- sqrt(1 / cbar_MR_RSA6) 
  loglik_MR_RSA6[s, ] <- dnorm(y_MR_RSA6, yloo_MR_RSA6[s, ], sdloo_MR_RSA6[s, ], log = TRUE)
}

log_ratios_MR_RSA6 <- -loglik_MR_RSA6
r_eff_MR_RSA6 <- relative_eff(exp(loglik_MR_RSA6), chain_id = rep(1:4, each = 800)) 
psis_result_MR_RSA6 <- psis(log_ratios_MR_RSA6, r_eff = r_eff_MR_RSA6)

plot(psis_result_MR_RSA6, label_points = TRUE)

(psis_loo_MR_RSA6_LG <- loo(loglik_MR_RSA6))

## comparisons 

Combined_model_comparison <- compare(psis_loo_MR_RSA1, psis_loo_MR_RSA2,
                                     psis_loo_MR_RSA3, psis_loo_MR_RSA4,
                                     psis_loo_MR_RSA5, psis_loo_MR_RSA6,
                                     psis_loo_MR_RSA3_LG, psis_loo_MR_RSA4_LG,
                                     psis_loo_MR_RSA5_LG, psis_loo_MR_RSA6_LG)

loo_list_combined <- list(psis_loo_MR_RSA1, psis_loo_MR_RSA2,
                          psis_loo_MR_RSA3, psis_loo_MR_RSA4,
                          psis_loo_MR_RSA5, psis_loo_MR_RSA6,
                          psis_loo_MR_RSA3_LG, psis_loo_MR_RSA4_LG,
                          psis_loo_MR_RSA5_LG, psis_loo_MR_RSA6_LG)

combined_mod_wt <- loo_model_weights(loo_list_combined)


MR_RSA3_mod_comparison <- loo_compare(psis_loo_MR_RSA3, psis_loo_MR_RSA3_LG)

MR_RSA4_mod_comparison <- loo_compare(psis_loo_MR_RSA4, psis_loo_MR_RSA4_LG)

MR_RSA5_mod_comparison <- loo_compare(psis_loo_MR_RSA5, psis_loo_MR_RSA5_LG)

MR_RSA6_mod_comparison <- loo_compare(psis_loo_MR_RSA6, psis_loo_MR_RSA6_LG)



## model with standardized predictors for comparing relative effect sizes ##

dat_STD <- list(N=nrow(Overlap.MR.RSA),
           RSA=Overlap.MR.RSA[,'LogRSAcm2'],
           MRSA=Overlap.MR.RSA[,'RSA.LogMeanMassG'],
           MR=Overlap.MR.RSA[, 'LogMeanWholeOrganismMRWatts'],
           MMR=Overlap.MR.RSA[, 'MR.LogMeanMassG'],
           Temp = Overlap.MR.RSA[, 'MR.InverseTemp'],
           ThermoStrat = Overlap.MR.RSA[, "ThermoStrat"],
           d_mat = d_mat,
           vcov_mat = vcov_mat)

# model C5 with z-score standardization

C_MOD5_1SD <- stan(file = here("C_MOD5_1SD.stan"),
                 data = dat_STD,
                 iter=10000,
                 warmup = 2000,
                 thin = 10,
                 chains=4,
                 control=list("adapt_delta"=0.81),
                 pars=c( "aMR",
                         "std_bMassMR",
                         "std_bMassRSA",
                         "std_bResid",
                         "std_resid",
                         "std_bTemp",
                         "bTherm",
                         "bIntMassTS",
                         "lambda_MR",
                 				 "std_sigma_MR"))


