# Bigman et al. metabolic rate and respiratory surface area scaling Science Advances

# 03 - 'Respiratory surface area models'

library(rstan)
library(beepr)
library(loo)

# set options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# create matrix for sigma
d_mat <- diag(1, 109, 109)


## model names correspond to Table S1 ##

## RSA1 ##
dat <- list(N = nrow(Overlap.MR.RSA),
            RSA = Overlap.MR.RSA[,'LogRSAcm2'],
            MRSA = Overlap.MR.RSA[,'LogCenteredMeanMassRSA'],
            d_mat = d_mat,
            vcov_mat = vcov_mat)

RSA_MOD1 <- stan(file = "RSA_Mod1.stan",
                                data = dat,
                                iter=5000,
                                warmup = 1000,
                                thin = 10,
                                chains=4,
                                control=list("adapt_delta"=0.81),
                                pars=c("a",
                                       "bMass",
                                       "resid",
                                       "lambda",
                                       "log_lik",
                                       "sigma",
                                       "sigma_total"))


## RSA2 ##

data_matrix_mod2 <- model.matrix( ~ LogCenteredMeanMassRSA + ThermoStrat,
                                    data = Overlap.MR.RSA)

dat <- list(N = nrow(Overlap.MR.RSA),
            RSA = Overlap.MR.RSA[,'LogRSAcm2'],
            data_matrix = data_matrix_mod2,
            K = ncol(data_matrix_mod2),
            d_mat = d_mat,
            vcov_mat = vcov_mat)

RSA_MOD2 <- stan(file = "RSA_data_matrix_mod.stan",
                 data = dat,
                 iter= 5000,
                 warmup = 1000,
                 thin = 10, 
                 chains=4,
                 control=list("adapt_delta"=0.81),
                 pars=c("beta",
                        "resid",
                        "lambda",
                        "sigma",
                        "sigma_total",
                        "log_lik"))

## RSA3 ##

data_matrix_mod3 <- model.matrix( ~ LogCenteredMeanMassRSA * ThermoStrat,
                                    data = Overlap.MR.RSA)


dat <- list(N = nrow(Overlap.MR.RSA),
            RSA = Overlap.MR.RSA[,'LogRSAcm2'],
            data_matrix = data_matrix_mod3,
            K = ncol(data_matrix_mod3),
            d_mat = d_mat,
            vcov_mat = vcov_mat)

RSA_MOD3 <- stan(file = "RSA_data_matrix_mod.stan",
                 data = dat,
                 iter= 5000,
                 warmup = 1000,
                 thin = 10,
                 chains=4,
                 control=list("adapt_delta"=0.81),
                 pars=c("beta",
                        "resid",
                        "lambda",
                        "log_lik",
                        "sigma_total",
                        "sigma"))



#### model comparison ####

## RSA1 ##

posterior_RSA1 <- as.data.frame(RSA_MOD1)
y_RSA1 <- Overlap.MR.RSA$LogRSAcm2
N_RSA1 <- length(y_RSA1)
S_RSA1 <- nrow(posterior_RSA1)
loglik_RSA1 <- yloo_RSA1 <- sdloo_RSA1 <- matrix(nrow = S_RSA1, ncol = N_RSA1)

st_RSA1 <- as.matrix(posterior_RSA1[grep("sigma_total", colnames(posterior_RSA1))])
str(st_RSA1)
dim(st_RSA1)

sigma_array_RSA1 <- array(st_RSA1, dim = c(nrow(st_RSA1), 109, 109))
str(sigma_array_RSA1)
dim(sigma_array_RSA1)

s1_RSA1 <- sigma_array_RSA1[2, ,]
s2_RSA1 <- st_RSA1[2, ]
s1_RSA1 == s2_RSA1

for (s in 1:S_RSA1) {
  p_RSA1 <- posterior_RSA1[s, ] 
  eta_RSA1 <- p_RSA1$a + p_RSA1$bMass * Overlap.MR.RSA$LogCenteredMeanMassRSA 
  Cinv_RSA1 <- solve(sigma_array_RSA1[s, , ]) 
  g_RSA1 <- Cinv_RSA1 %*% (y_RSA1 - eta_RSA1) 
  cbar_RSA1 <- diag(Cinv_RSA1) 
  yloo_RSA1[s, ] <- y_RSA1 - g_RSA1 / cbar_RSA1 
  sdloo_RSA1[s, ] <- sqrt(1 / cbar_RSA1) 
  loglik_RSA1[s, ] <- dnorm(y_RSA1, yloo_RSA1[s, ], sdloo_RSA1[s, ], log = TRUE)
}

log_ratios_RSA1 <- -loglik_RSA1
r_eff_RSA1 <- relative_eff(exp(loglik_RSA1), chain_id = rep(1:4, each = 400)) 
psis_result_RSA1 <- psis(log_ratios_RSA1, r_eff = r_eff_RSA1)

plot(psis_result_RSA1, label_points = TRUE)

(psis_loo_RSA1 <- loo(loglik_RSA1))


## RSA2 ##

posterior_RSA2 <- as.data.frame(RSA_MOD2)
posterior_RSA2_short <- posterior_RSA2 %>% 
                        dplyr::select("beta[1]", "beta[2]", "beta[3]") %>%
                        rename(a = "beta[1]",
                               bMass = "beta[2]",
                               bTherm = "beta[3]")

y_RSA2 <- Overlap.MR.RSA$LogRSAcm2
N_RSA2 <- length(y_RSA2)
S_RSA2 <- nrow(posterior_RSA2)
loglik_RSA2 <- yloo_RSA2 <- sdloo_RSA2 <- matrix(nrow = S_RSA2, ncol = N_RSA2)

st_RSA2 <- as.matrix(posterior_RSA2[grep("sigma_total", colnames(posterior_RSA2))])
str(st_RSA2)
dim(st_RSA2)

sigma_array_RSA2 <- array(st_RSA2, dim = c(nrow(st_RSA2), 109, 109))
str(sigma_array_RSA2)
dim(sigma_array_RSA2)

s1_RSA2 <- sigma_array_RSA2[2, ,]
s2_RSA2 <- st_RSA2[2, ]
s1_RSA2 == s2_RSA2

for (s in 1:S_RSA2) {
  p_RSA2 <- posterior_RSA2_short[s, ] 
  eta_RSA2 <- p_RSA2$a + p_RSA2$bMass * Overlap.MR.RSA$LogCenteredMeanMassRSA +
              p_RSA2$bTherm * Overlap.MR.RSA$ThermoStrat
  Cinv_RSA2 <- solve(sigma_array_RSA2[s, , ]) 
  g_RSA2 <- Cinv_RSA2 %*% (y_RSA2 - eta_RSA2) 
  cbar_RSA2 <- diag(Cinv_RSA2) 
  yloo_RSA2[s, ] <- y_RSA2 - g_RSA2 / cbar_RSA2 
  sdloo_RSA2[s, ] <- sqrt(1 / cbar_RSA2) 
  loglik_RSA2[s, ] <- dnorm(y_RSA2, yloo_RSA2[s, ], sdloo_RSA2[s, ], log = TRUE)
}

log_ratios_RSA2 <- -loglik_RSA2
r_eff_RSA2 <- relative_eff(exp(loglik_RSA2), chain_id = rep(1:4, each = 400)) 
psis_result_RSA2 <- psis(log_ratios_RSA2, r_eff = r_eff_RSA2)

plot(psis_result_RSA2, label_points = TRUE)

(psis_loo_RSA2 <- loo(loglik_RSA2))


## RSA3 ##

posterior_RSA3 <- as.data.frame(RSA_MOD3)
posterior_RSA3_short <- posterior_RSA3 %>% 
                        dplyr::select("beta[1]", "beta[2]", "beta[3]", "beta[4]") %>%
                        rename(a = "beta[1]",
                               bMass = "beta[2]",
                               bTherm = "beta[3]",
                               bMass_Therm = "beta[4]")

y_RSA3 <- Overlap.MR.RSA$LogRSAcm2
N_RSA3 <- length(y_RSA3)
S_RSA3 <- nrow(posterior_RSA3)
loglik_RSA3 <- yloo_RSA3 <- sdloo_RSA3 <- matrix(nrow = S_RSA3, ncol = N_RSA3)

st_RSA3 <- as.matrix(posterior_RSA3[grep("sigma_total", colnames(posterior_RSA3))])
str(st_RSA3)
dim(st_RSA3)

sigma_array_RSA3 <- array(st_RSA3, dim = c(nrow(st_RSA3), 109, 109))
str(sigma_array_RSA3)
dim(sigma_array_RSA3)

s1_RSA3 <- sigma_array_RSA3[2, ,]
s2_RSA3 <- st_RSA3[2, ]
s1_RSA3 == s2_RSA3

for (s in 1:S_RSA3) {
  p_RSA3 <- posterior_RSA3_short[s, ] 
  eta_RSA3 <- p_RSA3$a + p_RSA3$bMass * Overlap.MR.RSA$LogCenteredMeanMassRSA +
              p_RSA3$bTherm * Overlap.MR.RSA$ThermoStrat +
              p_RSA3$bMass_Therm * Overlap.MR.RSA$LogCenteredMeanMassRSA * 
                  Overlap.MR.RSA$ThermoStrat
  Cinv_RSA3 <- solve(sigma_array_RSA3[s, , ]) 
  g_RSA3 <- Cinv_RSA3 %*% (y_RSA3 - eta_RSA3) 
  cbar_RSA3 <- diag(Cinv_RSA3) 
  yloo_RSA3[s, ] <- y_RSA3 - g_RSA3 / cbar_RSA3 
  sdloo_RSA3[s, ] <- sqrt(1 / cbar_RSA3) 
  loglik_RSA3[s, ] <- dnorm(y_RSA3, yloo_RSA3[s, ], sdloo_RSA3[s, ], log = TRUE)
}

log_ratios_RSA3 <- -loglik_RSA3
r_eff_RSA3 <- relative_eff(exp(loglik_RSA3), chain_id = rep(1:4, each = 400)) 
psis_result_RSA3 <- psis(log_ratios_RSA3, r_eff = r_eff_RSA3)

plot(psis_result_RSA3, label_points = TRUE)

(psis_loo_RSA3 <- loo(loglik_RSA3))


## compare models ##
RSA_model_comparison <- loo_compare(psis_loo_RSA1, psis_loo_RSA2, psis_loo_RSA3)


# weights
loo_list_RSA_only <- list(psis_loo_RSA1, psis_loo_RSA2, psis_loo_RSA3)
RSA_only_mod_wt <- loo_model_weights(loo_list_RSA_only)


#### models with respiratory organ in place of thermoregulatory strategy ####

## RSA2_LG ##

data_matrix_mod2_LG <- model.matrix( ~ LogCenteredMeanMassRSA + LungsGills,
                                       data = Overlap.MR.RSA)

dat <- list(N = nrow(Overlap.MR.RSA),
            RSA = Overlap.MR.RSA[,'LogRSAcm2'],
            data_matrix = data_matrix_mod2_LG,
            K = ncol(data_matrix_mod2_LG),
            d_mat = d_mat,
            vcov_mat = vcov_mat)

RSA_MOD2_LG <- stan(file = "RSA_data_matrix_mod.stan",
                 data = dat,
                 iter=10000,
                 warmup = 2000,
                 chains=4,
                 control=list("adapt_delta"=0.81),
                 pars=c("beta",
                        "resid",
                        "lambda",
                        "sigma",
                        "sigma_total",
                        "log_lik"))

## RSA3_LG  ##
data_matrix_mod3_LG <- model.matrix( ~ LogCenteredMeanMassRSA * LungsGills,
                                    data = Overlap.MR.RSA)

dat <- list(N = nrow(Overlap.MR.RSA),
            RSA = Overlap.MR.RSA[,'LogRSAcm2'],
            data_matrix = data_matrix_mod3_LG,
            K = ncol(data_matrix_mod3_LG),
            d_mat = d_mat,
            vcov_mat = vcov_mat)

RSA_MOD3_LG <- stan(file = "RSA_data_matrix_mod.stan",
                 data = dat,
                 iter=10000,
                 warmup = 2000,
                 chains=4,
                 control=list("adapt_delta"=0.81),
                 pars=c("beta",
                        "resid",
                        "lambda",
                        "log_lik",
                        "sigma_total",
                        "sigma"))

#### model comparison ####

## RSA2_LG ##

posterior_RSA2 <- Thin(as.data.frame(RSA_MOD2_LG), By = 10) %>% as.data.frame()
posterior_RSA2_short <- posterior_RSA2 %>% 
                        dplyr::select("beta[1]", "beta[2]", "beta[3]") %>%
                        rename(a = "beta[1]",
                               bMass = "beta[2]",
                               bTherm = "beta[3]")

y_RSA2 <- Overlap.MR.RSA$LogRSAcm2
N_RSA2 <- length(y_RSA2)
S_RSA2 <- nrow(posterior_RSA2)
loglik_RSA2 <- yloo_RSA2 <- sdloo_RSA2 <- matrix(nrow = S_RSA2, ncol = N_RSA2)

st_RSA2 <- as.matrix(posterior_RSA2[grep("sigma_total", colnames(posterior_RSA2))])
str(st_RSA2)
dim(st_RSA2)

sigma_array_RSA2 <- array(st_RSA2, dim = c(nrow(st_RSA2), 109, 109))
str(sigma_array_RSA2)
dim(sigma_array_RSA2)

s1_RSA2 <- sigma_array_RSA2[2, ,]
s2_RSA2 <- st_RSA2[2, ]
s1_RSA2 == s2_RSA2

for (s in 1:S_RSA2) {
  p_RSA2 <- posterior_RSA2_short[s, ] 
  eta_RSA2 <- p_RSA2$a + p_RSA2$bMass * Overlap.MR.RSA$LogCenteredMeanMassRSA +
              p_RSA2$bTherm * Overlap.MR.RSA$ThermoStrat
  Cinv_RSA2 <- solve(sigma_array_RSA2[s, , ]) 
  g_RSA2 <- Cinv_RSA2 %*% (y_RSA2 - eta_RSA2) 
  cbar_RSA2 <- diag(Cinv_RSA2) 
  yloo_RSA2[s, ] <- y_RSA2 - g_RSA2 / cbar_RSA2 
  sdloo_RSA2[s, ] <- sqrt(1 / cbar_RSA2) 
  loglik_RSA2[s, ] <- dnorm(y_RSA2, yloo_RSA2[s, ], sdloo_RSA2[s, ], log = TRUE)
}

log_ratios_RSA2 <- -loglik_RSA2
r_eff_RSA2 <- relative_eff(exp(loglik_RSA2), chain_id = rep(1:4, each = 800)) 
psis_result_RSA2 <- psis(log_ratios_RSA2, r_eff = r_eff_RSA2)

plot(psis_result_RSA2, label_points = TRUE)

(psis_loo_RSA2_LG <- loo(loglik_RSA2))


## RSA3_LG ##

posterior_RSA3 <- Thin(as.data.frame(RSA_MOD3_LG), By = 10) %>% as.data.frame()

posterior_RSA3_short <- posterior_RSA3 %>% 
                        dplyr::select("beta[1]", "beta[2]", "beta[3]", "beta[4]") %>%
                        rename(a = "beta[1]",
                               bMass = "beta[2]",
                               bTherm = "beta[3]",
                               bMass_Therm = "beta[4]")

y_RSA3 <- Overlap.MR.RSA$LogRSAcm2
N_RSA3 <- length(y_RSA3)
S_RSA3 <- nrow(posterior_RSA3)
loglik_RSA3 <- yloo_RSA3 <- sdloo_RSA3 <- matrix(nrow = S_RSA3, ncol = N_RSA3)

st_RSA3 <- as.matrix(posterior_RSA3[grep("sigma_total", colnames(posterior_RSA3))])
str(st_RSA3)
dim(st_RSA3)

sigma_array_RSA3 <- array(st_RSA3, dim = c(nrow(st_RSA3), 109, 109))
str(sigma_array_RSA3)
dim(sigma_array_RSA3)

s1_RSA3 <- sigma_array_RSA3[2, ,]
s2_RSA3 <- st_RSA3[2, ]
s1_RSA3 == s2_RSA3

for (s in 1:S_RSA3) {
  p_RSA3 <- posterior_RSA3_short[s, ] 
  eta_RSA3 <- p_RSA3$a + p_RSA3$bMass * Overlap.MR.RSA$LogCenteredMeanMassRSA +
              p_RSA3$bTherm * Overlap.MR.RSA$ThermoStrat +
              p_RSA3$bMass_Therm * Overlap.MR.RSA$LogCenteredMeanMassRSA * 
                  Overlap.MR.RSA$ThermoStrat
  Cinv_RSA3 <- solve(sigma_array_RSA3[s, , ]) 
  g_RSA3 <- Cinv_RSA3 %*% (y_RSA3 - eta_RSA3) 
  cbar_RSA3 <- diag(Cinv_RSA3) 
  yloo_RSA3[s, ] <- y_RSA3 - g_RSA3 / cbar_RSA3 
  sdloo_RSA3[s, ] <- sqrt(1 / cbar_RSA3) 
  loglik_RSA3[s, ] <- dnorm(y_RSA3, yloo_RSA3[s, ], sdloo_RSA3[s, ], log = TRUE)
}

log_ratios_RSA3 <- -loglik_RSA3
r_eff_RSA3 <- relative_eff(exp(loglik_RSA3), chain_id = rep(1:4, each = 800)) 
psis_result_RSA3 <- psis(log_ratios_RSA3, r_eff = r_eff_RSA3)

plot(psis_result_RSA3, label_points = TRUE)

(psis_loo_RSA3_LG <- loo(loglik_RSA3))


RSA_model_comparison <- compare(psis_loo_RSA1, psis_loo_RSA2, psis_loo_RSA3,
                                psis_loo_RSA2_LG, psis_loo_RSA3_LG)

# weights
loo_list_RSA_only <- list(psis_loo_RSA1, psis_loo_RSA2, psis_loo_RSA3,
                                psis_loo_RSA2_LG, psis_loo_RSA3_LG)
RSA_only_mod_wt <- loo_model_weights(loo_list_RSA_only)

RSA2_mod_comparison <- compare(psis_loo_RSA2, psis_loo_RSA2_LG)
RSA3_mod_comparison <- compare(psis_loo_RSA3, psis_loo_RSA3_LG)
