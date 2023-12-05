# Bigman et al. metabolic rate and respiratory surface area scaling Science Advances

# 02 - 'Metabolic rate models'

## set stan options ##
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## create identity matrix for sigma ##
d_mat <- diag(1, 109, 109)

#### models, names correspond to Table S1 ####

## MR1 ##

dat <- list(N=nrow(Overlap.MR.RSA),
            MR=Overlap.MR.RSA[, 'LogMeanWholeOrganismMRWatts'],
            MMR=Overlap.MR.RSA[, 'LogCenteredMeanMassMR'],
            InTemp=Overlap.MR.RSA[, 'Scaled.Centered.MR.InverseTemp'],
            d_mat = d_mat,
            vcov_mat = vcov_mat)

MR_MOD1 <- stan(file = here("MR_MOD1.stan"),
                data = dat,
                iter = 5000,
                warmup = 1000,
                thin = 10,
                chains = 4,
                control = list("adapt_delta" = 0.81),
                pars = c("a",
                      "bMass",
                      "bTemp",
                      "lambda",
                      "log_lik",
                      "sigma",
                      "sigma_total"))

## MR2 ##

data_matrix_mod2 <- model.matrix( ~ LogCenteredMeanMassMR +
                                    Scaled.Centered.MR.InverseTemp +
                                    ThermoStrat,
                                    data = Overlap.MR.RSA)

dat <- list(N=nrow(Overlap.MR.RSA),
            MR=Overlap.MR.RSA[, 'LogMeanWholeOrganismMRWatts'],
            data_matrix = data_matrix_mod2,
            K = ncol(data_matrix_mod2),
            d_mat = d_mat,
            vcov_mat = vcov_mat)

MR_MOD2 <- stan(file =  here("MR_data_matrix_mod.stan"),
                data = dat,
                iter = 5000,
                warmup = 1000,
                chains = 4,
                thin = 10,
                control = list("adapt_delta" = 0.99,max_treedepth = 18),
                pars = c("beta",
                         "lambda",
                         "sigma",
                         "log_lik",
                         "sigma_total"))


## MR3 ##

data_matrix_mod3 <- model.matrix( ~ LogCenteredMeanMassMR *
                                    ThermoStrat +
                                    Scaled.Centered.MR.InverseTemp,
                                    data = Overlap.MR.RSA)

dat <- list(N=nrow(Overlap.MR.RSA),
            MR=Overlap.MR.RSA[, 'LogMeanWholeOrganismMRWatts'],
            data_matrix = data_matrix_mod3,
            K = ncol(data_matrix_mod3),
            d_mat = d_mat,
            vcov_mat = vcov_mat)

MR_MOD3 <- stan(file = here("MR_data_matrix_mod.stan"),
                  data=dat,
                  iter=5000,
                  warmup = 1000,
                  chains=4,
                  thin = 10,
                  control=list("adapt_delta"=0.999, max_treedepth = 18),
                  pars=c("beta",
                         "lambda",
                         "sigma",
                         "log_lik",
                         "mu_MR",
                         "sigma_total"))


#### model comparison ####

## MR1 ##

posterior_MR1 <- as.data.frame(MR_MOD1)
y_MR1 <- Overlap.MR.RSA$LogCenteredMeanMassMR
N_MR1 <- length(y_MR1)
S_MR1 <- nrow(posterior_MR1)
loglik_MR1 <- yloo_MR1 <- sdloo_MR1 <- matrix(nrow = S_MR1, ncol = N_MR1)

st_MR1 <- as.matrix(posterior_MR1[grep("sigma_total", colnames(posterior_MR1))])
str(st_MR1)
dim(st_MR1)

sigma_array_MR1 <- array(st_MR1, dim = c(nrow(st_MR1), 109, 109))
str(sigma_array_MR1)
dim(sigma_array_MR1)

s_MR1 <- sigma_array_MR1[2, ,]
s2_MR1 <- st_MR1[2, ]
s_MR1 == s2_MR1

for (s in 1:S_MR1) {
  p_MR1 <- posterior_MR1[s, ]
  eta_MR1 <- p_MR1$a + p_MR1$bMass * Overlap.MR.RSA$LogCenteredMeanMassMR +
             p_MR1$bTemp * Overlap.MR.RSA$Scaled.Centered.MR.InverseTemp
  Cinv_MR1 <- solve(sigma_array_MR1[s, , ])
  g_MR1 <- Cinv_MR1 %*% (y_MR1 - eta_MR1)
  cbar_MR1 <- diag(Cinv_MR1)
  yloo_MR1[s, ] <- y_MR1 - g_MR1 / cbar_MR1
  sdloo_MR1[s, ] <- sqrt(1 / cbar_MR1)
  loglik_MR1[s, ] <- dnorm(y_MR1, yloo_MR1[s, ], sdloo_MR1[s, ], log = TRUE)
}

log_ratios_MR1 <- -loglik_MR1
r_eff_MR1 <- relative_eff(exp(loglik_MR1), chain_id = rep(1:4, each = 400))
psis_result_MR1 <- psis(log_ratios_MR1, r_eff = r_eff_MR1)

plot(psis_result_MR1, label_points = TRUE)

(psis_loo_MR1 <- loo(loglik_MR1))


## MR2 ##

posterior_MR2 <- as.data.frame(MR_MOD2)
posterior_MR2_short <- posterior_MR2 %>% dplyr::select("beta[1]", "beta[2]",
                                          "beta[3]", "beta[4]") %>%
                                   rename(a = "beta[1]",
                                          bMass = "beta[2]",
                                          bTemp = "beta[3]",
                                          bTherm = "beta[4]"
                                          )

y_MR2 <- Overlap.MR.RSA$LogMeanWholeOrganismMRWatts
N_MR2 <- length(y_MR2)
S_MR2 <- nrow(posterior_MR2)
loglik_MR2 <- yloo_MR2 <- sdloo_MR2 <- matrix(nrow = S_MR2, ncol = N_MR2)

st_MR2 <- as.matrix(posterior_MR2[grep("sigma_total", colnames(posterior_MR2))])
str(st_MR2)
dim(st_MR2)

sigma_array_MR2 <- array(st_MR2, dim = c(nrow(st_MR2), 109, 109))
str(sigma_array_MR2)
dim(sigma_array_MR2)

s1_MR2 <- sigma_array_MR2[2, ,]
s2_MR2 <- st_MR2[2, ]
s1_MR2 == s2_MR2

for (s in 1:S_MR2) {
  p_MR2 <- posterior_MR2_short[s, ]
  eta_MR2 <- p_MR2$a + p_MR2$bMass * Overlap.MR.RSA$LogCenteredMeanMassMR +
             p_MR2$bTemp * Overlap.MR.RSA$Scaled.Centered.MR.InverseTemp +
             p_MR2$bTherm * Overlap.MR.RSA$ThermoStrat
  Cinv_MR2 <- solve(sigma_array_MR2[s, , ])
  g_MR2 <- Cinv_MR2 %*% (y_MR2 - eta_MR2)
  cbar_MR2 <- diag(Cinv_MR2)
  yloo_MR2[s, ] <- y_MR2 - g_MR2 / cbar_MR2
  sdloo_MR2[s, ] <- sqrt(1 / cbar_MR2)
  loglik_MR2[s, ] <- dnorm(y_MR2, yloo_MR2[s, ], sdloo_MR2[s, ], log = TRUE)
}

log_ratios_MR2 <- -loglik_MR2
r_eff_MR2 <- relative_eff(exp(loglik_MR2), chain_id = rep(1:4, each = 400))
psis_result_MR2  <- psis(log_ratios_MR2, r_eff = r_eff_MR2)

plot(psis_result_MR2, label_points = TRUE)

(psis_loo_MR2 <- loo(loglik_MR2))


## MR3 ##

posterior_MR3 <- as.data.frame(MR_MOD3)
posterior_MR3_short <- posterior_MR3 %>% dplyr::select(
                                "beta[1]", "beta[2]", "beta[3]", "beta[4]",
                                "beta[5]") %>%
                                rename(a = "beta[1]",
                                          bMass = "beta[2]",
                                          bTherm = "beta[3]",
                                          bTemp = "beta[4]",
                                          bMass_Therm = "beta[5]"
                                          )

y_MR3 <- Overlap.MR.RSA$LogMeanWholeOrganismMRWatts
N_MR3 <- length(y_MR3)
S_MR3 <- nrow(posterior_MR3)
loglik_MR3 <- yloo_MR3 <- sdloo_MR3 <- matrix(nrow = S_MR3, ncol = N_MR3)

st_MR3 <- as.matrix(posterior_MR3[grep("sigma_total", colnames(posterior_MR3))])
str(st_MR3)
dim(st_MR3)

sigma_array_MR3 <- array(st_MR3, dim = c(nrow(st_MR3), 109, 109))
str(sigma_array_MR3)
dim(sigma_array_MR3)

s1_MR3 <- sigma_array_MR3[2, ,]
s2_MR3 <- st_MR3[2, ]
s1_MR3 == s2_MR3

for (s in 1:S_MR3) {
  p_MR3 <- posterior_MR3_short[s, ]
  eta_MR3 <- p_MR3$a + p_MR3$bMass * Overlap.MR.RSA$LogCenteredMeanMassMR  +
             p_MR3$bTherm * Overlap.MR.RSA$ThermoStrat +
             p_MR3$bTemp * Overlap.MR.RSA$Scaled.Centered.MR.InverseTemp +
             p_MR3$bMass_Therm * Overlap.MR.RSA$LogCenteredMeanMassMR *
                 Overlap.MR.RSA$ThermoStrat
  Cinv_MR3 <- solve(sigma_array_MR3[s, , ])
  g_MR3 <- Cinv_MR3 %*% (y_MR3 - eta_MR3)
  cbar_MR3 <- diag(Cinv_MR3)
  yloo_MR3[s, ] <- y_MR3 - g_MR3 / cbar_MR3
  sdloo_MR3[s, ] <- sqrt(1 / cbar_MR3)
  loglik_MR3[s, ] <- dnorm(y_MR3, yloo_MR3[s, ], sdloo_MR3[s, ], log = TRUE)
}

log_ratios_MR3 <- -loglik_MR3
r_eff_MR3 <- relative_eff(exp(loglik_MR3), chain_id = rep(1:4, each = 400))
psis_result_MR3 <- psis(log_ratios_MR3, r_eff = r_eff_MR3)

plot(psis_result_MR3, label_points = TRUE)

(psis_loo_MR3 <- loo(loglik_MR3))


#### compare models ####

MR_only_mod_comparison <- loo_compare(psis_loo_MR1, psis_loo_MR2, psis_loo_MR3)

# weights
loo_list_MR_only <- list(psis_loo_MR1, psis_loo_MR2, psis_loo_MR3)
MR_only_mod_wt <- loo_model_weights(loo_list_MR_only)


### check diagnostics

library(shinystan)
launch_shinystan(MR_MOD3_summary)


#### models with respiratory organ in place of thermoregulatory strategy ####

## MR2_LG ##

data_matrix_mod2 <- model.matrix( ~ LogCenteredMeanMassMR + 
                                    Scaled.Centered.MR.InverseTemp + 
                                    LungsGills, 
                                    data = Overlap.MR.RSA)

dat <- list(N=nrow(Overlap.MR.RSA),
            MR=Overlap.MR.RSA[, 'LogMeanWholeOrganismMRWatts'],
            data_matrix = data_matrix_mod2,
            K = ncol(data_matrix_mod2),
            d_mat = d_mat,
            vcov_mat = vcov_mat)

MR_MOD2_LG <- stan(file = "MR_data_matrix_mod.stan",
                   data=dat,
                   iter=10000,
                   warmup = 2000,
                   thin= 10,
                   chains=4,
                   control=list("adapt_delta"=0.99,
                                 max_treedepth = 18),
                   pars=c("beta",
                          "lambda",
                          "sigma",
                          "log_lik",
                          "sigma_total"))


## MR3_LG ##

data_matrix_mod3 <- model.matrix( ~ LogCenteredMeanMassMR * 
                                    LungsGills +
                                    Scaled.Centered.MR.InverseTemp,
                                    data = Overlap.MR.RSA)

dat <- list(N=nrow(Overlap.MR.RSA),
            MR=Overlap.MR.RSA[, 'LogMeanWholeOrganismMRWatts'],
            data_matrix = data_matrix_mod3,
            K = ncol(data_matrix_mod3),
            d_mat = d_mat,
            vcov_mat = vcov_mat)

MR_MOD3_LG <- stan(file = "MR_data_matrix_mod.stan",
                   data=dat,
                   iter=10000,
                   warmup = 2000,
                   thin= 10,
                   chains=4,
                   control=list("adapt_delta"=0.99,
                                 max_treedepth = 18),
                   pars=c("beta",
                          "lambda",
                          "sigma",
                          "log_lik",
                          "sigma_total"))


#### model comparison for models with respiratory organ ####

## MR2_LG ##

posterior_MR2 <- as.data.frame(MR_MOD2_LG)
posterior_MR2 <- Thin(posterior_MR2, By = 10)
posterior_MR2 <- as.data.frame(posterior_MR2)

posterior_MR2_short <- posterior_MR2 %>% dplyr::select("beta[1]", "beta[2]", 
                                          "beta[3]", "beta[4]") %>%
                                   rename(a = "beta[1]",
                                          bMass = "beta[2]",
                                          bTemp = "beta[3]",
                                          bTherm = "beta[4]"
                                          )

y_MR2 <- Overlap.MR.RSA$LogMeanWholeOrganismMRWatts
N_MR2 <- length(y_MR2)
S_MR2 <- nrow(posterior_MR2)
loglik_MR2 <- yloo_MR2 <- sdloo_MR2 <- matrix(nrow = S_MR2, ncol = N_MR2)

st_MR2 <- as.matrix(posterior_MR2[grep("sigma_total", colnames(posterior_MR2))])
str(st_MR2)
dim(st_MR2)

sigma_array_MR2 <- array(st_MR2, dim = c(nrow(st_MR2), 109, 109))
str(sigma_array_MR2)
dim(sigma_array_MR2)

s1_MR2 <- sigma_array_MR2[2, ,]
s2_MR2 <- st_MR2[2, ]
s1_MR2 == s2_MR2

for (s in 1:S_MR2) {
  p_MR2 <- posterior_MR2_short[s, ] 
  eta_MR2 <- p_MR2$a + p_MR2$bMass * Overlap.MR.RSA$LogCenteredMeanMassMR + 
             p_MR2$bTemp * Overlap.MR.RSA$Scaled.Centered.MR.InverseTemp + 
             p_MR2$bTherm * Overlap.MR.RSA$ThermoStrat
  Cinv_MR2 <- solve(sigma_array_MR2[s, , ]) 
  g_MR2 <- Cinv_MR2 %*% (y_MR2 - eta_MR2) 
  cbar_MR2 <- diag(Cinv_MR2) 
  yloo_MR2[s, ] <- y_MR2 - g_MR2 / cbar_MR2 
  sdloo_MR2[s, ] <- sqrt(1 / cbar_MR2) 
  loglik_MR2[s, ] <- dnorm(y_MR2, yloo_MR2[s, ], sdloo_MR2[s, ], log = TRUE)
}

log_ratios_MR2 <- -loglik_MR2 
r_eff_MR2 <- relative_eff(exp(loglik_MR2), chain_id = rep(1:4, each = 8000)) 
psis_result_MR2  <- psis(log_ratios_MR2, r_eff = r_eff_MR2)

plot(psis_result_MR2, label_points = TRUE)

(psis_loo_MR2_LG <- loo(loglik_MR2))


## MR3_LG ##

posterior_MR3 <- as.data.frame(MR_MOD3_LG)
posterior_MR3 <- Thin(posterior_MR3, By = 10)
posterior_MR3 <- as.data.frame(posterior_MR3)

posterior_MR3_short <- posterior_MR3 %>% dplyr::select(
                                "beta[1]", "beta[2]", "beta[3]", "beta[4]", 
                                "beta[5]") %>%
                                rename(a = "beta[1]",
                                          bMass = "beta[2]",
                                          bTherm = "beta[3]",
                                          bTemp = "beta[4]",
                                          bMass_Therm = "beta[5]"
                                          )

y_MR3 <- Overlap.MR.RSA$LogMeanWholeOrganismMRWatts
N_MR3 <- length(y_MR3)
S_MR3 <- nrow(posterior_MR3)
loglik_MR3 <- yloo_MR3 <- sdloo_MR3 <- matrix(nrow = S_MR3, ncol = N_MR3)

st_MR3 <- as.matrix(posterior_MR3[grep("sigma_total", colnames(posterior_MR3))])
str(st_MR3)
dim(st_MR3)

sigma_array_MR3 <- array(st_MR3, dim = c(nrow(st_MR3), 109, 109))
str(sigma_array_MR3)
dim(sigma_array_MR3)

s1_MR3 <- sigma_array_MR3[2, ,]
s2_MR3 <- st_MR3[2, ]
s1_MR3 == s2_MR3

for (s in 1:S_MR3) {
  p_MR3 <- posterior_MR3_short[s, ] 
  eta_MR3 <- p_MR3$a + p_MR3$bMass * Overlap.MR.RSA$LogCenteredMeanMassMR  + 
             p_MR3$bTherm * Overlap.MR.RSA$ThermoStrat +  
             p_MR3$bTemp * Overlap.MR.RSA$Scaled.Centered.MR.InverseTemp + 
             p_MR3$bMass_Therm * Overlap.MR.RSA$LogCenteredMeanMassMR * 
                 Overlap.MR.RSA$ThermoStrat
  Cinv_MR3 <- solve(sigma_array_MR3[s, , ]) 
  g_MR3 <- Cinv_MR3 %*% (y_MR3 - eta_MR3) 
  cbar_MR3 <- diag(Cinv_MR3) 
  yloo_MR3[s, ] <- y_MR3 - g_MR3 / cbar_MR3 
  sdloo_MR3[s, ] <- sqrt(1 / cbar_MR3) 
  loglik_MR3[s, ] <- dnorm(y_MR3, yloo_MR3[s, ], sdloo_MR3[s, ], log = TRUE)
}

log_ratios_MR3 <- -loglik_MR3
r_eff_MR3 <- relative_eff(exp(loglik_MR3), chain_id = rep(1:4, each = 8000)) 
psis_result_MR3 <- psis(log_ratios_MR3, r_eff = r_eff_MR3)

plot(psis_result_MR3, label_points = TRUE)

(psis_loo_MR3_LG <- loo(loglik_MR3))

## compare models ##

MR_only_mod_comparison <- compare(psis_loo_MR1, psis_loo_MR2, psis_loo_MR3,
                                  psis_loo_MR2_LG, psis_loo_MR3_LG)

# weights
loo_list_MR_only <- list(psis_loo_MR1, psis_loo_MR2, psis_loo_MR3,
                         psis_loo_MR2_LG, psis_loo_MR3_LG)

MR_only_mod_wt <- loo_model_weights(loo_list_MR_only)

