library(rstan)

source("~/Desktop/PhD/papers/Cap4_MO_Gompertz/application/source_functions.R")
source("~/Desktop/PhD/papers/Cap4_MO_Gompertz/application/stan_other_cure_models.R")

# Data preparation and Model Fitting
wmu = as.matrix(test_data[, c("s2","s3")])
data_test = list(N = dim(test_data)[1], s = dim(wmu)[2], time = test_data$time, delta = test_data$delta, X = wmu)

# Fit and evaluate Weibull
weibull_fit = stan(file = 'weibull_app.stan', data = data_test, chains = 2, iter = 5000, warmup = 2500, init = 0, seed = 100)
summary(weibull_fit, pars = c("alpha", "mu_base","gamma0", "gamma"))$summary[,c("mean", "sd", "2.5%", "97.5%", "n_eff", "Rhat")]
comp_model_metrics(weibull_fit, log_lik_name = "log_lik")

# Fit and evaluate Log-Normal
lognormal_fit = stan(file = 'log_normal_app.stan', data = data_test, chains = 2, iter = 5000, warmup = 2500, seed = 100, cores = 2)
summary(lognormal_fit, pars = c("mu_ln", "sigma_ln","gamma0", "gamma"))$summary[,c("mean", "sd", "2.5%", "97.5%", "n_eff", "Rhat")]
comp_model_metrics(lognormal_fit, log_lik_name = "log_lik")


# Fit and evaluate Promotion Weibull
promotion_weibull_fit = stan(file = 'promotion_weibull_app.stan', data = data_test, chains = 2, iter = 5000, warmup = 2500, seed = 100, cores = 2)
summary(promotion_weibull_fit, pars = c("alpha", "mu_base","gamma0", "gamma"))$summary[,c("mean", "sd", "2.5%", "97.5%", "n_eff", "Rhat")]
comp_model_metrics(promotion_weibull_fit, log_lik_name = "log_lik")

# Fit and evaluate Promotion Log-Normal
promotion_lognormal_fit = stan(file = 'promotion_lognormal_app.stan', data = data_test, chains = 2, iter = 5000, warmup = 2500, seed = 100, cores = 2)
summary(promotion_lognormal_fit, pars = c("mu_ln", "sigma_ln","gamma0", "gamma"))$summary[,c("mean", "sd", "2.5%", "97.5%", "n_eff", "Rhat")]
comp_model_metrics(promotion_lognormal_fit, log_lik_name = "log_lik")


data_regmog = list(N = dim(test_data)[1], s = dim(wmu)[2]+1, time = test_data$time, delta = test_data$delta, xb = cbind(1,wmu))

# Fit and Evaluate Gompertz
gompertz_app_fit = stan(file = 'cod_gompertz_app.stan', data = data_regmog, chains = 2, iter = 5000, warmup = 2500, seed = 100, cores = 2)
traceplot(gompertz_app_fit, pars = c("alpha","bs"))
comp_model_metrics(gompertz_app_fit, log_lik_name = "log_lik")


# 2. Build exposure matrix D and interval index t_idx

data_semiparametric = prepare_semiparametric(test_data, n_intervals = 9, cov_formula = ~ (s2 + s3))

# --- Fit Promotion time Semiparametric Model ---
promotion_semiparametric_fit = stan(file = "cod_promotion_semiparametric.stan",
                                    data = data_semiparametric$stan_data,
                                    iter = 5000,
                                    warmup = 2500,
                                    chains = 2)

#traceplot(promotion_semiparametric_fit)
head(summary(promotion_semiparametric_fit)$summary, 1)
comp_model_metrics(promotion_semiparametric_fit, log_lik_name = "log_lik")


# --- Fit Mixture Cure Semiparametric Model ---
mc_semiparametric_fit = stan(file = "cod_mc_semiparametric.stan",
                                    data = data_semiparametric$stan_data,
                                    iter = 8000,
                                    warmup = 3000,
                                    chains = 2)

traceplot(mc_semiparametric_fit)
head(summary(mc_semiparametric_fit)$summary, 10)
comp_model_metrics(mc_semiparametric_fit, log_lik_name = "log_lik")
































