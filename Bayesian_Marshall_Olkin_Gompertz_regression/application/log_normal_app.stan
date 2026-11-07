
    data {
      int<lower=1> N;
      vector<lower=1e-10>[N] time;
      array[N] int<lower=0, upper=1> delta;
      int<lower=1> s;
      matrix[N, s] X;
    }
    
    parameters {
      real gamma0;                    // Intercept: cure probability
      vector[s] gamma;                // Covariates for cure probability
      real mu_ln;                     // Mean of log(T) for susceptible
      real<lower=0> sigma_ln;         // Standard deviation of log(T)
    }
    
    model {
      vector[N] logit_p_cure = gamma0 + X * gamma;
    
      gamma0 ~ normal(0, 10);
      gamma  ~ normal(0, 10);
      mu_ln  ~ normal(0, 10);
      sigma_ln ~ lognormal(0, 10);
    
      for (i in 1:N) {
        real p_cure = inv_logit(logit_p_cure[i]);
        real p_uncured = 1 - p_cure;
        real log_S_u = lognormal_lccdf(time[i] | mu_ln, sigma_ln);
        real log_f_u = lognormal_lpdf(time[i] | mu_ln, sigma_ln);
    
        if (delta[i] == 1) {
          target += log(p_uncured) + log_f_u;
        } else {
          target += log_sum_exp(log(p_cure), log(p_uncured) + log_S_u);
        }
      }
    }
    
    generated quantities {
      vector[N] log_lik;
      for (i in 1:N) {
        real p_cure = inv_logit(gamma0 + dot_product(X[i], gamma));
        real p_uncured = 1 - p_cure;
        real log_S_u = lognormal_lccdf(time[i] | mu_ln, sigma_ln);
        real log_f_u = lognormal_lpdf(time[i] | mu_ln, sigma_ln);
    
        if (delta[i] == 1) {
          log_lik[i] = log(p_uncured) + log_f_u;
        } else {
          log_lik[i] = log_sum_exp(log(p_cure), log(p_uncured) + log_S_u);
        }
      }
    }

