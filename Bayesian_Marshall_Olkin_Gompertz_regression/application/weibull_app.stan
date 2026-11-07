
    data {
      int<lower=1> N;
      vector[N] time;
      array[N] int<lower=0, upper=1> delta;
      int<lower=1> s;
      matrix[N, s] X;
    }
    
    parameters {
      real gamma0;                       // Intercept: cure probability
      vector[s] gamma;                   // Covariates for cure probability
      real<lower=0> alpha;               // Weibull shape
      real<lower=0> mu_base;             // Weibull scale parameter
    }
    
    model {
      vector[N] logit_p_cure = gamma0 + X * gamma;
    
      gamma0 ~ normal(0, 10);
      gamma  ~ normal(0, 10);
      alpha  ~ lognormal(0, 10);
      mu_base ~ lognormal(0, 10);
    
      for (i in 1:N) {
        real log_p_cure = log_inv_logit(logit_p_cure[i]);
        real log_p_uncured = log1m_inv_logit(logit_p_cure[i]);
        real log_S_u = -mu_base * pow(time[i], alpha);
        real log_f_u = log(mu_base) + log(alpha) + (alpha - 1) * log(time[i]) + log_S_u;
    
        if (delta[i] == 1) {
          target += log_p_uncured + log_f_u; // Event
        } else {
          target += log_sum_exp(log_p_cure, log_p_uncured + log_S_u); // Censored
        }
      }
    }
    
    generated quantities {
      vector[N] log_lik;
      for (i in 1:N) {
        real linear_pred = gamma0 + dot_product(X[i], gamma);
        real log_p_cure = log_inv_logit(linear_pred);
        real log_p_uncured = log1m_inv_logit(linear_pred);
        real log_S_u = -mu_base * pow(time[i], alpha);
        real log_f_u = log(mu_base) + log(alpha) + (alpha - 1) * log(time[i]) + log_S_u;
    
        if (delta[i] == 1) {
          log_lik[i] = log_p_uncured + log_f_u;
        } else {
          log_lik[i] = log_sum_exp(log_p_cure, log_p_uncured + log_S_u);
        }
      }
    }

