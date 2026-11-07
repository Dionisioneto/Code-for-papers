
    data {
      int<lower=1> N;
      vector<lower=1e-10>[N] time;
      array[N] int<lower=0, upper=1> delta;
      int<lower=1> s;
      matrix[N, s] X;
    }
    
    transformed data {
      vector[N] log_time = log(time); // Optimization
    }
    
    parameters {
      real gamma0;                    // Intercept for expected number of promoters
      vector[s] gamma;                // Covariates
      real<lower=0> alpha;            // Weibull shape
      real<lower=0> mu_base;          // Weibull scale
    }
    
    model {
      vector[N] log_lambda = gamma0 + X * gamma;
    
      gamma0 ~ normal(0, 10);
      gamma  ~ normal(0, 10);
      alpha   ~ lognormal(0, 5);
      mu_base ~ lognormal(0, 5);
    
      for (i in 1:N) {
        real lambda_i = exp(log_lambda[i]);
        real F0 = 1 - exp(-mu_base * pow(time[i], alpha));
        real log_f0 = log(mu_base) + log(alpha) + (alpha - 1) * log_time[i] - mu_base * pow(time[i], alpha);
    
        target += delta[i] * (log_lambda[i] + log_f0) - lambda_i * F0;
      }
    }
    
    generated quantities {
      vector[N] log_lik;
      vector[N] p_cure;
      for (i in 1:N) {
        real log_lambda = gamma0 + dot_product(X[i], gamma);
        real lambda_i = exp(log_lambda);
        real F0 = 1 - exp(-mu_base * pow(time[i], alpha));
        real log_f0 = log(mu_base) + log(alpha) + (alpha - 1) * log_time[i] - mu_base * pow(time[i], alpha);
    
        p_cure[i] = exp(-lambda_i); // Probability of zero promoters
        log_lik[i] = delta[i] * (log_lambda + log_f0) - lambda_i * F0;
      }
    }

