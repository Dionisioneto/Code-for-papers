
    data {
      int<lower=1> N;
      vector<lower=1e-10>[N] time;
      array[N] int<lower=0, upper=1> delta;
      int<lower=1> s;
      matrix[N, s] X;
    }
    
    parameters {
      real gamma0;
      vector[s] gamma;
      real mu_ln;
      real<lower=0> sigma_ln;
    }
    
    model {
      vector[N] log_lambda = gamma0 + X * gamma;
      gamma0 ~ normal(0, 10);
      gamma  ~ normal(0, 10);
      mu_ln  ~ normal(0, 10);
      sigma_ln ~ lognormal(0, 10);
    
      for (i in 1:N) {
        real lambda_i = exp(log_lambda[i]);
        real F0 = exp(lognormal_lcdf(time[i] | mu_ln, sigma_ln));
        real log_f0 = lognormal_lpdf(time[i] | mu_ln, sigma_ln);
        target += delta[i] * (log_lambda[i] + log_f0) - lambda_i * F0;
      }
    }
    
    generated quantities {
      vector[N] log_lik;
      vector[N] p_cure;
      for (i in 1:N) {
        real log_lambda = gamma0 + dot_product(X[i], gamma);
        real lambda_i = exp(log_lambda);
        real F0 = exp(lognormal_lcdf(time[i] | mu_ln, sigma_ln));
        real log_f0 = lognormal_lpdf(time[i] | mu_ln, sigma_ln);
        
        p_cure[i] = exp(-lambda_i);
        log_lik[i] = delta[i] * (log_lambda + log_f0) - lambda_i * F0;
      }
    }

