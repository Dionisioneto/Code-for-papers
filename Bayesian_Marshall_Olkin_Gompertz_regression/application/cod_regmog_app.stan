
    data {
      int<lower=0> N; 
      int<lower=0> s; 

      vector[N] time;  
      vector[N] delta; 

      matrix[N, s] xb; 
    }

    parameters {
      real alpha;
      vector[s] bs; 
      real<lower=0> lambda;
    }

    transformed parameters {
      vector[N] beta = exp(xb * bs);
    }

    model {
      alpha ~ normal(0,10);
      for (j in 1:s) {
        bs[j] ~ normal(0, 10);
      }
      lambda ~ gamma(0.1,0.1);

      // Custom log-likelihood for the survival model
      for (i in 1:N) {
        target += delta[i] * log(beta[i] * lambda * exp((beta[i] - beta[i] * exp(alpha * time[i])) / alpha + alpha * time[i])) -
                  delta[i] * log((1 - ((1-lambda) * exp((beta[i] - beta[i] * exp(alpha * time[i])) / alpha)))^2) +
                  (1 - delta[i]) * (log(lambda * exp(-(beta[i] / alpha) * (exp(alpha * time[i]) - 1)))) -
                  (1 - delta[i]) * (log(1 - ((1 - lambda) * exp((-beta[i] / alpha) * (exp(alpha * time[i]) - 1)))));
      }
    }
    
    generated quantities {
      vector[N] log_lik;
      // Log-likelihood for model comparison
      for (i in 1:N) {
         log_lik[i] = delta[i] * log(beta[i] * lambda * exp((beta[i] - beta[i] * exp(alpha * time[i])) / alpha + alpha * time[i])) -
                  delta[i] * log((1 - ((1-lambda) * exp((beta[i] - beta[i] * exp(alpha * time[i])) / alpha)))^2) +
                  (1 - delta[i]) * (log(lambda * exp(-(beta[i] / alpha) * (exp(alpha * time[i]) - 1)))) -
                  (1 - delta[i]) * (log(1 - ((1 - lambda) * exp((-beta[i] / alpha) * (exp(alpha * time[i]) - 1)))));
      }

      // Cure fraction per scenario (Stage I, II, III)
      real p01 = (lambda * exp(exp(bs[1]) / alpha)) / (1 - (1 - lambda) * exp(exp(bs[1]) / alpha));
      real p02 = (lambda * exp(exp(bs[1] + bs[2]) / alpha)) / (1 - (1 - lambda) * exp(exp(bs[1] + bs[2]) / alpha));
      real p03 = (lambda * exp(exp(bs[1] + bs[3]) / alpha)) / (1 - (1 - lambda) * exp(exp(bs[1] + bs[3]) / alpha));
    }

