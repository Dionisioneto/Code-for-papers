
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
      vector[N] beta = exp(xb * bs); // Cálculo de beta
    }

    model {
      // Priori alpha
      alpha ~ normal(0,10);

      // Priori (regression coef in beta)
      for (j in 1:s) {
        bs[j] ~ normal(0, 10);
      }

      // Priori lambda
      lambda ~  gamma(0.1,0.1);

      // log-likelihood
      for (i in 1:N) {
        target += delta[i] * log(beta[i] * lambda * exp((beta[i] - beta[i] * exp(alpha * time[i])) / alpha + alpha * time[i])) -
                  delta[i] * log((1 - ((1-lambda) * exp((beta[i] - beta[i] * exp(alpha * time[i])) / alpha)))^2) +
                  (1 - delta[i]) * (log(lambda * exp(-(beta[i] / alpha) * (exp(alpha * time[i]) - 1)))) -
                  (1 - delta[i]) * (log(1 - ((1 - lambda) * exp((-beta[i] / alpha) * (exp(alpha * time[i]) - 1)))));
      }
    }
      // Metrics for comparison
      generated quantities {
        vector[N] log_lik;
        for (i in 1:N) {
         log_lik[i] = delta[i] * log(beta[i] * lambda * exp((beta[i] - beta[i] * exp(alpha * time[i])) / alpha + alpha * time[i])) -
                  delta[i] * log((1 - ((1-lambda) * exp((beta[i] - beta[i] * exp(alpha * time[i])) / alpha)))^2) +
                  (1 - delta[i]) * (log(lambda * exp(-(beta[i] / alpha) * (exp(alpha * time[i]) - 1)))) -
                  (1 - delta[i]) * (log(1 - ((1 - lambda) * exp((-beta[i] / alpha) * (exp(alpha * time[i]) - 1)))));
        }
        
}


