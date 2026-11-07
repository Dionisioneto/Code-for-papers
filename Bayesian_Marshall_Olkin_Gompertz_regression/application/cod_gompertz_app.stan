
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
      }
      
      transformed parameters {
        vector[N] beta = exp(xb * bs);
      }
      
      model {
        alpha ~ normal(0, 10);
        bs ~ normal(0, 10); // Vectorized prior
      
        for (i in 1:N) {
          real log_S = -(beta[i] / alpha) * (exp(alpha * time[i]) - 1);
          if (delta[i] == 1) {
            target += log(beta[i]) + alpha * time[i] + log_S;
          } else {
            target += log_S;
          }
        }
      }
      
      generated quantities {
        vector[N] log_lik;
        real p01;  // Stage I cure fraction limit
        real p02;  // Stage II cure fraction limit
        real p03;  // Stage III cure fraction limit
        
        for (i in 1:N) {
          real log_S = -(beta[i] / alpha) * (exp(alpha * time[i]) - 1);
          if (delta[i] == 1) {
            log_lik[i] = log(beta[i]) + alpha * time[i] + log_S;
          } else {
            log_lik[i] = log_S;
          }
        }
        // Derived limits representing Gompertz cure fractions
        p01 = exp(exp(bs[1]) / alpha);
        p02 = exp(exp(bs[1] + bs[2]) / alpha);
        p03 = exp(exp(bs[1] + bs[3]) / alpha);
      }

