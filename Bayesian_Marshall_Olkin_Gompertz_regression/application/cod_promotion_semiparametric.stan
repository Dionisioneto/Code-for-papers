
          data {
        int<lower=1> N;                        // Number of patients
        int<lower=1> s;                        // Number of covariates (including intercept)
        matrix[N, s] X;                        // Design matrix
        
        vector<lower=0>[N] time;              // Observed times (failure or censoring)
        array[N] int<lower=0, upper=1> delta;  // Indicator (1 = failure, 0 = censoring)
        
        // Structure for the semiparametric baseline (Piecewise)
        int<lower=1> J;                        // Total number of time intervals
        array[N] int<lower=1, upper=J> t_idx;  // Index of the final interval where the observation fell
        matrix[N, J] D;                        // Exposure matrix: D[i,j] = time individual i spent in interval j
      }
      
      parameters {
        vector[s] beta;                        // Predictor effects on the cure fraction
        
        // Parameters for the non-parametric F(t) function
        vector[J] log_lambda;                  // Baseline log-hazard in each interval
        real<lower=0> sigma_lambda;            // Standard deviation for the smoothing prior (RW1)
      }
      
      transformed parameters {
        vector[N] theta;
        vector[N] F_t;
        vector[N] f_t;
        vector[J] lambda = exp(log_lambda);
        
        // Linear predictors of the cure rate
        theta = exp(X * beta); 
        
        // Exact calculation of F(t) and f(t) for each individual based on exposure
        for (i in 1:N) {
          real Lambda_t = 0;
          
          // Sums the cumulative hazard over the intervals the individual was at risk
          for (j in 1:t_idx[i]) {
            Lambda_t += lambda[j] * D[i, j];
          }
          
          // F(t) and f(t) derived from the hazard function
          F_t[i] = 1 - exp(-Lambda_t);
          f_t[i] = lambda[t_idx[i]] * exp(-Lambda_t);
        }
      }
      
      model {
        // 1. Priors for the parametric coefficients
        beta ~ normal(0, 5);
        
        // 2. Semiparametric Smoothing Prior (Random Walk of order 1)
        sigma_lambda ~ normal(0, 1); // Half-normal prior due to the lower bound <lower=0>
        log_lambda[1] ~ normal(-3, 5);
        for (j in 2:J) {
          log_lambda[j] ~ normal(log_lambda[j-1], sigma_lambda);
        }
        
        // 3. Marginal Likelihood
        for (i in 1:N) {
          if (delta[i] == 1) {
            // Failure: log(h(t)) + log(S(t)) => log(theta * f(t)) - theta * F(t)
            target += log(theta[i]) + log(f_t[i]) - theta[i] * F_t[i];
          } else {
            // Censoring: log(S(t)) => - theta * F(t)
            target += - theta[i] * F_t[i];
          }
        }
      }
      
      generated quantities {
        vector[N] log_lik; // Strict storage for WAIC
        
        for (i in 1:N) {
          if (delta[i] == 1) {
            log_lik[i] = log(theta[i]) + log(f_t[i]) - theta[i] * F_t[i];
          } else {
            log_lik[i] = - theta[i] * F_t[i];
          }
        }
      }

