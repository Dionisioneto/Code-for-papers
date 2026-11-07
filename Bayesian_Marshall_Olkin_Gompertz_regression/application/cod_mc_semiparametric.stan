
    data {
          int<lower=1> N;                        // Number of patients
          int<lower=1> s;                    // Number of covariates for the cure probability (logit)

          matrix[N, s] X;                    // Design matrix for the cure fraction (logistic)
          
          vector<lower=0>[N] time;              // Observed times
          vector[N] delta;                       // Indicator (1 = failure, 0 = censoring)
          
          // Structure for the semiparametric baseline (Piecewise) of the uncured
          int<lower=1> J;                        // Number of intervals
          array[N] int<lower=1, upper=J> t_idx;  // Index of the final interval
          matrix[N, J] D;                        // Exposure matrix
        }
        
        parameters {
          vector[s] gamma;                   // Logistic coefficients for the cure probability
          
          // Parameters for the semiparametric survival of the uncured (S_u)
          vector[J] log_lambda;                  // Baseline log-hazard of the uncured per interval
          real<lower=0> sigma_lambda;            // Standard deviation for smoothing (RW)
        }
        
        transformed parameters {
          vector[N] p_cure;                      // Cure probability
          vector[N] S_u;                         // Survival of the uncured
          vector[N] f_u;                         // Density of the uncured
          vector[J] lambda = exp(log_lambda);
          
          // 1. Incidence Submodel (Cure Fraction via Logistic)
          for (i in 1:N) {
            p_cure[i] = inv_logit(X[i] * gamma);
          }
          
          // 2. Latency Submodel (Uncured with Piecewise Constant)
          for (i in 1:N) {
            real H_u = 0; // Cumulative hazard of the uncured
            for (j in 1:t_idx[i]) {
              H_u += lambda[j] * D[i, j];
            }
            S_u[i] = exp(-H_u);
            f_u[i] = lambda[t_idx[i]] * S_u[i];
          }
        }
        
        model {
          // 1. Priors
          gamma ~ normal(0, 10);                  // Logistic coefficients
          sigma_lambda ~ normal(0, 5);           // Baseline hazard smoothing
          log_lambda[1] ~ normal(0, 5);
          
          // Random Walk (RW1) for the hazard of the uncured
          for (j in 2:J) {
            log_lambda[j] ~ normal(log_lambda[j-1], sigma_lambda);
          }
          
          // 2. Marginal Likelihood of the Mixture Model
          // For failure (delta = 1): f_pop(t) = (1 - p_cure) * f_u(t)
          // For censoring (delta = 0): S_pop(t) = p_cure + (1 - p_cure) * S_u(t)
          for (i in 1:N) {
            if (delta[i] == 1) {
              target += log(1 - p_cure[i]) + log(f_u[i]);
            } else {
              target += log(p_cure[i] + (1 - p_cure[i]) * S_u[i]);
            }
          }
        }
        
        generated quantities {
          vector[N] log_lik;
          
          for (i in 1:N) {
            if (delta[i] == 1) {
              log_lik[i] = log(1 - p_cure[i]) + log(f_u[i]);
            } else {
              log_lik[i] = log(p_cure[i] + (1 - p_cure[i]) * S_u[i]);
            }
          }
        }

