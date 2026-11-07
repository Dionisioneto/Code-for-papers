# --- Stan Code: Weibull Mixture Cure Model ---
weibull_app = "
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
"

# --- Stan Code: Log-Normal Mixture Cure Model ---
log_normal_app = "
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
"

writeLines(weibull_app, con = "weibull_app.stan")
writeLines(log_normal_app, con = "log_normal_app.stan")


# --- Stan Code: Promotion Weibull ---
promotion_weibull_app = "
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
"

# --- Stan Code: Promotion Log-Normal ---
promotion_lognormal_app = "
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
"

writeLines(promotion_weibull_app, con = "promotion_weibull_app.stan")
writeLines(promotion_lognormal_app, con = "promotion_lognormal_app.stan")

## semiparametric versions

cod_promotion_semiparametric = "
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
"


cod_mc_semiparametric = "
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
"

writeLines(cod_promotion_semiparametric, con = "cod_promotion_semiparametric.stan")
writeLines(cod_mc_semiparametric, con = "cod_mc_semiparametric.stan")

## Defective
cod_gompertz_app = "
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
"

writeLines(cod_gompertz_app, con = "cod_gompertz_app.stan")


