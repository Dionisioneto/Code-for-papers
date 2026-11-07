## Source for functions
library(loo)

comp_model_metrics = function(fit_model, log_lik_name = "log_lik") {
  
  # 1. Extract the log-likelihood matrix from the Stan fit
  log_lik_mat = tryCatch({
    loo::extract_log_lik(fit_model, parameter_name = log_lik_name, merge_chains = TRUE)
  }, error = function(e) {
    stop("Could not find '", log_lik_name, "' in the fit object. Check the 'generated quantities' block.")
  })
  
  # 2. Calculate LOO and WAIC (suppressing warnings for Pareto k > 0.7)
  suppressWarnings({
    loo_res = loo::loo(log_lik_mat)
    waic_res = loo::waic(log_lik_mat)
  })
  
  # 3. Calculate LPML (Expected log predictive density)
  lpml_val = loo_res$estimates["elpd_loo", "Estimate"]
  
  # 4. Calculate DIC (Deviance Information Criterion)
  deviance_samples = -2 * rowSums(log_lik_mat)
  D_bar = mean(deviance_samples)            # Posterior mean deviance
  pD = var(deviance_samples) / 2            # Effective number of parameters
  dic_val = D_bar + pD
  
  # 5. Extract absolute WAIC and LOOIC
  waic_val = waic_res$estimates["waic", "Estimate"]
  looic_val = loo_res$estimates["looic", "Estimate"]
  
  # Return structured results
  results = data.frame(
    Metric = c("WAIC", "PSIS-LOO", "LPML", "DIC"),
    Value = round(c(waic_val, looic_val, lpml_val, dic_val), 4)
  )
  
  return(results)
}


## for semiparametric cure models

prepare_semiparametric = function(data, n_intervals = 9, cov_formula = ~ (s2 + s3)) {
  # 1. Define cut points (knots) based on failure times
  # n_intervals requires (n_intervals + 1) cut points
  cuts = quantile(data$time[data$delta == 1], probs = seq(0, 1, length.out = n_intervals + 1))
  cuts[1] = 0
  cuts[length(cuts)] = max(data$time) * 1.01 # Guarantee inclusion of max time
  
  J = length(cuts) - 1
  N = nrow(data)
  
  # 2. Build exposure matrix D and interval index t_idx
  D = matrix(0, nrow = N, ncol = J)
  t_idx = numeric(N)
  
  for (i in 1:N) {
    t_i = data$time[i]
    
    if (t_i == 0) t_i = 1e-6 # Ensure t_i > 0 for indexing
    
    idx = sum(t_i > cuts[-(J+1)]) 
    t_idx[i] = max(1, idx)          # Enforce minimum index of 1
    
    # Fill time spent in each interval
    for(j in 1:t_idx[i]) {
      if (j < t_idx[i]) {
        D[i, j] = cuts[j+1] - cuts[j] # Spanned full interval
      } else {
        D[i, j] = t_i - cuts[j]       # Partial time in terminal interval
      }
    }
  }
  
  # Generate the model matrix once to avoid repeating code
  X_mat = model.matrix(cov_formula, data = data)
  
  # 3. Format data structures for Stan
  stan_data = list(
    N = N,
    s = ncol(X_mat), 
    X = X_mat,
    time = data$time,
    delta = data$delta,
    J = J,
    t_idx = t_idx,
    D = D
  )
  
  
  
  # Return both lists so they can be accessed easily
  return(list(
    stan_data = stan_data, 
    cuts = cuts
  ))
}






