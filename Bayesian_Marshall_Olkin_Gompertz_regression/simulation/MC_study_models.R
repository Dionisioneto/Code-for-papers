### -----------
### Monte Carlo Study 
### Model comparison
### Generate data and estimate via STAN
### -----------

library(rstan)
library(loo)

metrics_stan = function(fit_model, log_lik_name = "log_lik") {
  
  # 1. Extract the log-likelihood matrix
  log_lik_mat = tryCatch({
    loo::extract_log_lik(fit_model, parameter_name = log_lik_name, merge_chains = TRUE)
  }, error = function(e) {
    stop("It was not possible to find '", log_lik_name, "' in the fitted model. Verify 'generated quantities' in the Stan code.")
  })
  
  # 2. Compute LOO and WAIC using the loo package
  suppressWarnings({
    loo_res = loo::loo(log_lik_mat)
    waic_res = loo::waic(log_lik_mat)
  })
  
  # 3. Calculate the LPML (Log Pseudo Marginal Likelihood)
  # Corresponds to the elpd_loo (Expected Log Predictive Density) estimate
  lpml_val = loo_res$estimates["elpd_loo", "Estimate"]
  
  # 4. Calculate the DIC (Deviance Information Criterion)
  # Deviance at each MCMC iteration: D(theta) = -2 * log(p(y | theta))
  deviance_samples = -2 * rowSums(log_lik_mat)
  D_bar = mean(deviance_samples)           # Posterior mean of the deviance
  pD = var(deviance_samples) / 2           # Effective number of parameters (approximation)
  dic_val = D_bar + pD
  
  # 5. Extract absolute WAIC and LOOIC values
  waic_val = waic_res$estimates["waic", "Estimate"]
  looic_val = loo_res$estimates["looic", "Estimate"]
  
  # Structure the final result in a data frame to facilitate comparisons
  results = data.frame(
    Metrica = c("WAIC", "PSIS-LOO", "-2*LPML", "DIC"),
    Valor = round(c(waic_val, looic_val, -2*lpml_val, dic_val), 4)
  )
  
  return(results)
}


# Base functions

St_Gompertz = function(t,alpha,beta){
  st = exp(-(beta/alpha)*(exp(alpha*t)-1))
  return(st)
}

Ft_Gompertz = function(t,alpha,beta){
  Ft = 1 - exp(-(beta/alpha)*(exp(alpha*t)-1))
  return(Ft)
}

ft_Gompertz = function(t,alpha,beta){
  ft = beta*exp(alpha*t)*exp(-(beta/alpha)*(exp(alpha*t)-1))
  return(ft)
}

# MO Gompertz

ftmo_gompertz = function(t,alpha,beta,lambda){
  num = lambda*ft_Gompertz(t=t,alpha=alpha,beta=beta)
  dem = (1-(1-lambda)*St_Gompertz(t=t,alpha=alpha,beta=beta))^2
  return(num/dem)
}

Ftmo_gompertz = function(t,alpha,beta,lambda){
  num = lambda*St_Gompertz(t=t,alpha=alpha,beta=beta)
  dem = 1 - ((1-lambda)*St_Gompertz(t=t,alpha=alpha,beta=beta))
  return(1-(num/dem))
}

Stmo_gompertz = function(t,alpha,beta,lambda){
  num = lambda*St_Gompertz(t=t,alpha=alpha,beta=beta)
  dem = 1 - ((1-lambda)*St_Gompertz(t=t,alpha=alpha,beta=beta))
  return(num/dem)
}

## Generate data

## inserir xa e xb sem o vetor de 1's do intercepto
gen.cure.reg.mog = function(n,a,bs,xb,l){#, max_cen=20){
  finv.mog = function(t,alpha,beta,lambda,unif){Ftmo_gompertz(t=t,alpha=alpha,beta=beta,lambda=lambda) - unif}
  ger.mog = function(alpha,beta,lambda,unif){t=uniroot(finv.mog,c(0,1000),tol=0.01,alpha=alpha,beta=beta,lambda=lambda,unif=unif);return(t$root)}
  
  ## individual parameters
  
  
  b = exp(cbind(1,xb) %*% bs)
  p_base = exp(b/a)
  p_cure = (l*p_base)/(1-((1-l)*p_base))
  
  t = rep(NA,n)
  
  rmi = rbinom(n=n,size=1,prob=1-p_cure)
  uni = runif(n=n,0,max=1-p_cure)
  
  for(i in 1:n){
    t[i] = ifelse(rmi[i]==0, Inf, ger.mog(alpha=a,beta=b[i,],lambda=l,unif=uni[i]))
  }
  
  t_finite = ifelse(t==Inf,0,t)
  u2 = runif(n=n,0,max(t_finite))
  #u2 = runif(n=n,0,max_cen)
  t2 = pmin(t,u2); delta = ifelse(t<u2,1,0)
  
  df = data.frame(time=t2,delta=delta,xb=xb)
  return(df)
}


cod_regmog = "
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

"

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
      // lambda foi removido
    }

    transformed parameters {
      vector[N] beta = exp(xb * bs); // Cálculo de beta
    }

    model {
      // Priori alpha
      alpha ~ normal(0, 10);

      // Priori (regression coef in beta) de forma vetorizada
      bs ~ normal(0, 10);

      // log-likelihood
      for (i in 1:N) {
        // Sobrevivência base calculada uma única vez
        real log_S = -(beta[i] / alpha) * (exp(alpha * time[i]) - 1);
        
        if (delta[i] == 1) {
          // Evento
          target += log(beta[i]) + alpha * time[i] + log_S;
        } else {
          // Censura
          target += log_S;
        }
      }
    }

    generated quantities {
      vector[N] log_lik;
      
      // Cure fraction per scenario
      real p01;  // Stage 1
      real p02;  // Stage 2
      real p03;  // Stage 3
      
      for (i in 1:N) {
        real log_S = -(beta[i] / alpha) * (exp(alpha * time[i]) - 1);
        
        if (delta[i] == 1) {
          log_lik[i] = log(beta[i]) + alpha * time[i] + log_S;
        } else {
          log_lik[i] = log_S;
        }
      }

      // Reduzido da fórmula MOG para Gompertz com lambda = 1
      p01 = exp(exp(bs[1]) / alpha);
      
      p02 = exp(exp(bs[1] + bs[2]) / alpha);

      p03 = exp(exp(bs[1] + bs[3]) / alpha);
    }
"


## To stan file 
writeLines(cod_regmog, con = "cod_regmog.stan")
writeLines(cod_gompertz_app, con = "cod_gompertz_app.stan")

## test
n0 = c(50, 1000, 300, 500, 800, 1000, 2000)

## scenario 1 (real values for parameters, average cure and average delta)
a0=-0.517; l0 = 0.485; bs0 = c(-2, 0.701, 2.487) 

params = c(a0,l0,bs0)

## scenario 2 (high cure and delta)
#a0=-0.2; l0 = 0.5; bs0 = c(-2, 0.701, 2.487) 

## scenario 3 (low cure and delta)
#a0=-2; l0 = 1.5; bs0 = c(-2, 0.701, 2.487)  

# Supondo um dataframe 'dados' com 'tempo', 'status', e covariável 'trt'
# 1. Definir os pontos de corte (knots)




cod_mc_semiparametric = "
      data {
      int<lower=1> N;                        // Número de pacientes
      int<lower=1> s;                    // Número de covariáveis para a probabilidade de cura (logit)
      int<lower=1> lat_cov;                  // Número de covariáveis para o grupo não curado (opcional, aqui omitido por simplicidade ou pode usar X)
      
      matrix[N, s] Z;                    // Matriz de desenho para a fração de cura (logística)
      
      vector<lower=0>[N] t_obs;              // Tempos observados
      vector[N] delta;                       // Indicador (1 = falha, 0 = censura)
      
      // Estrutura para o baseline semiparamétrico (Piecewise) dos não-curados
      int<lower=1> J;                        // Número de intervalos
      array[N] int<lower=1, upper=J> t_idx;  // Índice do intervalo final
      matrix[N, J] D;                        // Matriz de exposição
    }
    
    parameters {
      vector[s] gamma;                   // Coeficientes logísticos para a probabilidade de cura
      
      // Parâmetros para a sobrevivência semiparamétrica dos não-curados (S_u)
      vector[J] log_lambda;                  // Log-risco basal dos não-curados por intervalo
      real<lower=0> sigma_lambda;            // Desvio padrão da suavização (RW)
    }
    
    transformed parameters {
      vector[N] p_cure;                      // Probabilidade de cura
      vector[N] S_u;                         // Sobrevivência dos não-curados
      vector[N] f_u;                         // Densidade dos não-curados
      vector[J] lambda = exp(log_lambda);
      
      // 1. Submodelo de Incidência (Fração de Cura via Logística)
      for (i in 1:N) {
        p_cure[i] = inv_logit(Z[i] * gamma);
      }
      
      // 2. Submodelo de Latência (Não-curados com Piecewise Constant)
      for (i in 1:N) {
        real H_u = 0; // Risco cumulativo dos não-curados
        for (j in 1:t_idx[i]) {
          H_u += lambda[j] * D[i, j];
        }
        S_u[i] = exp(-H_u);
        f_u[i] = lambda[t_idx[i]] * S_u[i];
      }
    }
    
    model {
      // 1. Prioris
      gamma ~ normal(0, 10);                  // Coeficientes da logística
      sigma_lambda ~ normal(0, 10);           // Suavização do risco basal
      log_lambda[1] ~ normal(0, 10);
      
      // Random Walk (RW1) para o risco dos não-curados
      for (j in 2:J) {
        log_lambda[j] ~ normal(log_lambda[j-1], sigma_lambda);
      }
      
      // 2. Verossimilhança Marginal do Modelo de Mistura
      // Para falha (delta = 1): f_pop(t) = (1 - p_cure) * f_u(t)
      // Para censura (delta = 0): S_pop(t) = p_cure + (1 - p_cure) * S_u(t)
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

writeLines(cod_mc_semiparametric, con = "cod_mc_semiparametric.stan")

get_D_tindx = function(time, delta, breaks){
  cutss = quantile(time[delta == 1], probs = seq(0, 1, length.out = breaks))
  cutss[1] = 0
  cutss[length(cutss)] = max(time) * 1.01 
  
  J = length(cutss) - 1
  N = length(time)
  
  D = matrix(0, nrow = N, ncol = J)
  t_idx = numeric(N)
  
  for (i in 1:N) {
    t_i = time[i]
    
    if (t_i == 0) t_i = 1e-6 
    
    idx = sum(t_i > cutss[-(J+1)]) 
    t_idx[i] = max(1, idx)         
    
    for(j in 1:t_idx[i]) {
      if (j < t_idx[i]) {
        D[i, j] = cutss[j+1] - cutss[j]  
      } else {
        D[i, j] = t_i - cutss[j]         
      }
    }
  }
  
  # CORREÇÃO 2: Retornar como uma lista
  return(list(
    D = D, 
    t_idx = t_idx,
    J = J,
    cutss = cutss
  ))
}


# values=get_D_tindx(time=dataset$time, delta=dataset$delta, breaks=10)

# 3. Lista pronta para o Stan


# Saving information

model_names = c("MOG", "DG", "MC_SP_5", "MC_SP_10")
niter=1000


sim_WAIC = array(data=NA,dim=c(niter,4,length(n0)), dimnames = list(NULL,model_names,NULL))
sim_LPML = array(data=NA,dim=c(niter,4,length(n0)), dimnames = list(NULL,model_names,NULL))
censprop = array(data=NA,dim=c(niter,1,length(n0)))


comp_mogreg = rstan::stan_model(file = "cod_regmog.stan")
comp_gompertz_app = rstan::stan_model(file = "cod_gompertz_app.stan")
comp_mc_semiparametric = rstan::stan_model(file = "cod_mc_semiparametric.stan")


for(q in 1:length(n0)){
     for(k in 1:niter){
       print(paste("Iteration ", k, "for sample size: ", n0[q]))

        stages = sample(x = c(1,2,3), size = n0[q], replace = T, prob = c(1/3, 1/3, 1/3))
        W_mat = cbind(stage2 = as.numeric(stages == 2), stage3 = as.numeric(stages == 3))

        dataset = gen.cure.reg.mog(n=n0[q],a=a0, bs = bs0, xb = W_mat, l = l0)

        data_cr_mogreg = list(N = n0[q],
                               s = dim(W_mat)[2]+1,
                               time = dataset$time,
                               delta = dataset$delta,
                               xb = cbind(1,W_mat))
  

        mod_mogreg = rstan::sampling(object = comp_mogreg, 
                                      data   = data_cr_mogreg,
                                      warmup = 2500,                  
                                      iter   = 5000,
                                      chains = 2,
                                      cores  = getOption("mc.cores", 2),
                                      refresh = 1000) 
        
        mod_gompertz = rstan::sampling(object = comp_gompertz_app, 
                                        data   = data_cr_mogreg,
                                        warmup = 2500,                  
                                        iter   = 5000,
                                        chains = 2,
                                        cores  = getOption("mc.cores", 2),
                                        refresh = 1000) 
        
        #3 getting values for the semiparamtric models
        
        values_5=get_D_tindx(time=dataset$time, delta=dataset$delta, breaks=5)
        values_10=get_D_tindx(time=dataset$time, delta=dataset$delta, breaks=10)
        
        semiparametric_data_5 = list( 
          N = nrow(dataset),
          s = ncol(model.matrix(~ (xb.stage2+xb.stage3), data = dataset)), 
          lat_cov = 1, 
          Z = model.matrix(~ (xb.stage2+xb.stage3), data = dataset), 
          
          t_obs = dataset$time,
          delta = dataset$delta,
          J = values_5$J,
          t_idx = values_5$t_idx,
          D = values_5$D
        )
        
        semiparametric_data_10 = list( 
          N = nrow(dataset),
          s = ncol(model.matrix(~ (xb.stage2+xb.stage3), data = dataset)), 
          lat_cov = 1, 
          Z = model.matrix(~ (xb.stage2+xb.stage3), data = dataset), 
          
          t_obs = dataset$time,
          delta = dataset$delta,
          J = values_10$J,
          t_idx = values_10$t_idx,
          D = values_10$D
        )
        
        
        mod_semiparametric_5 =  rstan::sampling(object = comp_mc_semiparametric,
                                  data = semiparametric_data_5,
                                  warmup = 2500,                  
                                  iter   = 5000,
                                  chains = 2,
                                  cores  = getOption("mc.cores", 2),
                                  refresh = 0
                                )
        
        mod_semiparametric_10 =  rstan::sampling(object = comp_mc_semiparametric,
                                     data = semiparametric_data_10,
                                     warmup = 2500,                  
                                     iter   = 5000,
                                     chains = 2,
                                     cores  = getOption("mc.cores", 2),
                                     refresh = 0
        )
        
        
      # metrics_stan(comp_mogreg)$Valor[1]
      # metrics_stan(comp_gompertz)
      # metrics_stan(comp_semiparametric_5)
      # metrics_stan(comp_semiparametric_10)
      
      sim_WAIC[k,1:4,q]=c(metrics_stan(mod_mogreg)$Valor[1],
                        metrics_stan(mod_gompertz)$Valor[1],
                        metrics_stan(mod_semiparametric_5)$Valor[1],
                        metrics_stan(mod_semiparametric_10)$Valor[1])
      
      sim_LPML[k,1:4,q]=c(metrics_stan(mod_mogreg)$Valor[3],
                        metrics_stan(mod_gompertz)$Valor[3],
                        metrics_stan(mod_semiparametric_5)$Valor[3],
                        metrics_stan(mod_semiparametric_10)$Valor[3])
      
      censprop[k,1,q] = 1-mean(dataset$delta)
      
     }
    write.csv2(sim_WAIC[,,q], file = paste("1sim_WAIC_", n0[q], ".csv", sep = ""))
    write.csv2(sim_LPML[,,q], file = paste("1sim_LPML_", n0[q], ".csv", sep = ""))
    write.csv2(censprop[,1,q], file = paste("11censprop_", n0[q], ".csv", sep = ""))
  }





# 
# k=1; j = 5
# 
# params = c(a0,l0,bs0)
# 
# print(paste("Iteration ", k, "for sample size: ", n0[j]))
# 
# stages = sample(x = c(1,2,3), size = n0[j], replace = T, prob = c(1/3, 1/3, 1/3))
# W_mat = cbind(stage2 = as.numeric(stages == 2), stage3 = as.numeric(stages == 3))
# 
# dataset = gen.cure.reg.mog(n=n0[j],a=a0, bs = bs0, xb = W_mat, l = l0)
# 
# data_cr_mogreg = list(N = n0[j],
#                       s = dim(W_mat)[2]+1,
#                       time = dataset$time,
#                       delta = dataset$delta,
#                       xb = cbind(1,W_mat))
# 
# comp_mogreg = rstan::stan_model(file = "cod_regmog.stan")
# comp_gompertz_app = rstan::stan_model(file = "cod_gompertz_app.stan")
# 
# comp_mogreg = rstan::sampling(object = comp_mogreg, 
#                               data   = data_cr_mogreg,
#                               warmup = 2500,                  
#                               iter   = 5000,
#                               chains = 2,
#                               cores  = getOption("mc.cores", 2),
#                               refresh = 0) 
# 
# comp_gompertz = rstan::sampling(object = comp_gompertz_app, 
#                                 data   = data_cr_mogreg,
#                                 warmup = 2500,                  
#                                 iter   = 5000,
#                                 chains = 2,
#                                 cores  = getOption("mc.cores", 2),
#                                 refresh = 0) 
# 
# 












