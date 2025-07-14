### -----------
### Monte Carlo Kumaraswamy Gompertz defective distribution
### Data Generation and Bayesian estimation
### -----------

library(rstan)

## Reparametrized version of kumarawamy gompertz defective distribution
## k=1 (identifiable model)
## regresion coeffficents in the cuper parameter

ftkg = function(t,alpha,p,psi){
  kappa=1
  beta=alpha*log(1-(1-p^(1/kappa))^(1/psi))
  ft = kappa*psi*beta*exp(alpha*t)*exp(-beta/alpha*(exp(alpha*t)-1))*(1-exp(-beta/alpha*(exp(alpha*t)-1)))^(psi-1)*(1-(1-exp((-beta/alpha)*(exp(alpha*t)-1)))^psi)^(kappa-1)
  return(ft)
}

Stkg = function(t,alpha,p,psi){
  kappa=1
  beta=alpha*log(1-(1-p^(1/kappa))^(1/psi))
  st = (1-(1-exp((-beta/alpha)*(exp(alpha*t)-1)))^psi)^kappa
  return(st)
}

## generate data
gen.cure.reg.kg = function(n,alpha,bs,xs,psi){
  finv.kg = function(t,alpha,p,psi,unif){(1-Stkg(t=t,alpha=alpha,p=p,psi=psi)) - unif}
  ger.kg = function(alpha,p,psi,unif){t=uniroot(finv.kg,c(0,1000),tol=0.001,alpha=alpha,p=p,psi=psi,unif=unif);return(t$root)}
  
  t = rep(NA,n)
  ps = 1/(1+exp(-as.matrix(cbind(1,xs)%*%bs)))
  
  for(i in 1:n){
    rmi = rbinom(n=1,size=1,prob=1-ps[i])
    uni = runif(n=1,0,max=1-ps[i])
    t[i] = ifelse(rmi==0, Inf, ger.kg(alpha=alpha,p=ps[i],psi=psi,unif=uni))
  }
  
  t_finite = ifelse(t==Inf,0,t)
  u2 = runif(n=n,0,max(t_finite))
  t2 = pmin(t,u2); delta = ifelse(t<u2,1,0)
  
  return(as.data.frame(cbind(t2,delta,xs)))
}

cod_regkg = "
data {
  int<lower=0> N; // Número de observações
  int<lower=0> q; // Número de covariáveis para a cura

  vector[N] time;  // Tempo de falha
  vector[N] delta; // Indicador de censura
  
  matrix[N, q] X; // Matriz de covariáveis para cura (N x q)
}

parameters {
  vector[q] bs; // Coeficientes para p
  real alpha; 
  real<lower=0> psi;
}

transformed parameters {
  vector[N] p = inv_logit(X * bs);
}

model {
  // Priori para os coeficientes em p
  for (j in 1:q) {
    bs[j] ~ normal(0, 10);
  }

  // Priori para os outros parâmetros
  alpha ~ normal(-1, 10);
  psi ~ gamma(0.01, 0.01);

  // Log-verossimilhança (sem alterar sua forma original)
  for (i in 1:N) {
    real beta_i = alpha * log(1 - pow(1 - pow(p[i], 1.0), 1.0 / psi));

    target += delta[i] * (
      log(1 * psi * beta_i * exp(alpha * time[i]) *
      exp((beta_i - beta_i * exp(alpha * time[i])) / alpha) *
      pow(1 - exp((beta_i - beta_i * exp(alpha * time[i])) / alpha), psi - 1) *
      pow(1 - pow(1 - exp((beta_i - beta_i * exp(alpha * time[i])) / alpha), psi), 1 - 1))
    ) +
    (1 - delta[i]) * (
      log(pow(1 - pow(1 - exp((beta_i - beta_i * exp(alpha * time[i])) / alpha), psi), 1))
    );
  }
}
"

## Transcrever o código escrito para um file stan 
writeLines(cod_regkg, con = "cod_regkg.stan")

## 4. Tables to save data ====
n_replicas = 1000                      ## number of Monte Carlo replicates  
ncols_mc = 4                             ## number of information in each replicate (mean,sd,bias,cp)
n_amostral = c(100,300,500,1000)           ## sample sizes to study

## Tables

## bs = (b0, b1, b2, alpha, psi)
kg_mc_bs = array(data=0, dim = c(n_replicas, ncols_mc*5, length(n_amostral)),
                      dimnames = list(1:n_replicas,rep(c("mean", "sd", "bias", "cp"),5),n_amostral))


## Censoring
kg_mc_cens = array(data=0, dim = c(n_replicas, 1, length(n_amostral)))


## ---
## Scenarios: 
## ---

# 1 (45% censoring)
# a0=-2
# bs0 = c(-1,0.5,0.5)
# psi0=2
# params=c(bs0,a0,psi0)

# 2 (% censoring)
a0=-0.12
bs0 = c(-1.78,-1.49, -1.75)
psi0=0.69
params=c(bs0,a0,psi0)


for(j in 1:length(n_amostral)) {
  for(i in 1:n_replicas) {
    
    repeat {
      erro_ao_rodar <- FALSE
      
      tryCatch({
        print(paste("Realização Monte Carlo: ", i, ", tamanho amostral:", n_amostral[j]))
        
        ## Gerar dados
        x1 = rnorm(n_amostral[j])
        x2 = rbinom(n = n_amostral[j], size = 1, prob = 0.5)
        xs0 = cbind(x1, x2)
        dados.reg.kg = gen.cure.reg.kg(n = n_amostral[j], alpha = a0, bs = bs0, xs = xs0, psi = psi0)
        
        ## Dados e modelo Stan
        data_regkg = list(
          N = dim(dados.reg.kg)[1],
          q = dim(cbind(1, dados.reg.kg[,3:4]))[2],
          time = dados.reg.kg[,1],
          delta = dados.reg.kg[,2],
          X = cbind(1, dados.reg.kg[,3:4])
        )
        
        ## Rodar modelo
        regkgfit = stan(file = 'cod_regkg.stan', data = data_regkg, 
                        chains = 1, iter = 5000, warmup = 1000)
        
        ## Extração dos resultados
        post_summary = summary(regkgfit)$summary
        post_mean = post_summary[1:5,1]
        post_sd = post_summary[1:5,3]
        rbias = (post_mean - params) / params * 100
        cpost95 = as.integer((post_summary[1:5,4] <= params) & (params <= post_summary[1:5,8]))
        
        ## Armazenar resultados
        kg_mc_bs[i,1,j] = post_mean[1]; kg_mc_bs[i,2,j] = post_sd[1]; kg_mc_bs[i,3,j] = rbias[1]; kg_mc_bs[i,4,j] = cpost95[1]
        kg_mc_bs[i,5,j] = post_mean[2]; kg_mc_bs[i,6,j] = post_sd[2]; kg_mc_bs[i,7,j] = rbias[2]; kg_mc_bs[i,8,j] = cpost95[2]
        kg_mc_bs[i,9,j] = post_mean[3]; kg_mc_bs[i,10,j] = post_sd[3]; kg_mc_bs[i,11,j] = rbias[3]; kg_mc_bs[i,12,j] = cpost95[3]
        kg_mc_bs[i,13,j] = post_mean[4]; kg_mc_bs[i,14,j] = post_sd[4]; kg_mc_bs[i,15,j] = rbias[4]; kg_mc_bs[i,16,j] = cpost95[4]
        kg_mc_bs[i,17,j] = post_mean[5]; kg_mc_bs[i,18,j] = post_sd[5]; kg_mc_bs[i,19,j] = rbias[5]; kg_mc_bs[i,20,j] = cpost95[5]
        
        kg_mc_cens[i,,j] = 1 - (sum(dados.reg.kg[,2]) / n_amostral[j])
        
        if(i == n_replicas) {
          write.csv2(kg_mc_bs[,,j], file = paste("2kg_mc_bs_", "rep_", n_replicas, "_n_", n_amostral[j], ".csv", sep = ""))
          write.csv2(kg_mc_cens[,,j], file = paste("2kg_mc_cens_", "rep_", n_replicas, "_n_", n_amostral[j], ".csv", sep = ""))
        }
        
      }, error = function(e) {
        message("Erro em alguma parte da iteração: ", conditionMessage(e))
        erro_ao_rodar <<- TRUE
      })
      
      if (!erro_ao_rodar) break  # se deu tudo certo, sai do repeat
    }
  }
}







































