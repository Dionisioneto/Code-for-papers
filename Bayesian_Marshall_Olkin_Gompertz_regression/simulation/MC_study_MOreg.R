### -----------
###Monte Carlo Study 
### Bayesian regression MO Gompertz
### Generate data and estimate via STAN
### -----------

library(pacman)
p_load(survival,rstan)

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
gen.cure.reg.mog = function(n,as,bs,xa,xb,l){#, max_cen=20){
  finv.mog = function(t,alpha,beta,lambda,unif){Ftmo_gompertz(t=t,alpha=alpha,beta=beta,lambda=lambda) - unif}
  ger.mog = function(alpha,beta,lambda,unif){t=uniroot(finv.mog,c(0,1000),tol=0.01,alpha=alpha,beta=beta,lambda=lambda,unif=unif);return(t$root)}
  
  ## individual parameters
  
  a = cbind(1,xa) %*% as
  b = exp(cbind(1,xb) %*% bs)
  p_base = exp(b/a)
  p_cure = (l*p_base)/(1-((1-l)*p_base))
  
  t = rep(NA,n)
  
  rmi = rbinom(n=n,size=1,prob=1-p_cure)
  uni = runif(n=n,0,max=1-p_cure)
  
  for(i in 1:n){
    t[i] = ifelse(rmi[i]==0, Inf, ger.mog(alpha=a[i,],beta=b[i,],lambda=l,unif=uni[i]))
  }
  
  t_finite = ifelse(t==Inf,0,t)
  u2 = runif(n=n,0,max(t_finite))
  #u2 = runif(n=n,0,max_cen)
  t2 = pmin(t,u2); delta = ifelse(t<u2,1,0)
  
  df = data.frame(time=t2,delta=delta,xa=xa,xb=xb)
  return(df)
}


cod_regmog = "
data {
  int<lower=0> N; // Número de observações
  int<lower=0> p; // Número de covariáveis para alpha
  int<lower=0> s; // Número de covariáveis para beta

  vector[N] time;  // Tempo de falha
  vector[N] delta; // Indicador de censura
  
  matrix[N, p] xa; // Matriz de covariáveis para alpha (N x p)
  matrix[N, s] xb; // Matriz de covariáveis para beta (N x s)
}

parameters {
  vector[p] as; // Coeficientes para alpha
  vector[s] bs; // Coeficientes para beta
  real<lower=0> lambda; 
}

transformed parameters {
  vector[N] alpha = xa * as; // Cálculo de alpha
  vector[N] beta = exp(xb * bs); // Cálculo de beta
}

model {
  // Priori para os coeficientes em alpha
  for (j in 1:p) {
    as[j] ~ normal(0, 10);
  }
  
  // Priori para os coeficientes em beta
  for (j in 1:s) {
    bs[j] ~ normal(0, 10);
  }
  
  // Priori para lambda
  lambda ~ uniform(0,50);
  
  // Definição da função log-verossimilhança
  for (i in 1:N) {
    target += delta[i] * log(beta[i] * lambda * exp((beta[i] - beta[i] * exp(alpha[i] * time[i])) / alpha[i] + alpha[i] * time[i])) -
              delta[i] * log((1 - ((1-lambda) * exp((beta[i] - beta[i] * exp(alpha[i] * time[i])) / alpha[i])))^2) +
              (1 - delta[i]) * (log(lambda * exp(-(beta[i] / alpha[i]) * (exp(alpha[i] * time[i]) - 1)))) -
              (1 - delta[i]) * (log(1 - ((1 - lambda) * exp((-beta[i] / alpha[i]) * (exp(alpha[i] * time[i]) - 1)))));
  }
}
"


## Transcrever o código escrito para um file stan 
writeLines(cod_regmog, con = "cod_regmog.stan")



## Choosing parameters and setting tables 

as0 = c(-0.4, 0.1, -0.1)
bs0 = c(-2, 2, 0.2)
l0 =0.6

params = c(as0,bs0,l0)

## Tables to save simulation results
niter=1
ns = c(100, 300, 500, 700, 1000, 2000)

simres = array(data=NA,dim=c(niter,4*length(params),length(ns)))
censprop = array(data=NA,dim=c(niter,1,length(ns)))
#rhatvals = array(data=NA,dim=c(niter,1,length(ns)))


for(j in 1:length(ns)){
  for(k in 1:niter){
    print(paste("Iteration ", k, "for sample size: ", ns[j]))
    
    x1 = rbinom(n=ns[j],size=1,prob=0.5)
    x2 = runif(n=ns[j],0,1)
    xs0 = cbind(x1,x2)
    
    dataset =  gen.cure.reg.mog(n=ns[j],as=as0,bs=bs0,xa=xs0,xb=xs0,l=l0)
    
    # organizanizing the data for stan
    
    data_cr_mogreg = list(N = dim(xs0)[1],
                         p = ncol(xs0)+1,
                         s = ncol(xs0)+1,
                         time = dataset$time,
                         delta = dataset$delta,
                         xa = cbind(1,xs0),
                         xb = cbind(1,xs0))
    
    
    m4 = stan(file   = "cod_regmog.stan", 
              data   =  data_cr_mogreg,
              warmup = 5000,                 
              iter   = 8000,
              chains = 2,
              seed   = 4,
              cores  = getOption("mc.cores",1)) 
    
    summary(m4)$summary[1:length(params), ]
    params
    
    postmean=summary(m4)$summary[1:length(params), "mean"]
    postsd=summary(m4)$summary[1:length(params), "sd"]
    ci = cbind(summary(m4)$summary[1:length(params), "2.5%"],summary(m4)$summary[1:length(params), "97.5%"])
    cp = params >= ci[,1] & params <= ci[,2] ## coverage probability
    
    bias = (params - postmean)
    
    censoring = 1-mean(dataset$delta)

    #mean,sd,bias,cp
    simres[k,1:length(params),j]=postmean
    simres[k,(length(params)+1):(2*length(params)),j]=postsd
    simres[k,(2*length(params)+1):(3*length(params)),j]=cp
    simres[k,(3*length(params)+1):(4*length(params)),j]=bias
    
    censprop[k,1,j] = censoring
    #rhatvals[k,1,j] = rhat
  }
  write.csv2(simres[,,j], file = paste("simres_", ns[j], ".csv", sep = ""))
  write.csv2(censprop[k,1,j], file = paste("censprop_", ns[j], ".csv", sep = ""))
  
}








