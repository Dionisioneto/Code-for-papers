### -----------
### Monte Carlo Study 
### Bayesian regression MO Gompertz
### Generate data and estimate via STAN
### -----------
# 
# library(pacman)
# p_load(survival,rstan)

library(rstan)

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

## To stan file 
writeLines(cod_regmog, con = "cod_regmog.stan")

## test
n0 = c(50, 1000, 300, 500, 800, 1000, 2000)

## scenario 1 (real values for parameters, average cure and average delta)
a0=-0.517; l0 = 0.485; bs0 = c(-2, 0.701, 2.487) 

## scenario 2 (high cure and delta)
#a0=-0.2; l0 = 0.5; bs0 = c(-2, 0.701, 2.487) 

## scenario 3 (low cure and delta)
#a0=-2; l0 = 1.5; bs0 = c(-2, 0.701, 2.487)  

params = c(a0,l0,bs0)

# j=4
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
# 
# m4 = stan(file   = "cod_regmog.stan", 
#           data   =  data_cr_mogreg,
#           warmup = 3000,                 
#           iter   = 5000,
#           chains = 1,
#           seed   = 4,
#           cores  = getOption("mc.cores",1)) 
# 
# traceplot(m4, pars = c("alpha", "lambda", "bs"))
# summary(m4, pars = c("alpha", "lambda", "bs"))$summary[,c("mean", "sd", "2.5%", "97.5%","n_eff", "Rhat")]
# params


params = c(a0,l0,bs0)

niter=1000

simres = array(data=NA,dim=c(niter,4*length(params),length(n0)))
censprop = array(data=NA,dim=c(niter,1,length(n0)))
hatvals = array(data=NA,dim=c(niter,4*length(params),length(n0)))

comp_m4 = rstan::stan_model(file = "cod_regmog.stan")

for(j in 1:length(n0)){
  for(k in 1:niter){
    print(paste("Iteration ", k, "for sample size: ", n0[j]))
    
    stages = sample(x = c(1,2,3), size = n0[j], replace = T, prob = c(1/3, 1/3, 1/3))
    W_mat = cbind(stage2 = as.numeric(stages == 2), stage3 = as.numeric(stages == 3))

    dataset = gen.cure.reg.mog(n=n0[j],a=a0, bs = bs0, xb = W_mat, l = l0)
    
    data_cr_mogreg = list(N = n0[j],
                          s = dim(W_mat)[2]+1,
                          time = dataset$time,
                          delta = dataset$delta,
                          xb = cbind(1,W_mat))
    

    m4 = rstan::sampling(object = comp_m4, 
                         data   = data_cr_mogreg,
                         warmup = 5000,                  
                         iter   = 8000,
                         chains = 2,
                         cores  = getOption("mc.cores", 2),
                         refresh = 0) # refresh=0 limpa o output no console do cluster
    
    sum_m4 = summary(m4, pars = c("alpha", "lambda", "bs"))$summary[,c("mean", "sd", "2.5%", "97.5%","n_eff", "Rhat")]
    
    postmean=sum_m4[, "mean"]
    postsd=sum_m4[, "sd"]
    cp = params >= sum_m4[, "2.5%"] & params <= sum_m4[, "97.5%"] ## coverage probability
    
    bias = (params - postmean)
    
    censoring = 1-mean(dataset$delta)

    #mean,sd,bias,cp, Rhat
    simres[k,1:length(params),j]=postmean
    simres[k,(length(params)+1):(2*length(params)),j]=postsd
    simres[k,(2*length(params)+1):(3*length(params)),j]=cp
    simres[k,(3*length(params)+1):(4*length(params)),j]=bias
    
    hatvals[k,1:length(params),j]=as.numeric(sum_m4[,"Rhat"])
    
    censprop[k,1,j] = censoring

  }
  write.csv2(simres[,,j], file = paste("simres_", n0[j], ".csv", sep = ""))
  write.csv2(hatvals[,,j], file = paste("Rhatvals_", n0[j], ".csv", sep = ""))
  write.csv2(censprop[,1,j], file = paste("censprop_", n0[j], ".csv", sep = ""))
  
}








