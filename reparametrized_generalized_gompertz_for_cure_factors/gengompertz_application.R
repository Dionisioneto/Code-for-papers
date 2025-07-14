## ---
## Bayesian estimation [K Gompertz ]
## ---

rm(list = ls())
setwd("C:/Users/dioni/OneDrive - University of São Paulo/Doutorado em Estatística/03_semestres/2025.1/topicos_avancados_em_pesquisa_2/application")

utero = read.csv2("utero.csv"); utero = utero[,-1]

head(utero)
dim(utero)

## gráfico da curva de sobrevivência

library(survival)
library(survminer)

kmfit = survfit(Surv(time,delta)~1,data=utero)

ggsurvplot(kmfit, data = utero,
           conf.int = TRUE,       # Intervalo de confiança
           pval = TRUE,           # Valor-p do teste log-rank
           legend.title = "Grupo",
           xlab = "\n Time (years)",
           ylab = "Survival Probability \n ",
           color="black")




library(rstan)

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
    vector[N] p = 1/(1 + exp(-X * bs));
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
      real beta_i = alpha * log(1 - pow(1 - p[i], 1.0 / psi));
  
      target += delta[i] * (
        log(psi * beta_i * exp(alpha * time[i]) *
        exp((beta_i - beta_i * exp(alpha * time[i])) / alpha) *
        pow(1.0 - exp((beta_i - beta_i * exp(alpha * time[i])) / alpha), psi - 1) 
      )) +
      (1 - delta[i]) * (
        log(1.0 - pow(1.0 - exp((beta_i - beta_i * exp(alpha * time[i])) / alpha), psi)));
    }
  }
  
  
  // Para computar métricas de comparação
  generated quantities {
    vector[N] log_lik;
    for (i in 1:N) {
      real beta_i =  alpha * log(1 - pow(1 - p[i], 1.0 / psi));
      
      log_lik[i] =  delta[i] * (
        log(psi * beta_i * exp(alpha * time[i]) *
        exp((beta_i - beta_i * exp(alpha * time[i])) / alpha) *
        pow(1.0 - exp((beta_i - beta_i * exp(alpha * time[i])) / alpha), psi - 1))) +
        (1 - delta[i]) * (
        log(1.0 - pow(1.0 - exp((beta_i - beta_i * exp(alpha * time[i])) / alpha), psi)));
      
    }
  }
"

## Transcrever o código escrito para um file stan 
writeLines(cod_regkg, con = "cod_regkg.stan")

# Data and stan function

#cols_cov = c("catidade", "cirurgia", "quimio", "hormonio") #[WORKING]
cols_cov = c("catidade", "RECDIST","cirurgia", "quimio", "hormonio") # retirar hormonio

X_cov = cbind(1,utero[,cols_cov])

data_regkg = list(N = dim(utero)[1],
                  q = dim(cbind(1,utero[,cols_cov]))[2],
                  time = utero[,"time"],
                  delta = utero[,"delta"],
                  X = X_cov)

# Run the Bayesian model
inic=Sys.time()
regkgfit = stan(file = 'cod_regkg.stan', data = data_regkg,
                chains = 4, iter=5000, warmup=1000,cores=4, seed = 1234)
end=Sys.time()
end-inic

# It took 9 mins to run using 4 cores

# summary
my_parameters = rownames(summary(regkgfit)$summary)[1:8]
sumprint = summary(regkgfit)$summary[my_parameters,]

round(sumprint[,c(1,3,4,8)],4)

write.csv2(sumprint, file = "summprint_oficial.csv")
dir()


## Computing model selection measurements

## [WAIC] 
library(loo)
loglik_regkg = extract_log_lik(regkgfit, parameter_name = "log_lik")
waic_regkg = waic(loglik_regkg)
waic_regkg

# [LOOIC] (alternativa melhor que WAIC na maioria dos casos)
loo_regkg = loo(loglik_regkg)
loo_regkg$estimates[3,1]

#par(mfrow=c(1,2))

plot(loo_regkg, main = "", xlab = "index") ## uma observação influente identificada

which(loo_regkg$diagnostics$pareto_k>0.7)


library(ggplot2)
library(loo)

# Supondo que você já rodou o loo:
# loo_regkg = loo(loglik_regkg)

# Extraindo os valores do pareto_k
pareto_k = loo_regkg$diagnostics$pareto_k

# Criando um dataframe para o ggplot
df_loo = data.frame(
  obs = 1:length(pareto_k),
  pareto_k = pareto_k
)

# Gráfico com destaque e rótulo para observações influentes
library(ggrepel)
ggplot(df_loo, aes(x = obs, y = pareto_k)) +
  geom_point(aes(color = pareto_k > 0.7), size = 2, show.legend = FALSE) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  geom_text_repel(
    data = subset(df_loo, pareto_k > 0.7),
    aes(label = obs),
    size = 3,
    color = "red",
    max.overlaps = Inf
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "red") +
  labs(
    title = "",
    x = "\n Index",
    y = "Pareto shape k \n"
  ) +
  coord_cartesian(ylim = c(-0.5, 1)) +
  theme_minimal()



## [LPLM]

## cpo
lik = exp(extract(regkgfit)$log_lik)

cpo = sapply(1:ncol(lik), function(i) {
  inv_mean <- mean(1 / lik[, i])
  1 / inv_mean
})

#plot(1:length(cpo),log(cpo), ylab="Log (CPO)", xlab="Data point")

-2*sum(log(cpo))

# Supondo que 'cpo' já está no ambiente
df_cpo = data.frame(
  obs = 1:length(cpo),
  log_cpo = log(cpo)
)

ggplot(df_cpo, aes(x = obs, y = log_cpo)) +
  geom_point(size = 2, color = "black") +
  labs(
    title = "",
    x = "\n Index",
    y = "Log (CPO) \n"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal()





## [DIC]

log_lik_regkgfit = extract(regkgfit)$log_lik  # dim: (n_iter x N)
mean_dev_regkgfit = -2 * mean(rowSums(log_lik_regkgfit))

posterior_means = summary(regkgfit)$summary
alpha_bar = posterior_means["alpha", "mean"]
psi_bar = posterior_means["psi", "mean"]
bs_bar = posterior_means[as.vector(paste0("bs[", 1:(length(cols_cov)+1), "]")),"mean"]

# Calculando p médio
p_bar = 1 / (1 + exp(-as.matrix(X_cov) %*% bs_bar))

# Calcular log-verossimilhança com os parâmetros médios
N = nrow(utero)
loglik_bar_regkgfit = numeric(N)

N = dim(utero)[1]
time = utero[,"time"]
delta = utero[,"delta"]


for (i in 1:N) {
  beta_i = alpha_bar * log(1-(1-p_bar[i])^(1/psi_bar))
  
  loglik_bar_regkgfit[i] = delta[i]*log(psi_bar * beta_i * exp(alpha_bar * time[i]) *
                                 exp((beta_i - beta_i * exp(alpha_bar * time[i])) / alpha_bar) *
                                  (1.0 - exp((beta_i - beta_i * exp(alpha_bar * time[i])) / alpha_bar))^(psi_bar-1)) +
                          (1-delta[i])*(log(1-(1 - exp((beta_i - beta_i * exp(alpha_bar * time[i])) / alpha_bar))^(psi_bar)))
}


deviance_loglik_bar_regkgfit = -2 * sum(loglik_bar_regkgfit)

DIC = (2*mean_dev_regkgfit) - deviance_loglik_bar_regkgfit
DIC


## Calculando os resíduos martingales e deviance

surv_kgomp = function(alpha,psi,betas,time,X){

  pcure = 1 / (1 + exp(-as.matrix(X) %*% betas))
  beta_i = alpha * log(1-(1-pcure)^(1/psi))
  st= 1-(1 - exp((beta_i - beta_i * exp(alpha * time)) / alpha))^(psi)
  return(st)
}

surv=surv_kgomp(alpha=alpha_bar,psi=psi_bar,betas=bs_bar,time=utero$time,X=X_cov)

rmi = delta + log(surv)

dev=sign(rmi)*sqrt(-2*(rmi+delta*log(delta-rmi)))


par(mfrow=c(1,2))
plot(rmi, xlab= "Data point", ylab="Martingale residulas")

plot(dev, xlab= "Data point", ylab="Deviance residulas",ylim=c(-4,4))
abline(h=3,col='red',lwd=2)
abline(h=-3,col='red',lwd=2)

dev.off()

df_rmi = data.frame(
  obs = 1:length(rmi),
  rmi = rmi
)

ggplot(df_rmi, aes(x = obs, y = rmi)) +
  geom_point(size = 2, color = "black") +
  labs(
    title = "",
    x = "\n Index",
    y = "Martingale residual \n"
  ) +
  theme_minimal()


# Criação do dataframe correto
df_dev = data.frame(
  obs = 1:length(dev),
  rmi = dev
)

# Gráfico com identificação dos pontos críticos
ggplot(df_dev, aes(x = obs, y = rmi)) +
  geom_point(aes(color = (rmi < -3 | rmi > 3)), size = 2, show.legend = FALSE) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  geom_text_repel(
    data = subset(df_dev, rmi < -3 | rmi > 3),
    aes(label = obs),
    size = 3,
    color = "red",vjust = -1.2,
    max.overlaps = Inf
  ) +
  geom_hline(yintercept = -3, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 3, linetype = "dashed", color = "red") +
  coord_cartesian(ylim = c(-4, 4)) +
  labs(
    title = "",
    x = "\n Index",
    y = "Deviance residual\n"
  ) +
  theme_minimal()





## Estiamtivas sem pontos influentes ou ouliers
 
which(loo_regkg$diagnostics$pareto_k>0.7) ## verificando o ponto que teve alto pareto k 
which(dev>3) ## verificando o ponto que teve alta deviance


utero2_pareto = utero[-which(loo_regkg$diagnostics$pareto_k>0.7),]
X_cov2 = cbind(1,utero2_pareto[,cols_cov])

data_regkg2 = list(N = dim(utero2_pareto)[1],
                  q = dim(cbind(1,utero2_pareto[,cols_cov]))[2],
                  time = utero2_pareto[,"time"],
                  delta = utero2_pareto[,"delta"],
                  X = X_cov2)

# Run the Bayesian model
inic=Sys.time()
regkgfit2 = stan(file = 'cod_regkg.stan', data = data_regkg2,
                chains = 4, iter=5000, warmup=1000,cores=4,
                seed = 1234)
end=Sys.time()
end-inic

loglik_regkg2 = extract_log_lik(regkgfit2, parameter_name = "log_lik")
loo_regkg2 = loo(loglik_regkg2)
loo_regkg2$estimates[3,1]


plot(loo_regkg2, main = "", xlab = "index") ## uma observação influente identificada


round(summary(regkgfit2)$summary[1:8,c(1,3,4,8)],4)



## 3
utero3_pareto = utero[-which(dev>3),]
X_cov3 = cbind(1,utero3_pareto[,cols_cov])

data_regkg3 = list(N = dim(utero3_pareto)[1],
                   q = dim(cbind(1,utero3_pareto[,cols_cov]))[2],
                   time = utero3_pareto[,"time"],
                   delta = utero3_pareto[,"delta"],
                   X = X_cov3)

# Run the Bayesian model
inic=Sys.time()
regkgfit3 = stan(file = 'cod_regkg.stan', data = data_regkg3,
                 chains = 4, iter=5000, warmup=1000,cores=4,
                 seed = 1234)
end=Sys.time()
end-inic


my_parameters = rownames(summary(regkgfit)$summary)[1:8]
sumprint3 = summary(regkgfit3)$summary[my_parameters,]

round(sumprint3[,c(1,3,4,8)],4)
round(summary(regkgfit3)$summary[1:8,c(1,3,4,8)],4)



write.csv2(round(sumprint3[,c(1,3,4,8)],4), file = "summprint3_oficial.csv")
dir()


loglik_regkg3 = extract_log_lik(regkgfit3, parameter_name = "log_lik")
loo_regkg3 = loo(loglik_regkg3)
loo_regkg3$estimates[3,1]


plot(loo_regkg3, main = "", xlab = "index") ## uma observação influente identificada

## sem as duas observações

utero4_pareto = utero[-c(which(dev>3),which(loo_regkg$diagnostics$pareto_k>0.7)),]
X_cov4 = cbind(1,utero4_pareto[,cols_cov])

data_regkg4 = list(N = dim(utero4_pareto)[1],
                   q = dim(cbind(1,utero4_pareto[,cols_cov]))[2],
                   time = utero4_pareto[,"time"],
                   delta = utero4_pareto[,"delta"],
                   X = X_cov4)

# Run the Bayesian model
inic=Sys.time()
regkgfit4 = stan(file = 'cod_regkg.stan', data = data_regkg4,
                 chains = 4, iter=5000, warmup=1000,cores=4,
                 seed = 1234)
end=Sys.time()
end-inic


my_parameters = rownames(summary(regkgfit)$summary)[1:8]
sumprint4 = summary(regkgfit4)$summary[my_parameters,]

round(sumprint4[,c(1,3,4,8)],4)

write.csv2(round(sumprint4[,c(1,3,4,8)],4), file = "summprint4_oficial.csv")
dir()

loglik_regkg4 = extract_log_lik(regkgfit4, parameter_name = "log_lik")
loo_regkg4 = loo(loglik_regkg4)
loo_regkg4$estimates[3,1]


plot(loo_regkg4, main = "", xlab = "index") ## uma observação influente identificada

## Calcular medidas de diagnóstico

# [LOOIC] (alternativa melhor que WAIC na maioria dos casos)
loo_regkg4 = loo(loglik_regkg4)
loo_regkg4$estimates[3,1]





## [LPLM]

## cpo
lik4 = exp(extract(regkgfit4)$log_lik)

cpo4 = sapply(1:ncol(lik4), function(i) {
  inv_mean <- mean(1 / lik4[, i])
  1 / inv_mean
})

plot(1:length(cpo4),log(cpo4), ylab="Log (CPO)", xlab="Data point")

-2*sum(log(cpo4))



# DIC
## [DIC]

log_lik_regkgfit4 = extract(regkgfit4)$log_lik  # dim: (n_iter x N)
mean_dev_regkgfit4 = -2 * mean(rowSums(log_lik_regkgfit4))

posterior_means4 = summary(regkgfit4)$summary
alpha_bar4 = posterior_means["alpha", "mean"]
psi_bar4 = posterior_means["psi", "mean"]
bs_bar4 = posterior_means[as.vector(paste0("bs[", 1:(length(cols_cov)+1), "]")),"mean"]

# Calculando p médio
p_bar4 = 1 / (1 + exp(-as.matrix(X_cov) %*% bs_bar4))

# Calcular log-verossimilhança com os parâmetros médios
N = nrow(utero)
loglik_bar_regkgfit4 = numeric(N)

N = dim(utero)[1]
time = utero[,"time"]
delta = utero[,"delta"]


for (i in 1:N) {
  beta_i = alpha_bar4 * log(1-(1-p_bar4[i])^(1/psi_bar4))
  
  loglik_bar_regkgfit4[i] = delta[i]*log(psi_bar4 * beta_i * exp(alpha_bar4 * time[i]) *
                                          exp((beta_i - beta_i * exp(alpha_bar4 * time[i])) / alpha_bar4) *
                                          (1.0 - exp((beta_i - beta_i * exp(alpha_bar4 * time[i])) / alpha_bar4))^(psi_bar4-1)) +
    (1-delta[i])*(log(1-(1 - exp((beta_i - beta_i * exp(alpha_bar4 * time[i])) / alpha_bar4))^(psi_bar4)))
}


deviance_loglik_bar_regkgfit4 = -2 * sum(loglik_bar_regkgfit4)

DIC = (2*mean_dev_regkgfit4) - deviance_loglik_bar_regkgfit4
DIC






## [Modelo de regressão sob a distribuição Gompertz]

cod_g = "
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
    }
    
    transformed parameters {
      vector[N] p = 1 / (1 + exp(-X * bs)); // Probabilidade de cura
    }
    
    model {
      // Priori para os coeficientes em p
      bs ~ normal(0, 10);
    
      // Priori para alpha
      alpha ~ normal(-1, 10);
    
      // Log-verossimilhança para psi = 1
      for (i in 1:N) {
        real beta_i =  alpha * log(p[i]);
        
        if (delta[i] == 1) {
          target+= log(beta_i * exp(alpha * time[i]) *
          exp((beta_i - beta_i * exp(alpha * time[i])) / alpha));
        } else {
          target+= log(exp((beta_i - beta_i * exp(alpha * time[i])) / alpha));
        }
      }
    }
    
    generated quantities {
    vector[N] log_lik;
    for (i in 1:N) {
        real beta_i =  alpha * log(p[i]);
        
        if (delta[i] == 1) {
          log_lik[i] = log(beta_i * exp(alpha * time[i]) *
          exp((beta_i - beta_i * exp(alpha * time[i])) / alpha));
        } else {
          log_lik[i] = log(exp((beta_i - beta_i * exp(alpha * time[i])) / alpha));
        }
      }
  }
"

## Transcrever o código escrito para um file stan 
writeLines(cod_g, con = "cod_g.stan")


## Run the Bayesian model

inic_reggfit = Sys.time()
reggfit = stan(file = 'cod_g.stan', data = data_regkg,
               chains = 4, iter=5000, warmup=1000, cores = 4)
end_reggfit = Sys.time()

end_reggfit-inic_reggfit

summary(reggfit)$summary[1:7,c(1,3,4,8)]


## Computando os valores de comparação

## [WAIC] 
library(loo)
loglik_regg = extract_log_lik(reggfit, parameter_name = "log_lik")
waic_regg = waic(loglik_regg)
#waic_regg

# [LOOIC] (alternativa melhor que WAIC na maioria dos casos)
loo_regg = loo(loglik_regg)
loo_regg
loo_regg$estimates[3,1]

#loo_compare(loo_regg,loo_regkg)

## [LPLM]

## cpo
lik_regg = exp(extract(reggfit)$log_lik)

cpo_regg = sapply(1:ncol(lik_regg), function(i) {
  inv_mean <- mean(1 / lik_regg[, i])
  1 / inv_mean
})

plot(1:length(cpo_regg),log(cpo_regg))

-2*sum(log(cpo_regg))

## [DIC]

## (i) calcular dev media = -2 * mean(rowSums(log_lik_regkgfit))

loglik_regg = extract(reggfit)$log_lik
meandev_regg = -2*mean(rowSums(loglik_regg))

sum_regg = summary(reggfit)$summary
alpha_bar_regg = sum_regg["alpha", "mean"]
bs_bar_regg = sum_regg[as.vector(paste0("bs[", 1:(length(cols_cov)+1), "]")),"mean"]

# Calculando p médio
p_bar_regg = 1 / (1 + exp(-as.matrix(X_cov) %*% bs_bar_regg))

# Calcular log-verossimilhança com os parâmetros médios
N = nrow(utero)
loglik_bar_regg = numeric(N)

time = utero[,"time"]
delta = utero[,"delta"]

for (i in 1:N) {
  beta_i = alpha_bar_regg * log(p_bar_regg[i])
  
  if (delta[i] == 1) {
    loglik_bar_regg[i] = log(beta_i * exp(alpha_bar_regg * time[i]) *
                                 exp((beta_i - beta_i * exp(alpha_bar_regg * time[i])) / alpha_bar_regg));
  } else {
    loglik_bar_regg[i] = log(exp((beta_i - beta_i * exp(alpha_bar_regg * time[i])) / alpha_bar_regg));
  }
}

deviance_bar_regg = -2 * sum(loglik_bar_regg)

DIC_regg = (2*meandev_regg) - deviance_bar_regg
DIC_regg



## plotar a curva de sobrevivencia














