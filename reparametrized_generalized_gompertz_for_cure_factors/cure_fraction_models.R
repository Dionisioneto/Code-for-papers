## ------ 
## Script aplication 
## testing cure fraction models
## utero
## ------

rm(list = ls())
setwd("C:/Users/dioni/OneDrive - University of São Paulo/Doutorado em Estatística/03_semestres/2025.1/topicos_avancados_em_pesquisa_2/application")

library(rstan)

utero = read.csv2("utero.csv"); utero = utero[,-1]

head(utero)
dim(utero)


## [Comparar com outras estuturas]
## (a) Modelo de mistura padrão Weibull
## (b) Modelo base da Gompertz


## (a) [Modelo de mistrua padrão da weibull]

cod_weibull = "
    data {
      int<lower=1> N;                  // número de indivíduos
      vector<lower=0>[N] time;         // tempos observados
      int<lower=0, upper=1> delta[N];  // indicador de evento (1=evento, 0=censura)
      int<lower=1> K;                  // número de covariáveis para cura
      matrix[N, K] Z;                  // matriz de covariáveis para cura
    }
    
    parameters {
      vector[K] betas;                 // coeficientes do modelo de cura
      real<lower=0> lambda;           // parâmetro de escala da Weibull
      real<lower=0> alpha;            // parâmetro de forma da Weibull
    }
    
    model {
      // Priori
      betas ~ normal(0, 10);
      lambda ~ gamma(0.01, 0.01);           // prior fraca
      alpha ~ gamma(0.01, 0.01);
    
      for (i in 1:N) {
        real pi_i = inv_logit(Z[i] * betas);        // probabilidade de cura
        real S0 = exp(-(time[i]/lambda)^alpha); // sobrevivência dos suscetíveis
        real f0 = (alpha/lambda)*(time[i]/lambda)^(alpha-1)*exp(-(time[i]/lambda)^(alpha)); // densidade dos suscetíveis
    
        if (delta[i] == 1) {
          // evento observado: apenas suscetíveis podem falhar
          target += log(pi_i*f0);
        } else {
          // censura: mistura entre curados e suscetíveis sobreviventes
          target += log(pi_i + (1 - pi_i) * S0);
        }
      }
    }
    
    // Para computar métricas de comparação
    generated quantities {
    vector[N] log_lik;
    for (i in 1:N) {
      real pi_i = inv_logit(Z[i] * betas);        // probabilidade de cura
        real S0 = exp(-(time[i]/lambda)^alpha); // sobrevivência dos suscetíveis
        real f0 = (alpha/lambda)*(time[i]/lambda)^(alpha-1)*exp(-(time[i]/lambda)^(alpha)); // densidade dos suscetíveis
      
      if (delta[i] == 1) {
        log_lik[i] = log(pi_i*f0);
      } else {
        log_lik[i] = log(pi_i + (1 - pi_i) * S0);
      }
    }
  }
"

cols_cov = c("catidade", "RECDIST","cirurgia", "quimio", "hormonio") # retirar hormonio

## Transcrever o código escrito para um file stan 
writeLines(cod_weibull, con = "cod_weibull.stan")

weibull_data = list(
  N = dim(utero)[1],
  time = utero[,"time"],
  delta = utero[,"delta"],
  K = ncol(cbind(1,utero[,cols_cov])),
  Z = cbind(1,utero[,cols_cov])
)

weibullfit = stan("cod_weibull.stan", data = weibull_data,
                  chains = 4, iter = 5000, warmup = 1000, cores=4, seed=1234)


sum_wb = round(summary(weibullfit)$summary[1:8,c(1,3,4,8)],4)
sum_wb
write.csv2(sum_wb, file = "summprint_wb.csv")


weibull_samples = extract(weibullfit)

# plot(weibull_samples$betas[,1], type='l')
# plot(weibull_samples$betas[,2], type='l')
# plot(weibull_samples$betas[,3], type='l')
# plot(weibull_samples$betas[,4], type='l')
# plot(weibull_samples$betas[,5], type='l')
# 
# plot(weibull_samples$lambda, type='l')
# plot(weibull_samples$alpha, type='l')


## Calculando métricas de comparação

## [WAIC] 
library(loo)
loglik_wb = extract_log_lik(weibullfit, parameter_name = "log_lik")
waic_wb = waic(loglik_wb)
waic_wb


# [LOOIC] (alternativa melhor que WAIC na maioria dos casos)
loo_wb = loo(loglik_wb)
loo_wb$estimates[3,1]


plot(loo_wb, main="")


library(ggplot2)
library(loo)

# Supondo que você já rodou o loo:
# loo_regkg = loo(loglik_regkg)

# Extraindo os valores do pareto_k
pareto_k = loo_regkg$diagnostics$pareto_k

# Criando um dataframe para o ggplot
df_loo_wb = data.frame(
  obs = 1:length(loo_wb$diagnostics$pareto_k),
  pareto_k = loo_wb$diagnostics$pareto_k
)

# Gráfico com destaque e rótulo para observações influentes
library(ggrepel)
ggplot(df_loo_wb, aes(x = obs, y = pareto_k)) +
  geom_point(aes(color = pareto_k > 0.7), size = 2, show.legend = FALSE) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  geom_text_repel(
    data = subset(df_loo_wb, pareto_k > 0.7),
    aes(label = obs),
    size = 3,
    color = "red",
    max.overlaps = Inf
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 1.0, linetype = "twodash", color = "red") +
  labs(
    title = "",
    x = "\n Index",
    y = "Pareto shape k \n"
  ) +
  coord_cartesian(ylim = c(-0.5, 2.5)) +
  theme_minimal()





## [LPLM]

lik_wb = exp(extract(weibullfit)$log_lik)

cpo_wb = sapply(1:ncol(lik_wb), function(i) {
  inv_mean <- mean(1 / lik_wb[, i])
  1 / inv_mean
})

plot(1:length(cpo_wb),log(cpo_wb), ylab="Log (CPO)", xlab="Data point")

-2*sum(log(cpo_wb))


# Supondo que 'cpo' já está no ambiente
df_cpo_wb = data.frame(
  obs = 1:length(cpo_wb),
  log_cpo = log(cpo_wb)
)

ggplot(df_cpo_wb, aes(x = obs, y = log_cpo)) +
  geom_point(size = 2, color = "black") +
  labs(
    title = "",
    x = "\n Index",
    y = "Log (CPO) \n"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal()









# Calculando o DIC


log_lik_wbfit = extract(weibullfit)$log_lik  # dim: (n_iter x N)
mean_dev_wbfit = -2 * mean(rowSums(log_lik_wbfit))



post_meanswb = summary(weibullfit)$summary
lambda_bar_wb = post_meanswb["lambda", "mean"]
alpha_bar_wb = post_meanswb["alpha", "mean"]
bs_bar_wb = post_meanswb[as.vector(paste0("betas[", 1:(length(cols_cov)+1), "]")),"mean"]

# Calculando p médio
p_bar_wb = 1 / (1 + exp(-as.matrix(X_cov) %*% bs_bar_wb))

# Calcular log-verossimilhança com os parâmetros médios
N = nrow(utero)
loglik_bar_wb = numeric(N)

time = utero[,"time"]
delta = utero[,"delta"]

for (i in 1:N) {
  S0 = exp(-(time[i]/lambda_bar_wb)^alpha_bar_wb)
  f0 = (alpha_bar_wb/lambda_bar_wb)*(time[i]/lambda_bar_wb)^(alpha_bar_wb-1)*exp(-(time[i]/lambda_bar_wb)^(alpha_bar_wb))
  
  if (delta[i] == 1) {
    loglik_bar_wb[i] = log(p_bar_wb[i]*f0)
  } else {
    loglik_bar_wb[i] = log(p_bar_wb[i] + (1 - p_bar_wb[i]) * S0);
  }
}

deviance_bar_wb = -2 * sum(loglik_bar_wb)

DIC_wb = (2*mean_dev_wbfit) - deviance_bar_wb
DIC_wb

## analise dos residuos martingales weibull  


surv_wb = exp(-(time/lambda_bar_wb)^alpha_bar_wb)
rmi_wb = delta + log(surv_wb)
dev_wb = sign(rmi_wb)*sqrt(-2*(rmi_wb+delta*log(delta-rmi_wb)))


par(mfrow=c(1,2))
plot(rmi_wb, xlab= "Data point", ylab="Martingale residulas")

plot(dev_wb, xlab= "Data point", ylab="Deviance residulas",ylim=c(-4,4))
abline(h=3,col='red',lwd=2)
abline(h=-3,col='red',lwd=2)



df_rmi_wb = data.frame(
  obs = 1:length(rmi_wb),
  rmi = rmi_wb
)

ggplot(df_rmi_wb, aes(x = obs, y = rmi_wb)) +
  geom_point(size = 2, color = "black") +
  labs(
    title = "",
    x = "\n Index",
    y = "Martingale residual \n"
  ) +
  theme_minimal()


# Criação do dataframe correto
df_dev_wb = data.frame(
  obs = 1:length(dev_wb),
  rmi = dev_wb
)

# Gráfico com identificação dos pontos críticos
ggplot(df_dev_wb, aes(x = obs, y = rmi)) +
  geom_point(aes(color = (rmi < -3 | rmi > 3)), size = 2, show.legend = FALSE) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  geom_text_repel(
    data = subset(df_dev_wb, rmi < -3 | rmi > 3),
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



dev.off()

























