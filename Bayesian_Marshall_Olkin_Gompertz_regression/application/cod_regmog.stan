
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
  // Para computar métricas de comparação
  generated quantities {
    vector[N] log_lik;
    for (i in 1:N) {
     log_lik[i] = delta[i] * log(beta[i] * lambda * exp((beta[i] - beta[i] * exp(alpha[i] * time[i])) / alpha[i] + alpha[i] * time[i])) -
              delta[i] * log((1 - ((1-lambda) * exp((beta[i] - beta[i] * exp(alpha[i] * time[i])) / alpha[i])))^2) +
              (1 - delta[i]) * (log(lambda * exp(-(beta[i] / alpha[i]) * (exp(alpha[i] * time[i]) - 1)))) -
              (1 - delta[i]) * (log(1 - ((1 - lambda) * exp((-beta[i] / alpha[i]) * (exp(alpha[i] * time[i]) - 1)))));
  }
}

