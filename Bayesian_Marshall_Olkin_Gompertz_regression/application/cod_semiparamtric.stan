
    data {
      int<lower=1> N;                        // Número de pacientes
      int<lower=1> s;                        // Número de covariáveis (incluindo intercepto)
      matrix[N, s] X;                        // Matriz de planejamento
      
      vector<lower=0>[N] t_obs;              // Tempos observados (falha ou censura)
      vector[N]  delta;  // Indicador (1 = falha, 0 = censura)
      
      // Estrutura para o baseline semiparamétrico (Piecewise)
      int<lower=1> J;                        // Número total de intervalos de tempo
      array[N] int<lower=1, upper=J> t_idx;  // Índice do intervalo final em que a observação caiu
      matrix[N, J] D;                        // Matriz de exposição: D[i,j] = tempo que o indivíduo i passou no intervalo j
    }
    
    parameters {
      vector[s] beta;                        // Efeitos preditores sobre a fração de cura
      
      // Parâmetros para a função F(t) não-paramétrica
      vector[J] log_lambda;                  // Log-risco basal em cada intervalo
      real<lower=0> sigma_lambda;            // Desvio padrão para a priori de suavização (RW1)
    }
    
    transformed parameters {
      vector[N] theta;
      vector[N] F_t;
      vector[N] f_t;
      vector[J] lambda = exp(log_lambda);
      
      // Preditores lineares da taxa de cura
      theta = exp(X * beta); 
      
      // Cálculo exato de F(t) e f(t) para cada indivíduo com base na exposição
      for (i in 1:N) {
        real Lambda_t = 0;
        
        // Soma o risco cumulativo pelos intervalos em que o indivíduo esteve em risco
        for (j in 1:t_idx[i]) {
          Lambda_t += lambda[j] * D[i, j];
        }
        
        // F(t) e f(t) derivados da função de risco
        F_t[i] = 1 - exp(-Lambda_t);
        f_t[i] = lambda[t_idx[i]] * exp(-Lambda_t);
      }
    }
    
    model {
      // 1. Prioris para os coeficientes paramétricos
      beta ~ normal(0, 10);
      
      // 2. Priori Semiparamétrica de Suavização (Random Walk de ordem 1)
      sigma_lambda ~ normal(0, 5); // Priori half-normal devido ao limite inferior <lower=0>
      log_lambda[1] ~ normal(0, 5);
      for (j in 2:J) {
        log_lambda[j] ~ normal(log_lambda[j-1], sigma_lambda);
      }
      
      // 3. Verossimilhança Marginal
      for (i in 1:N) {
        if (delta[i] == 1) {
          // Falha: log(h(t)) + log(S(t)) => log(theta * f(t)) - theta * F(t)
          target += log(theta[i]) + log(f_t[i]) - theta[i] * F_t[i];
        } else {
          // Censura: log(S(t)) => - theta * F(t)
          target += - theta[i] * F_t[i];
        }
      }
    }
    
    generated quantities {
      vector[N] log_lik; // Armazenamento estrito para o WAIC
      
      for (i in 1:N) {
        if (delta[i] == 1) {
          log_lik[i] = log(theta[i]) + log(f_t[i]) - theta[i] * F_t[i];
        } else {
          log_lik[i] = - theta[i] * F_t[i];
        }
      }
    }

