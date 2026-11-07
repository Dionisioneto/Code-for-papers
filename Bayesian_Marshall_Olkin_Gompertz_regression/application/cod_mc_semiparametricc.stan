
      data {
      int<lower=1> N;                        // Número de pacientes
      int<lower=1> p_cov;                    // Número de covariáveis para a probabilidade de cura (logit)
      int<lower=1> lat_cov;                  // Número de covariáveis para o grupo não curado (opcional, aqui omitido por simplicidade ou pode usar X)
      
      matrix[N, p_cov] Z;                    // Matriz de desenho para a fração de cura (logística)
      
      vector<lower=0>[N] t_obs;              // Tempos observados
      vector[N] delta;                       // Indicador (1 = falha, 0 = censura)
      
      // Estrutura para o baseline semiparamétrico (Piecewise) dos não-curados
      int<lower=1> J;                        // Número de intervalos
      array[N] int<lower=1, upper=J> t_idx;  // Índice do intervalo final
      matrix[N, J] D;                        // Matriz de exposição
    }
    
    parameters {
      vector[p_cov] gamma;                   // Coeficientes logísticos para a probabilidade de cura
      
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
      gamma ~ normal(0, 3);                  // Coeficientes da logística
      sigma_lambda ~ normal(0, 1);           // Suavização do risco basal
      log_lambda[1] ~ normal(-2, 3);
      
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

