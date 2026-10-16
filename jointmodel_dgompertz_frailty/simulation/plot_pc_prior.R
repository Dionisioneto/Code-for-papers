# Carregar pacotes necessários
library(ggplot2)
library(INLA) # Opcional, apenas se quiser usar inla.pc.dsd()
library(dplyr)

# 1. Definir os vetores de parâmetros para comparação
u_vals <- c(0.5, 1.0, 2.0)
alpha_vals <- c(0.01, 0.05, 0.1)

# Gerar a sequência de valores para sigma (eixo X)
sigma_vals <- seq(0, 5, length.out = 500)

# 2. Criar um grid com todas as combinações de u, alpha e sigma
df_grid <- expand.grid(
  u = u_vals, 
  alpha = alpha_vals, 
  sigma = sigma_vals
)

# 3. Calcular a taxa lambda e criar a coluna alpha_label
df <- df_grid %>%
  mutate(
    lambda = -log(alpha) / u,
    densidade = lambda * exp(-lambda * sigma),
    u_factor = as.factor(u),
    # Crie uma STRING (texto) com a sintaxe matemática ("==" cria o símbolo de igualdade)
    alpha_label = paste0("b[xi] == ", alpha) 
  )

# 4. Plotar o gráfico com os painéis
ggplot(df, aes(x = sigma, y = densidade, color = u_factor)) +
  geom_line(linewidth = 1.2) +
  # O xintercept precisa estar no aes() para acompanhar as diferentes cores de 'u'
  #geom_vline(aes(xintercept = u, color = u_factor), linetype = "dashed") +
  facet_wrap(~ alpha_label, labeller = label_parsed) +
  labs(
    x = expression(xi),
    y = expression(pi(xi)),
    color = expression(a[xi]) # No labs(), usar expression() é a forma correta
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(size = 12),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  ) +
  theme(
    legend.margin = margin(t = 5, r = 10, b = 5, l = 10),
    legend.title = element_text(size = 20),
    legend.text  = element_text(size = 22),
    strip.text   = element_text(size = 22),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 22, margin = margin(t = 15)), 
    axis.title.y = element_text(size = 22,  margin = margin(r = 15)), 
    
    axis.text.x = element_text(size = 20, color = "black"), 
    axis.text.y = element_text(size = 20, color = "black"),
    
    axis.ticks = element_line(linewidth = 1, color = "black")
  )





