setwd("~/Desktop/PhD/01_Joint_Model_Defective_Frailty/sim_MC_frailtydgomp/cen1")
library(dplyr)
#library(tidyverse)
library(patchwork) 
library(ggplot2)
library(tidyr)
library(tibble)
library(scales)


cen_1_jmres_100 = read.csv2("cen_1_jmres_rep_1000_n_100.csv", row.names = 1)
cen_1_jmres_300 = read.csv2("cen_1_jmres_rep_1000_n_300.csv", row.names = 1)
cen_1_jmres_500 = read.csv2("cen_1_jmres_rep_1000_n_500.csv", row.names = 1)
cen_1_jmres_1000 = read.csv2("cen_1_jmres_rep_1000_n_1000.csv", row.names = 1)

cen_1_0jmres_100 = read.csv2("cen_1_0jmres_rep_1000_n_100.csv", row.names = 1)
cen_1_0jmres_300 = read.csv2("cen_1_0jmres_rep_1000_n_300.csv", row.names = 1)
cen_1_0jmres_500 = read.csv2("cen_1_0jmres_rep_1000_n_500.csv", row.names = 1)
cen_1_0jmres_1000 = read.csv2("cen_1_0jmres_rep_1000_n_1000.csv", row.names = 1)

head(cen_1_jmres_100)
tail(cen_1_jmres_100)

head(cen_1_jmres_300)
tail(cen_1_jmres_300)

head(cen_1_jmres_500)
tail(cen_1_jmres_500)

head(cen_1_jmres_1000)
tail(cen_1_jmres_1000)


colMeans(cen_1_jmres_100)
colMeans(cen_1_jmres_300)
colMeans(cen_1_jmres_500)
colMeans(cen_1_jmres_1000)

parameters = c("alpha","gamma0","gamma1", "gamma2","psi1", "sd.frailty0",
               "beta0","beta1","beta2",
               "sd.e10", "sd.int","sd.slo","rho")

parameters0 = c("alpha","gamma0","gamma1", "gamma2","psi1",
               "beta0","beta1","beta2",
               "sd.e10", "sd.int","sd.slo","rho")

text_measures = c("mean", "sd", "bias", "cp")

names_parameters = paste0(
                    rep(parameters, times = length(text_measures)),
                    "_",
                     rep(text_measures, each = length(parameters)))

names_parameters0 = paste0( rep(parameters0, times = length(text_measures)), "_", rep(text_measures, each = length(parameters0)))

summary_mc = cbind(colMeans(cen_1_jmres_100),
              colMeans(cen_1_jmres_300),
              colMeans(cen_1_jmres_500),
              colMeans(cen_1_jmres_1000))

rownames(summary_mc) = names_parameters 
colnames(summary_mc) = c(100, 300, 500, 1000)

summary_mc0 = cbind(colMeans(cen_1_0jmres_100),
                   colMeans(cen_1_0jmres_300),
                   colMeans(cen_1_0jmres_500),
                   colMeans(cen_1_0jmres_1000))

rownames(summary_mc0) = names_parameters0 
colnames(summary_mc0) = c(100, 300, 500, 1000)

summary_mc0[25:48,]
summary_mc[27:52,]


# Carregue o pacote tidyverse (que inclui dplyr, tidyr e ggplot2)
library(tidyverse)

# 1. Transformar os dados do formato largo (wide) para longo (long/pivot)
tabela_pivot <- summary_mc[27:52, ] %>%
  as.data.frame() %>%
  # Transforma os nomes das linhas (alpha_bias, etc.) em uma coluna real
  rownames_to_column(var = "variavel") %>%
  # Pivota todas as colunas de amostra (100, 300, 500, 1000) para o formato longo
  pivot_longer(
    cols = -variavel, 
    names_to = "tamanho_amostra", 
    values_to = "valor"
  ) %>%
  # Separa "variavel" em duas colunas: "parametro" (ex: alpha) e "metrica" (ex: bias ou cp)
  # O sep = "_" funciona perfeitamente para o seu padrão de nomes
  separate(variavel, into = c("parametro", "metrica"), sep = "_", extra = "merge") %>%
  # Garante que o tamanho da amostra seja tratado como número (para o eixo X contínuo)
  mutate(tamanho_amostra = as.numeric(tamanho_amostra))

# 2. Visualizar a tabela transformada (opcional)
head(tabela_pivot)

# 3. Criar o gráfico com ggplot2
ggplot(tabela_pivot, aes(x = tamanho_amostra, y = valor, color = parametro)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  # Divide o gráfico em dois painéis: um para "bias" e outro para "cp"
  facet_wrap(~ metrica, scales = "free_y", ncol = 1) +
  # Linha de referência no Y dependendo da métrica (0 para bias, 0.95 para cp)
  geom_hline(data = data.frame(metrica = c("bias", "cp"), ref = c(0, 0.95)),
             aes(yintercept = ref), linetype = "dashed", color = "black", alpha = 0.5) +
  labs(
    title = "Comportamento de Bias e CP por Tamanho de Amostra",
    x = "Tamanho da Amostra (N)",
    y = "Valor",
    color = "Parâmetro"
  ) +
  theme_minimal() +
  theme(legend.position = "right")




# Instale o patchwork se ainda não tiver: install.packages("patchwork")


# 1. Preparação e transformação dos dados
tabela_pivot <- summary_mc[27:52, ] %>%
  as.data.frame() %>%
  rownames_to_column(var = "variavel") %>%
  pivot_longer(
    cols = -variavel, 
    names_to = "tamanho_amostra", 
    values_to = "valor"
  ) %>%
  separate(variavel, into = c("parametro", "metrica"), sep = "_", extra = "merge") %>%
  mutate(tamanho_amostra = as.numeric(tamanho_amostra))

# 2. Mapear os parâmetros para expressões matemáticas (letras gregas)
# Isso permite que o ggplot renderize alfa, beta_0, etc., igual à imagem.
tabela_pivot <- tabela_pivot %>%
  mutate(
    param_label = case_when(
      parametro == "alpha" ~ "alpha",
      parametro == "gamma0" ~ "gamma[0]",
      parametro == "gamma1" ~ "gamma[01]",
      parametro == "gamma2" ~ "gamma[11]",
      parametro == "psi1" ~ "psi[1]",
      parametro == "sd.frailty0" ~ "omega",
      parametro == "beta0" ~ "beta[0]",
      parametro == "beta1" ~ "beta[1]",
      parametro == "beta2" ~ "beta[2]",
      parametro == "sd.e10" ~ "sigma[epsilon]",
      parametro == "sd.int" ~ "sigma[b[0]]",
      parametro == "sd.slo" ~ "sigma[b[1]]",
      parametro == "rho" ~ "rho",
      TRUE ~ parametro
    )
  )


# Define a ordem correta para não ficar em ordem alfabética
tabela_pivot$param_label <- factor(tabela_pivot$param_label, levels = unique(tabela_pivot$param_label))

# 3. Criar um tema base limpo (igual ao da imagem)
tema_artigo <- theme_bw() +
  theme(
    strip.background = element_blank(), # Remove a caixa cinza do título
    strip.text = element_text(color = "black", size = 11), # Texto do parâmetro
    panel.grid.minor = element_blank(), # Deixa apenas as grades principais
    plot.caption = element_text(hjust = 0.5, face = "bold", size = 12), # Centraliza (a) e (b)
    axis.title = element_text(size = 10)
  )

# 4. Gráfico da Esquerda: BIAS (a)
p_bias <- tabela_pivot %>%
  filter(metrica == "bias") %>%
  ggplot(aes(x = tamanho_amostra, y = valor)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") + # Linha tracejada no 0
  geom_line(color = "#377eb8") + 
  geom_point(color = "#377eb8", size = 1.5) +
  # labeller = label_parsed é o que transforma os textos em letras gregas
  facet_wrap(~ param_label, ncol = 1, scales = "free_y", labeller = label_parsed) +
  scale_x_continuous(breaks = c(100, 300, 500, 1000)) +
  labs(
    x = "Sample size (n)",
    y = NULL, # Sem texto no eixo Y para o Bias, igual à imagem
  ) +
  tema_artigo

# 5. Gráfico da Direita: COVERAGE PROBABILITY (b)
p_cp <- tabela_pivot %>%
  filter(metrica == "cp") %>%
  ggplot(aes(x = tamanho_amostra, y = valor)) +
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "gray50") + # Linha pontilhada no 0.95
  geom_line(color = "#377eb8") + 
  geom_point(color = "#377eb8", size = 1.5) +
  facet_wrap(~ param_label, ncol = 1, scales = "free_y", labeller = label_parsed) +
  scale_x_continuous(breaks = c(100, 300, 500, 1000)) +
  labs(
    x = "Sample size (n)",
    y = "Coverage Probability",
  ) +
  tema_artigo

# 6. Juntar os gráficos lado a lado usando o patchwork

p_bias
p_cp

grafico_final <- p_bias | p_cp

# Exibir o gráfico
grafico_final





############

# 1. Preparar o primeiro objeto (summary_mc)
df1 <- summary_mc[27:52, ] %>%
  as.data.frame() %>%
  rownames_to_column(var = "variavel") %>%
  pivot_longer(cols = -variavel, names_to = "tamanho_amostra", values_to = "valor") %>%
  separate(variavel, into = c("parametro", "metrica"), sep = "_", extra = "merge") %>%
  mutate(
    tamanho_amostra = as.numeric(tamanho_amostra),
    modelo = "Frailty" # Identificador do modelo
  )

# 2. Preparar o segundo objeto (summary_mc0)
df0 <- summary_mc0[25:48, ] %>%
  as.data.frame() %>%
  rownames_to_column(var = "variavel") %>%
  pivot_longer(cols = -variavel, names_to = "tamanho_amostra", values_to = "valor") %>%
  separate(variavel, into = c("parametro", "metrica"), sep = "_", extra = "merge") %>%
  mutate(
    tamanho_amostra = as.numeric(tamanho_amostra),
    modelo = "No Frailty" # Identificador do segundo modelo
  )

# 3. Unir os dois conjuntos de dados
tabela_completa <- bind_rows(df1, df0)

# 4. Mapear os parâmetros para expressões matemáticas (letras gregas)
tabela_completa <- tabela_completa %>%
  mutate(
    param_label = case_when(
      parametro == "alpha" ~ "alpha",
      parametro == "gamma0" ~ "gamma[0]",
      parametro == "gamma1" ~ "gamma[0*1]", 
      parametro == "gamma2" ~ "gamma[11]",
      parametro == "psi1" ~ "psi[1]",
      parametro == "sd.frailty0" ~ "varphi",
      parametro == "beta0" ~ "beta[0]",
      parametro == "beta1" ~ "beta[1]",
      parametro == "beta2" ~ "beta[2]",
      parametro == "sd.e10" ~ "sigma[epsilon]",
      parametro == "sd.int" ~ "sigma[b[0]]",
      parametro == "sd.slo" ~ "sigma[b[1]]",
      parametro == "rho" ~ "rho",
      TRUE ~ parametro
    )
  )

# Garantir que a ordem dos painéis seja mantida (usando a ordem do df1)
ordem_paineis <- unique(df1$parametro)
ordem_labels <- tabela_completa$param_label[match(ordem_paineis, tabela_completa$parametro)]
tabela_completa$param_label <- factor(tabela_completa$param_label, levels = ordem_labels)

# 5. Criar um tema base limpo
tema_artigo <- theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 11),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_blank(), # Remove o título da legenda para ficar mais limpo
    legend.position = "bottom",
    legend.text = element_text(size = 14), 
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14), 
    
    axis.text.x = element_text(size = 10, color = "black"), 
    axis.text.y = element_text(size = 10, color = "black"),
  )

# theme(
#   legend.position = "top",
#   legend.title = element_text(size = 22),
#   legend.text  = element_text(size = 22),
#   strip.text   = element_text(size = 22),
#   panel.grid.minor = element_blank(),
#   
#   axis.title.x = element_text(size = 22, margin = margin(t = 15)), 
#   axis.title.y = element_text(size = 22,  margin = margin(r = 15)), 
#   
#   axis.text.x = element_text(size = 20, color = "black"), 
#   axis.text.y = element_text(size = 20, color = "black"),
#   
#   axis.ticks = element_line(linewidth = 1, color = "black")                   
# )

# 6. Gráfico da Esquerda: BIAS
p_bias <- tabela_completa %>%
  filter(metrica == "bias") %>%
  ggplot(aes(x = tamanho_amostra, y = valor, color = modelo, shape = modelo, linetype = modelo)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray3") + 
  geom_line(linewidth = 0.7) +  
  geom_point(size = 1.8) +
  # A alteração principal foi feita na linha abaixo (scales = "free")
  facet_wrap(~ param_label, ncol = 1, scales = "free", labeller = label_parsed) + 
  scale_x_continuous(breaks = c(100, 300, 500, 1000)) +
  scale_y_continuous(labels = function(x) sprintf("%.3f", x),
                     n.breaks = 4) +
  scale_color_manual(values = c("#377eb8", "darkorange")) + 
  labs(
    x = "\n Sample size (n)",
    y = "Bias \n"
  ) +
  tema_artigo

# 7. Gráfico da Direita: COVERAGE PROBABILITY (b)
p_cp <- tabela_completa %>%
  filter(metrica == "cp") %>%
  ggplot(aes(x = tamanho_amostra, y = valor, color = modelo, shape = modelo, linetype = modelo)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "gray3") + 
  geom_line(linewidth = 0.7) + 
  geom_point(size = 1.8) +
  # A alteração principal foi feita na linha abaixo (scales = "free")
  facet_wrap(~ param_label, ncol = 1, scales = "free", labeller = label_parsed) +
  scale_x_continuous(breaks = c(100, 300, 500, 1000)) +
  scale_y_continuous(labels = function(x) sprintf("%.3f", x),
                     limits = c(NA, 1), n.breaks = 4) +
  scale_color_manual(values = c("#377eb8", "darkorange")) +
  labs(
    x = "\n Sample size (n)",
    y = "Coverage Probability\n "
  ) +
  tema_artigo

p_bias
p_cp


# 8. Juntar os gráficos lado a lado e unificar a legenda no fundo
grafico_final <- (p_bias | p_cp) + 
  plot_layout(guides = "collect") & # Agrupa a legenda para não aparecer duplicada
  theme(legend.position = "bottom")

# Exibir o gráfico
grafico_final

p_bias
ggsave("p_bias.pdf", width = 9, height = 40, units = "cm")

p_cp
ggsave("p_cp.pdf", width = 9, height = 40, units = "cm")



# 1. Preparar o primeiro objeto (summary_mc)
df1 = summary_mc[27:52, ] %>%
  as.data.frame() %>%
  rownames_to_column(var = "variavel") %>%
  pivot_longer(cols = -variavel, names_to = "tamanho_amostra", values_to = "valor") %>%
  separate(variavel, into = c("parametro", "metrica"), sep = "_", extra = "merge") %>%
  mutate(
    tamanho_amostra = as.numeric(tamanho_amostra),
    modelo = "Frailty" # Identificador do modelo
  )

# 2. Preparar o segundo objeto (summary_mc0)
df0 = summary_mc0[25:48, ] %>%
  as.data.frame() %>%
  rownames_to_column(var = "variavel") %>%
  pivot_longer(cols = -variavel, names_to = "tamanho_amostra", values_to = "valor") %>%
  separate(variavel, into = c("parametro", "metrica"), sep = "_", extra = "merge") %>%
  mutate(
    tamanho_amostra = as.numeric(tamanho_amostra),
    modelo = "No Frailty" # Identificador do segundo modelo
  )

# 3. Unir os dois conjuntos de dados
tabela_completa <- bind_rows(df1, df0)

# 4. Mapear os parâmetros para expressões matemáticas (letras gregas)
tabela_completa <- tabela_completa %>%
  mutate(
    param_label = case_when(
      parametro == "alpha" ~ "alpha",
      parametro == "gamma0" ~ "gamma[0]",
      parametro == parametro == "gamma1" ~ "gamma[01]",,
      parametro == "gamma2" ~ "gamma[2]",
      parametro == "psi1" ~ "psi[1]",
      parametro == "sd.frailty0" ~ "varphi",
      parametro == "beta0" ~ "beta[0]",
      parametro == "beta1" ~ "beta[1]",
      parametro == "beta2" ~ "beta[2]",
      parametro == "sd.e10" ~ "sigma[varepsilon]",
      parametro == "sd.int" ~ "sigma[b[0]]",
      parametro == "sd.slo" ~ "sigma[b[1]]",
      parametro == "rho" ~ "rho",
      TRUE ~ parametro
    )
  )



#########

cen_1_oneminus_100 = read.csv2("cen_1_jmres_oneminus_crep_1000_n_100.csv", row.names = 1)
cen_1_oneminus_300 = read.csv2("cen_1_jmres_oneminus_crep_1000_n_300.csv", row.names = 1)
cen_1_oneminus_500 = read.csv2("cen_1_jmres_oneminus_crep_1000_n_500.csv", row.names = 1)
cen_1_oneminus_1000 = read.csv2("cen_1_jmres_oneminus_crep_1000_n_1000.csv", row.names = 1)

cen_1_0oneminus_100 = read.csv2("cen_1_0jmres_rep_1000_n_100.csv", row.names = 1)
cen_1_0oneminus_300 = read.csv2("cen_1_0jmres_rep_1000_n_300.csv", row.names = 1)
cen_1_0oneminus_500 = read.csv2("cen_1_0jmres_rep_1000_n_500.csv", row.names = 1)
cen_1_0oneminus_1000 = read.csv2("cen_1_0jmres_rep_1000_n_1000.csv", row.names = 1)


cen_1_oneminus_300=cen_1_oneminus_300[abs(cen_1_oneminus_300$bias.10) < 50,]
cen_1_oneminus_500=cen_1_oneminus_500[abs(cen_1_oneminus_500$bias.10) < 50,]

summary_mc_oneminus = cbind(colMeans(cen_1_oneminus_100),
                             colMeans(cen_1_oneminus_300),
                             colMeans(cen_1_oneminus_500),
                             colMeans(cen_1_oneminus_1000))


summary_mc_oneminus

colnames(summary_mc_oneminus) = c(100, 300, 500, 1000)


tc_vals = c(0.5, 1, 3, 5, 10, 15, 20, 25, 30, 35, 40)


names_parameters_mc_oneminus = paste0(
                                rep(tc_vals, times = length(text_measures)),
                                "_",
                                rep(text_measures, each = length(tc_vals)))


rownames(summary_mc_oneminus) = names_parameters_mc_oneminus


summary_mc_oneminus[23:33,]






# Instale os pacotes caso não os tenha: install.packages(c("ggplot2", "tidyr", "dplyr", "scales", "tibble"))

# 1. Função que esconde o texto abaixo de 1 e formata sem decimais os de 1 para cima
formatador_x_oculto <- function(x) {
  ifelse(is.na(x), NA, 
         ifelse(x >= 1, 
                sprintf("%.0f", x),  # Mostra o número (sem decimais) se for >= 1
                ""))                 # Retorna vazio (esconde o texto) se for < 1
}

# 2. Seu tema (mantido igual)
tema_artigo2 = theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 18),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_blank(), 
    legend.position = "bottom",
    legend.text = element_text(size = 22), 
    axis.title.x = element_text(size = 22), 
    axis.title.y = element_text(size = 22), 
    axis.text.x = element_text(size = 20, color = "black"), 
    axis.text.y = element_text(size = 20, color = "black")
  )

# 3. Gráfico
ggplot(dados_df, aes(x = tempo, y = vies, color = tamanho_amostra, group = tamanho_amostra)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 2.5) +
  scale_y_continuous(trans = pseudo_log_trans(base = 15)) +
  scale_x_continuous(
    breaks = tc_vals,              # Mantém todas as marcações (risquinhos)
    labels = formatador_x_oculto   # Esconde o texto dos que são < 1
  ) +
  scale_color_discrete(labels = c("n=100", "n=300", "n=500", "n=1000")) +
  labs(
    x = "\n Time (Years)",
    y = "Bias \n",
    color = "Sample size (n) \n"
  ) +
  tema_artigo2



































tabela_pivot_mc_oneminus = summary_mc_oneminus[23:33,] %>%
              as.data.frame() %>%
              rownames_to_column(var = "variavel") %>%
              pivot_longer(
                cols = -variavel, 
                names_to = "tamanho_amostra", 
                values_to = "valor"
              ) %>%
              separate(variavel, into = c("parametro", "metrica"), sep = "_", extra = "merge") %>%
              mutate(tamanho_amostra = as.numeric(tamanho_amostra))


tabela_completa <- tabela_completa %>%
  mutate(
    param_label = case_when(
      parametro == "alpha" ~ "alpha",
      parametro == "gamma0" ~ "gamma[0]",
      parametro == "gamma1" ~ "gamma[1]",
      parametro == "gamma2" ~ "gamma[2]",
      parametro == "psi1" ~ "psi[1]",
      parametro == "sd.frailty0" ~ "omega",
      parametro == "beta0" ~ "beta[0]",
      parametro == "beta1" ~ "beta[1]",
      parametro == "beta2" ~ "beta[2]",
      parametro == "sd.e10" ~ "sigma[epsilon]",
      parametro == "sd.int" ~ "sigma[b[0]]",
      parametro == "sd.slo" ~ "sigma[b[1]]",
      parametro == "rho" ~ "rho",
      TRUE ~ parametro
    )
  )


tabela_pivot_mc_oneminus










summary_oneminus = colMeans(cen_1_oneminus_100)

# 1. Preparar o primeiro objeto (summary_mc)
df1 = summary_mc[27:52, ] %>%
  as.data.frame() %>%
  rownames_to_column(var = "variavel") %>%
  pivot_longer(cols = -variavel, names_to = "tamanho_amostra", values_to = "valor") %>%
  separate(variavel, into = c("parametro", "metrica"), sep = "_", extra = "merge") %>%
  mutate(
    tamanho_amostra = as.numeric(tamanho_amostra),
    modelo = "Frailty" # Identificador do modelo
  )

# 2. Preparar o segundo objeto (summary_mc0)
df0 = summary_mc0[25:48, ] %>%
  as.data.frame() %>%
  rownames_to_column(var = "variavel") %>%
  pivot_longer(cols = -variavel, names_to = "tamanho_amostra", values_to = "valor") %>%
  separate(variavel, into = c("parametro", "metrica"), sep = "_", extra = "merge") %>%
  mutate(
    tamanho_amostra = as.numeric(tamanho_amostra),
    modelo = "No Frailty" # Identificador do segundo modelo
  )

# 3. Unir os dois conjuntos de dados
tabela_completa <- bind_rows(df1, df0)

# 4. Mapear os parâmetros para expressões matemáticas (letras gregas)
tabela_completa <- tabela_completa %>%
  mutate(
    param_label = case_when(
      parametro == "alpha" ~ "alpha",
      parametro == "gamma0" ~ "gamma[0]",
      parametro == "gamma1" ~ "gamma[1]",
      parametro == "gamma2" ~ "gamma[2]",
      parametro == "psi1" ~ "psi[1]",
      parametro == "sd.frailty0" ~ "omega",
      parametro == "beta0" ~ "beta[0]",
      parametro == "beta1" ~ "beta[1]",
      parametro == "beta2" ~ "beta[2]",
      parametro == "sd.e10" ~ "sigma[epsilon]",
      parametro == "sd.int" ~ "sigma[b[0]]",
      parametro == "sd.slo" ~ "sigma[b[1]]",
      parametro == "rho" ~ "rho",
      TRUE ~ parametro
    )
  )


























































head(cen_1_oneminus_100)
colMeans(cen_1_oneminus_100)


plot(1:dim(cen_1_oneminus_100)[1], cen_1_oneminus_100$bias, type='l')
plot(1:dim(cen_1_oneminus_100)[1], cen_1_0oneminus_100$bias, type='l')


colMeans(cen_1_oneminus_100)

plot(1:dim(cen_1_oneminus_100)[1], cen_1_oneminus_100$bias.1, type='l')
plot(1:dim(cen_1_oneminus_100)[1], cen_1_oneminus_100$bias.2, type='l')
plot(1:dim(cen_1_oneminus_100)[1], cen_1_oneminus_100$bias.3, type='l')
plot(1:dim(cen_1_oneminus_100)[1], cen_1_oneminus_100$bias.4, type='l')
plot(1:dim(cen_1_oneminus_100)[1], cen_1_oneminus_100$bias.5, type='l')
plot(1:dim(cen_1_oneminus_100)[1], cen_1_oneminus_100$bias.6, type='l')
plot(1:dim(cen_1_oneminus_100)[1], cen_1_oneminus_100$bias.7, type='l')
plot(1:dim(cen_1_oneminus_100)[1], cen_1_oneminus_100$bias.8, type='l')
plot(1:dim(cen_1_oneminus_100)[1], cen_1_oneminus_100$bias.9, type='l')
plot(1:dim(cen_1_oneminus_100)[1], cen_1_oneminus_100$bias.10, type='l')




tabela_completa %>%
  filter(metrica == "bias") %>%
  ggplot(aes(x = tamanho_amostra, y = valor, color = modelo, shape = modelo, linetype = modelo)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray3") + 
  geom_line(linewidth = 0.7) +  
  geom_point(size = 1.8) +
  facet_wrap(~ param_label, ncol = 1, scales = "free_y", labeller = label_parsed) +
  scale_x_continuous(breaks = c(100, 300, 500, 1000)) +
  scale_y_continuous(labels = function(x) sprintf("%.3f", x),
                     n.breaks = 4) +
  scale_color_manual(values = c("#377eb8", "darkorange")) + # Azul para um, vermelho para o outro
  labs(
    x = "\n Sample size (n)",
    y = "Bias \n"
  ) +
  tema_artigo










