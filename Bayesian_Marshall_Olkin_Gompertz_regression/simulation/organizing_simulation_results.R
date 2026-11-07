
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggh4x) 


setwd("~/Code-for-papers/Bayesian_Marshall_Olkin_Gompertz_regression/simulation/simulation_results/scen1")

scen1_50 = read.csv2("1simres_50.csv")[,-1]
scen1_100 = read.csv2("1simres_100.csv")[,-1]
scen1_300 = read.csv2("1simres_300.csv")[,-1]
scen1_500 = read.csv2("1simres_500.csv")[,-1]
scen1_1000 = read.csv2("1simres_1000.csv")[,-1]

# head(scen1_50)
# dim(scen1_50)
# 
# head(scen1_1000)
# dim(scen1_1000)

setwd("~/Code-for-papers/Bayesian_Marshall_Olkin_Gompertz_regression/simulation/simulation_results/scen2")

scen2_50 = read.csv2("2simres_50.csv")[,-1]
scen2_100 = read.csv2("2simres_100.csv")[,-1]
scen2_300 = read.csv2("2simres_300.csv")[,-1]
scen2_500 = read.csv2("2simres_500.csv")[,-1]
scen2_1000 = read.csv2("2simres_1000.csv")[,-1]


sumres1 = rbind(colMeans(scen1_50)[11:20],
                colMeans(scen1_100)[11:20],
                colMeans(scen1_300)[11:20],
                colMeans(scen1_500)[11:20],
                colMeans(scen1_1000)[11:20])


sumres2 = rbind(colMeans(scen2_50)[11:20],
                colMeans(scen2_100)[11:20],
                colMeans(scen2_300)[11:20],
                colMeans(scen2_500)[11:20],
                colMeans(scen2_1000)[11:20])

colnames(sumres1) = colnames(sumres2) = c("alpha.cp", "Upsilon.cp","mu0.cp", "mu1.cp", "mu2.cp",
                                          "alpha.bias", "Upsilon.bias","mu0.bias", "mu1.bias", "mu2.bias")


rownames(sumres2) = rownames(sumres1) = c("n50", "n100", "n300", "n500", "n1000")



df1 = as.data.frame(sumres1)
df1$cenario = "Scenario 1"

df2 = as.data.frame(sumres2)
df1$n = df2$n = as.numeric(gsub("n", "", rownames(df2)))
df2$cenario = "Scenario 2"


df_all = bind_rows(df1, df2)
df_long = df_all %>%
  pivot_longer(
    cols = -c(n, cenario),
    names_to = c("parameter", "metric"),
    names_sep = "\\.",
    values_to = "value"
  )


df_long$parameter = factor(
  df_long$parameter, 
  levels = c("alpha", "Upsilon", "mu0", "mu1", "mu2"),
  labels = c(
    "alpha", "Upsilon", "psi[mu[0]]",  "psi[mu[1]]", "psi[mu[2]]"))


df_cp = df_long %>% filter(metric == "cp")
df_bias = df_long %>% filter(metric == "bias")

colors = c("Scenario 1" = "green4", "Scenario 2" = "steelblue")


# =====================================================================
# GRAPH 1: Coverage Probability
# =====================================================================
plot_cp = ggplot(df_cp, aes(x = n, y = value, color = cenario, group = cenario)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black", linewidth = 1.2) +
  geom_line(linewidth = 1.4) +
  geom_point(size = 2) +
  scale_color_manual(values = colors) +
  scale_x_continuous(breaks = unique(df_all$n)) +
  scale_y_continuous(limits = c(0.88, 1)) +
  facet_wrap(~ parameter, ncol = 1, scales = "free_x", labeller = label_parsed) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 22, face = "bold"),
    legend.position = "top",
    legend.text = element_text(size = 22),
    legend.title = element_blank(),
    
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 22, margin = margin(t = 15)), 
    axis.title.y = element_text(size = 22,  margin = margin(r = 15)), 
    
    axis.text.x = element_text(size = 13.4, color = "black"), 
    axis.text.y = element_text(size = 15, color = "black"),
    
    axis.ticks = element_line(linewidth = 1, color = "black")
  ) +
  labs(x = "Sample size (n)", y = "Coverage Probability")

print(plot_cp)


# =====================================================================
# GRAPH 2: Bias 
# =====================================================================


limites_bias = data.frame(
  parameter = factor(c("alpha", "Upsilon", "psi[mu[1]]", "psi[mu[2]]"), 
                     levels = c("alpha", "Upsilon", "psi[mu[0]]", "psi[mu[1]]", "psi[mu[2]]")),
  value = 0.05,
  n = 50)


plot_bias = ggplot(df_bias, aes(x = n, y = value, color = cenario, group = cenario)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1.2) +
  geom_blank(data = limites_bias, aes(x = n, y = value), inherit.aes = FALSE) +
  geom_line(linewidth = 1.4) +
  geom_point(size = 2) +
  scale_color_manual(values = colors) +
  scale_x_continuous(breaks = unique(df_all$n)) +
  scale_y_continuous(
    labels = function(x) sprintf("%.3f", x)
  ) +
  facet_wrap(~ parameter, ncol = 1, scales = "free", labeller = label_parsed) +
  theme_bw(base_size = 16) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 22, face = "bold"),
    legend.position = "top",
    legend.text = element_text(size = 22),
    legend.title = element_blank(),
    
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 22, margin = margin(t = 15)), 
    axis.title.y = element_text(size = 22,  margin = margin(r = 15)), 
    
    axis.text.x = element_text(size = 13.4, color = "black"), 
    axis.text.y = element_text(size = 15, color = "black"),
    
    axis.ticks = element_line(linewidth = 1, color = "black")
  ) +
  labs(x = "Sample size (n)", y = "Bias")


print(plot_bias)








