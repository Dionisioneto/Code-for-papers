###### Organizing the data from priors sensibility study ###### 

library(tidyr)
library(dplyr)
library(ggh4x) 

## Two s and three kinds of prior configuration

# setwd("~/Code-for-papers/Bayesian_Marshall_Olkin_Gompertz_regression/simulation/simulation_results/prior_sensibility/11")
# 
# scen11_50 = read.csv2("1simres_50.csv")[,-1]
# scen11_100 = read.csv2("1simres_100.csv")[,-1]
# scen11_300 = read.csv2("1simres_300.csv")[,-1]
# scen11_500 = read.csv2("1simres_500.csv")[,-1]
# scen11_1000 = read.csv2("1simres_1000.csv")[,-1]
# 
# setwd("~/Code-for-papers/Bayesian_Marshall_Olkin_Gompertz_regression/simulation/simulation_results/prior_sensibility/13")


# =====================================================================
# Second scenario
# =====================================================================

setwd("~/Code-for-papers/Bayesian_Marshall_Olkin_Gompertz_regression/simulation/simulation_results/prior_sensibility/21")
scen21_50 = read.csv2("1simres_50.csv")[,-1]
scen21_100 = read.csv2("1simres_100.csv")[,-1]
scen21_300 = read.csv2("1simres_300.csv")[,-1]
scen21_500 = read.csv2("1simres_500.csv")[,-1]
scen21_1000 = read.csv2("1simres_1000.csv")[,-1]

sumres21 = rbind(colMeans(scen21_50)[11:20],
                colMeans(scen21_100)[11:20],
                colMeans(scen21_300)[11:20],
                colMeans(scen21_500)[11:20],
                colMeans(scen21_1000)[11:20])


setwd("~/Code-for-papers/Bayesian_Marshall_Olkin_Gompertz_regression/simulation/simulation_results/prior_sensibility/22")
scen22_50 = read.csv2("1simres_50.csv")[,-1]
scen22_100 = read.csv2("1simres_100.csv")[,-1]
scen22_300 = read.csv2("1simres_300.csv")[,-1]
scen22_500 = read.csv2("1simres_500.csv")[,-1]
scen22_1000 = read.csv2("1simres_1000.csv")[,-1]

sumres22 = rbind(colMeans(scen22_50)[11:20],
                 colMeans(scen22_100)[11:20],
                 colMeans(scen22_300)[11:20],
                 colMeans(scen22_500)[11:20],
                 colMeans(scen22_1000)[11:20])



setwd("~/Code-for-papers/Bayesian_Marshall_Olkin_Gompertz_regression/simulation/simulation_results/prior_sensibility/23")
scen23_50 = read.csv2("1simres_50.csv")[,-1]
scen23_100 = read.csv2("1simres_100.csv")[,-1]
scen23_300 = read.csv2("1simres_300.csv")[,-1]
scen23_500 = read.csv2("1simres_500.csv")[,-1]
scen23_1000 = read.csv2("1simres_1000.csv")[,-1]

sumres23 = rbind(colMeans(scen23_50)[11:20],
                 colMeans(scen23_100)[11:20],
                 colMeans(scen23_300)[11:20],
                 colMeans(scen23_500)[11:20],
                 colMeans(scen23_1000)[11:20])



rownames(sumres21) =  rownames(sumres22) = rownames(sumres23) = c("n50", "n100", "n300", "n500", "n1000")

colnames(sumres21) = colnames(sumres22) =  colnames(sumres23) = c("alpha.cp", "Upsilon.cp","mu0.cp", "mu1.cp", "mu2.cp",
                                                                  "alpha.bias", "Upsilon.bias","mu0.bias", "mu1.bias", "mu2.bias")



# Organizing dataset
df21 = as.data.frame(sumres21)
df21$scenario = "Prior 1"

df22 = as.data.frame(sumres22)
df22$scenario = "Prior 2"

df23 = as.data.frame(sumres23)
df23$scenario = "Prior 3"

df21$n = df22$n = df23$n  = as.numeric(gsub("n", "", rownames(df21)))



df_all_cen2 = bind_rows(df21, df22, df23)
df_long_cen2 = df_all_cen2 %>%
  pivot_longer(
    cols = -c(n, scenario),
    names_to = c("parameter", "metric"),
    names_sep = "\\.",
    values_to = "value"
  )


df_long_cen2$parameter = factor(
  df_long_cen2$parameter, 
  levels = c("alpha", "Upsilon", "mu0", "mu1", "mu2"),
  labels = c(
    "alpha", "Upsilon", "psi[mu[0]]", "psi[mu[1]]","psi[mu[2]]"    
  )
)


df_cp_cen2 = df_long_cen2 %>% filter(metric == "cp")
df_bias_cen2 = df_long_cen2 %>% filter(metric == "bias")

colors = c("Prior 1" = "green4", "Prior 2" = "steelblue", "Prior 3" = "purple")



# =====================================================================
# Graph 1: Coverage Probability
# =====================================================================

plot_cp22 = ggplot(df_cp_cen2, aes(x = n, y = value, color = scenario, group = scenario)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black", linewidth = 1.2) +
  geom_line(linewidth = 1.4) +
  geom_point(size = 2) +
  scale_color_manual(values = colors) +
  scale_x_continuous(breaks = unique(df_cp_cen2$n)) +
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
  )+
  labs(x = "Sample size (n)", y = "Coverage Probability")

print(plot_cp22)


# =====================================================================
# Graph 2: Bias 
# =====================================================================

limits_bias = data.frame(
  parameter = factor(c("alpha", "Upsilon", "psi[mu[1]]", "psi[mu[2]]"), 
                     levels = c("alpha", "Upsilon", "psi[mu[0]]", "psi[mu[1]]", "psi[mu[2]]")),
  value = 0.05,
  n = 50 
)


plot_bias22 = ggplot(df_bias_cen2, aes(x = n, y = value, color = scenario, group = scenario)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1.2) +
  geom_blank(data = limits_bias, aes(x = n, y = value), inherit.aes = FALSE) +
  geom_line(linewidth = 1.4) +
  geom_point(size = 2) +
  scale_color_manual(values = colors) +
  scale_x_continuous(breaks = unique(df_bias_cen2$n)) +
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
  labs(x = "Sample size (n)", y = "Bias", color = "Cenário:")


print(plot_bias22)


# =====================================================================
# First scenario
# =====================================================================


setwd("~/Code-for-papers/Bayesian_Marshall_Olkin_Gompertz_regression/simulation/simulation_results/prior_sensibility/11")
scen11_50 = read.csv2("1simres_50.csv")[,-1]
scen11_100 = read.csv2("1simres_100.csv")[,-1]
scen11_300 = read.csv2("1simres_300.csv")[,-1]
scen11_500 = read.csv2("1simres_500.csv")[,-1]
scen11_1000 = read.csv2("1simres_1000.csv")[,-1]


sumres11 = rbind(colMeans(scen11_50)[11:20],
                 colMeans(scen11_100)[11:20],
                 colMeans(scen11_300)[11:20],
                 colMeans(scen11_500)[11:20],
                 colMeans(scen11_1000)[11:20])


colnames(sumres11) = c("alpha.cp", "Upsilon.cp","mu0.cp", "mu1.cp", "mu2.cp",
                       "alpha.bias", "Upsilon.bias","mu0.bias", "mu1.bias", "mu2.bias")

## Simulation study 1 is the same of prior sensibility study when we onsider the 11  of priors
setwd("~/Code-for-papers/Bayesian_Marshall_Olkin_Gompertz_regression/simulation/simulation_results/scen1")

scen12_50 = read.csv2("1simres_50.csv")[,-1]
scen12_100 = read.csv2("1simres_100.csv")[,-1]
scen12_300 = read.csv2("1simres_300.csv")[,-1]
scen12_500 = read.csv2("1simres_500.csv")[,-1]
scen12_1000 = read.csv2("1simres_1000.csv")[,-1]

sumres12 = rbind(colMeans(scen12_50)[11:20],
                 colMeans(scen12_100)[11:20],
                 colMeans(scen12_300)[11:20],
                 colMeans(scen12_500)[11:20],
                 colMeans(scen12_1000)[11:20])


colnames(sumres12) = c("alpha.cp", "Upsilon.cp","mu0.cp", "mu1.cp", "mu2.cp",
                       "alpha.bias", "Upsilon.bias","mu0.bias", "mu1.bias", "mu2.bias")

setwd("~/Code-for-papers/Bayesian_Marshall_Olkin_Gompertz_regression/simulation/simulation_results/prior_sensibility/13")

scen13_50 = read.csv2("1simres_50.csv")[,-1]
scen13_100 = read.csv2("1simres_100.csv")[,-1]
scen13_300 = read.csv2("1simres_300.csv")[,-1]

scen13_500_1 = read.csv2("11simres_500.csv")[,-1]
scen13_500_2 = read.csv2("21simres_500.csv")[,-1]
scen13_500 = rbind(scen13_500_1, scen13_500_2)

scen13_1000_1 = read.csv2("31simres_1000.csv")[,-1]
scen13_1000_2 = read.csv2("41simres_1000.csv")[,-1]
scen13_1000_3 = read.csv2("51simres_1000.csv")[,-1]
scen13_1000_4 = read.csv2("61simres_1000.csv")[,-1]
scen13_1000 = rbind(scen13_1000_1, scen13_1000_2,
                    scen13_1000_3, scen13_1000_4)


sumres13 = rbind(colMeans(scen13_50)[11:20],
                 colMeans(scen13_100)[11:20],
                 colMeans(scen13_300)[11:20],
                 colMeans(scen13_500)[11:20],
                 colMeans(scen13_1000)[11:20])

rownames(sumres11) = rownames(sumres12) = rownames(sumres13) = c("n50","n100","n300","n500","n1000")

colnames(sumres13) = c("alpha.cp", "Upsilon.cp","mu0.cp", "mu1.cp", "mu2.cp",
                       "alpha.bias", "Upsilon.bias","mu0.bias", "mu1.bias", "mu2.bias")



df11 = as.data.frame(sumres11)
df11$scenario = "Prior 1"

df12 = as.data.frame(sumres12)
df12$scenario = "Prior 2"

df13 = as.data.frame(sumres13)
df13$scenario = "Prior 3"

df11$n = df12$n = df13$n = as.numeric(gsub("n", "", rownames(df11)))


df_all_cen3 = bind_rows(df11, df12, df13)
df_long_cen3 = df_all_cen3 %>%
  pivot_longer(
    cols = -c(n, scenario),
    names_to = c("parameter", "metric"),
    names_sep = "\\.",
    values_to = "value"
  )

df_long_cen3$parameter = factor(
  df_long_cen3$parameter, 
  levels = c("alpha", "Upsilon", "mu0", "mu1", "mu2"),
  labels = c(
    "alpha",   
    "Upsilon", 
    "psi[mu[0]]",  
    "psi[mu[1]]",   
    "psi[mu[2]]"))


df_cp_cen3 = df_long_cen3 %>% filter(metric == "cp")
df_bias_cen3 = df_long_cen3 %>% filter(metric == "bias")
colors = c("Prior 1" = "green4", "Prior 2" = "steelblue", "Prior 3" = "purple")


# =====================================================================
# Graph 1: Coverage Probability
# =====================================================================
plot_cp11 = ggplot(df_cp_cen3, aes(x = n, y = value, color = scenario, group = scenario)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black", linewidth = 1.2) +
  geom_line(linewidth = 1.4) +
  geom_point(size = 2) +
  scale_color_manual(values = colors) +
  scale_x_continuous(breaks = unique(df_cp_cen3$n)) +
  scale_y_continuous(labels = function(x) sprintf("%.3f", x)) + # <-- Alteração feita aqui
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
  )+
  labs(x = "Sample size (n)", y = "Coverage Probability", color = "Cenário:")

print(plot_cp11)


# =====================================================================
# GRÁFICO 2: Bias (Com limite superior de 0.1 para todos exceto mu0)
# =====================================================================


plot_bias11 = ggplot(df_bias_cen3, aes(x = n, y = value, color = scenario, group = scenario)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1.2) +
  geom_line(linewidth = 1.4) +
  geom_point(size = 2) +
  scale_color_manual(values = colors) +
  scale_x_continuous(breaks = unique(df_bias_cen3$n)) +
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
  labs(x = "Sample size (n)", y = "Bias", color = "Cenário:")


print(plot_bias11)















