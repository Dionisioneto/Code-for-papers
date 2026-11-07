library(rstan)
library(ggplot2)
library(survminer)
library(survival)
library(dplyr)

#load source data
#source("~/Desktop/PhD/papers/Cap4_MO_Gompertz/application/source_functions.R")

# Load data
test_data = read.csv2("test_data.csv")

# Create a 'Stage' categorical variable based on binary indicators s2 and s3
test_data$Stage = ifelse(test_data$s2 == 0 & test_data$s3 == 0, "I",
                         ifelse(test_data$s2 == 1 & test_data$s3 == 0, "II", 
                                ifelse(test_data$s2 == 0 & test_data$s3 == 1, "III", NA)))

# Convert to factor with a specific ordered level for plotting
test_data$Stage = factor(test_data$Stage, levels = c("I", "II", "III"))

# Fit the Kaplan-Meier model by Stage
kmfit = survfit(Surv(time, delta) ~ Stage, data = test_data)

# Create the customized survival plot
plot = ggsurvplot(kmfit, data = test_data,
                  conf.int = TRUE,       
                  xlab = "\n Time (Years)",
                  legend = "top",
                  ylab = "Survival probability \n ",
                  palette = c("steelblue", "yellow4", "red4"), 
                  censor = FALSE,
                  legend.title = "" 
) 

# Adjust plot theme and font sizes
pp  = plot$plot + theme(
  axis.text.x = element_text(size = 24),
  axis.text.y = element_text(size = 24),
  axis.title.x = element_text(size = 26),
  axis.title.y = element_text(size = 26),
  legend.title = element_blank(),
  legend.text = element_text(size = 24)
)

pp

# Define the Stan model as a string
cod_regmog_app = "
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
      vector[N] beta = exp(xb * bs);
    }

    model {
      alpha ~ normal(0,10);
      for (j in 1:s) {
        bs[j] ~ normal(0, 10);
      }
      lambda ~ gamma(0.1,0.1);

      // Custom log-likelihood for the survival model
      for (i in 1:N) {
        target += delta[i] * log(beta[i] * lambda * exp((beta[i] - beta[i] * exp(alpha * time[i])) / alpha + alpha * time[i])) -
                  delta[i] * log((1 - ((1-lambda) * exp((beta[i] - beta[i] * exp(alpha * time[i])) / alpha)))^2) +
                  (1 - delta[i]) * (log(lambda * exp(-(beta[i] / alpha) * (exp(alpha * time[i]) - 1)))) -
                  (1 - delta[i]) * (log(1 - ((1 - lambda) * exp((-beta[i] / alpha) * (exp(alpha * time[i]) - 1)))));
      }
    }
    
    generated quantities {
      vector[N] log_lik;
      // Log-likelihood for model comparison
      for (i in 1:N) {
         log_lik[i] = delta[i] * log(beta[i] * lambda * exp((beta[i] - beta[i] * exp(alpha * time[i])) / alpha + alpha * time[i])) -
                  delta[i] * log((1 - ((1-lambda) * exp((beta[i] - beta[i] * exp(alpha * time[i])) / alpha)))^2) +
                  (1 - delta[i]) * (log(lambda * exp(-(beta[i] / alpha) * (exp(alpha * time[i]) - 1)))) -
                  (1 - delta[i]) * (log(1 - ((1 - lambda) * exp((-beta[i] / alpha) * (exp(alpha * time[i]) - 1)))));
      }

      // Cure fraction per scenario (Stage I, II, III)
      real p01 = (lambda * exp(exp(bs[1]) / alpha)) / (1 - (1 - lambda) * exp(exp(bs[1]) / alpha));
      real p02 = (lambda * exp(exp(bs[1] + bs[2]) / alpha)) / (1 - (1 - lambda) * exp(exp(bs[1] + bs[2]) / alpha));
      real p03 = (lambda * exp(exp(bs[1] + bs[3]) / alpha)) / (1 - (1 - lambda) * exp(exp(bs[1] + bs[3]) / alpha));
    }
"

# Write the Stan string to a file
writeLines(cod_regmog_app, con = "cod_regmog_app.stan")

# Prepare the data for Stan (extract covariates)
wmu = as.matrix(test_data[, c("s2","s3")])
data_regmog = list(N = dim(test_data)[1],
                   s = dim(wmu)[2]+1,
                   time = test_data$time,
                   delta = test_data$delta,
                   xb = cbind(1,wmu))

# Run the Bayesian model
regmogfit = stan(file = 'cod_regmog_app.stan', data = data_regmog,
                 chains = 2, iter=8000, warmup=5000, seed=100, cores=2)

# Summarize and extract posterior samples
sumregmog = summary(regmogfit, pars = c("alpha", "lambda","bs"))$summary[,c("mean", "sd", "2.5%","50%", "97.5%", "n_eff", "Rhat")]
round(sumregmog, digits=3)
traceplot(regmogfit, pars = c("alpha", "lambda","bs"))
post = rstan::extract(regmogfit)

comp_model_metrics(fit_model = regmogfit, log_lik_name = "log_lik")

# Calculate the sample proportions for each stage
weigths_stage = round(prop.table(table(test_data$Stage)), 4)

# Create a weighted global posterior for the cure fraction
post_global = (post$p01 * weigths_stage[1]) + (post$p02 * weigths_stage[2]) + (post$p03 * weigths_stage[3])

mean(post$p01)
mean(post$p02)
mean(post$p03)

plot(kmfit)
abline(h=mean(post$p01))

abline(h=mean(post$p02))

abline(h=mean(post$p03))

# Fit Kapaln-Meier estimator
kmfit = survfit(Surv(time, delta) ~ Stage, data = test_data)

plot = ggsurvplot(kmfit, data = test_data,
                   conf.int = F,       
                   xlab = "\n Time (Years)",
                   legend = "top",
                   ylab = "Survival probability \n ",
                   palette = c("black", "black", "black"), 
                   censor = FALSE,
                   legend.title = "") 

pp = plot$plot + 
  
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = ci_p01[1], ymax = ci_p01[2], fill = "steelblue", alpha = 0.15) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = ci_p02[1], ymax = ci_p02[2], fill = "yellow4", alpha = 0.15) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = ci_p03[1], ymax = ci_p03[2], fill = "red4", alpha = 0.15) +
  
  geom_hline(yintercept = mean_p01, color = "steelblue", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = mean_p02, color = "yellow4", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = mean_p03, color = "red4", linetype = "dashed", linewidth = 1) +
  
  guides(color = guide_legend(override.aes = list(
    color = c("steelblue", "yellow4", "red4"),
    alpha = 0.4,       
    linewidth = 5      
  ))) +
  theme(
    axis.text.x = element_text(size = 24),
    axis.text.y = element_text(size = 24),
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    legend.title = element_blank(),
    legend.text = element_text(size = 24)
  ) 

pp

# Calculate mean and 95% credible intervals for the global cure fraction
wpcure = mean(post_global)
wpcure_lw = quantile(post_global, probs = 0.025)
wpcure_up = quantile(post_global, probs = 0.975)
expr_global = sprintf("p == '%.3f (%.3f, %.3f)'", wpcure, wpcure_lw, wpcure_up)

# Fit unstratified KM model
kmfit2 = survfit(Surv(time, delta) ~ 1, data = test_data)

# Plot overall survival curve with global cure fraction overlay
plot2 = ggsurvplot(kmfit2, data = test_data, conf.int = FALSE, xlab = "\n Time (Years)", ylab = "Survival probability \n ", censor = FALSE, legend = "none", palette = "black") 

plot2$plot + 
  annotate("rect", xmin = -1, xmax = max(kmfit2$time), ymin = wpcure_lw, ymax = wpcure_up, fill = "purple4", alpha = 0.15) +
  geom_hline(yintercept = wpcure, color = "purple4", linetype = "dashed", linewidth = 1.2) +
  theme(axis.text.x = element_text(size = 24), axis.text.y = element_text(size = 24), axis.title.x = element_text(size = 26), axis.title.y = element_text(size = 26)) + 
  annotate("text", x = max(kmfit2$time) * 0.65, y = 0.25, label = expr_global, parse = TRUE, size = 9, color = "purple4")

# Calculate medians and 95% CIs for each stage
ci_01 = quantile(post$p01, probs = c(0.025, 0.5, 0.975))
ci_02 = quantile(post$p02, probs = c(0.025, 0.5, 0.975))
ci_03 = quantile(post$p03, probs = c(0.025, 0.5, 0.975))

# Consolidate into a dataframe for plotting
dados_plot = data.frame(
  Estagio = factor(c("I", "II", "III"), levels = c("I", "II", "III")),
  Mediana = c(ci_01[2], ci_02[2], ci_03[2]),
  Limite_Inf = c(ci_01[1], ci_02[1], ci_03[1]),
  Limite_Sup = c(ci_01[3], ci_02[3], ci_03[3])
)

# Format text labels for the plot
expr_01 = sprintf("p[Stage~I] == '%.3f (%.3f, %.3f)'", ci_01[2], ci_01[1], ci_01[3])
expr_02 = sprintf("p[Stage~II] == '%.3f (%.3f, %.3f)'", ci_02[2], ci_02[1], ci_02[3])
expr_03 = sprintf("p[Stage~III] == '%.3f (%.3f, %.3f)'", ci_03[2], ci_03[1], ci_03[3])

# Plot estimates with error bars
ggplot(dados_plot, aes(x = Estagio, y = Mediana)) +
  geom_line(aes(group = 1), color = "gray30", linewidth = 1) +
  geom_errorbar(aes(ymin = Limite_Inf, ymax = Limite_Sup, color = Estagio), width = 0.5, linewidth = 0.9) +
  geom_point(aes(color = Estagio), size = 4) +
  scale_color_manual(values = c("I" = "steelblue", "II" = "yellow4", "III" = "red4")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  annotate("text", x = 0.5, y = c(0.33, 0.23, 0.13), label = c(expr_01, expr_02, expr_03), color = c("steelblue", "yellow4", "red4"), parse = TRUE, hjust = 0, size = 7) +
  labs(x = "\n Stage", y = "Cure probability \n") + 
  theme_minimal() +
  theme(axis.text.x = element_text(size = 24), axis.text.y = element_text(size = 24), axis.title.x = element_text(size = 26), axis.title.y = element_text(size = 26), legend.position="none")


n_amostras = length(post$p01)

# Construct dataframe for absolute differences
df_diff_post = data.frame(
  valor = c(post$p02 - post$p01, post$p03 - post$p01,  post$p03 - post$p02),
  grupo = factor(rep(c("p02 - p01", "p03 - p01", "p03 - p02"), each = n_amostras))
)

# Plot density of absolute differences
ggplot(df_diff_post, aes(x = valor, fill = grupo)) +
  geom_density(alpha = 0.5, color = "black") + 
  theme_minimal() +
  labs(title = "Overlapping Density of Posterior Differences", x = "Difference", y = "Density", fill = "Comparison") +
  scale_fill_manual(values = c("#1f77b4", "#d62728", "green")) 

# Construct dataframe for ratios
df_ratio_post = data.frame(
  valor = c(post$p02/post$p01, post$p03/post$p01,  post$p03/post$p02),
  grupo = factor(rep(c("p02/p01", "p03/p01", "p03/p02"), each = n_amostras))
)

# Plot density of ratios
ggplot(df_ratio_post, aes(x = valor, fill = grupo)) +
  geom_density(alpha = 0.5, color = "black") +  
  theme_minimal() +
  labs(title = "Overlapping Density of Posterior Ratios", x = "Ratio", y = "Density", fill = "Comparison") +
  scale_fill_manual(values = c("#1f77b4", "#d62728", "green")) 

# Summarize the differences and ratios into tables
sum_ratio = df_ratio_post %>%
  group_by(grupo) %>%
  summarise(
    mean = mean(valor),
    ic_inf = quantile(valor, 0.025),
    ic_sup = quantile(valor, 0.975)
  )

sum_diff = df_diff_post %>%
  group_by(grupo) %>%
  summarise(
    mean = mean(valor),
    ic_inf = quantile(valor, 0.025),
    ic_sup = quantile(valor, 0.975)
  )

sum_ratio
sum_diff















