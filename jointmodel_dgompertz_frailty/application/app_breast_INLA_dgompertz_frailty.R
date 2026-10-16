rm(list = ls())
library(dplyr)
library(INLA)
library(INLAjoint)
library(survival)
setwd("~/Desktop/PhD/01_Joint_Model_Defective_Frailty/application")

ca153_long = read.csv2("ca153_long.csv")
cea_long = read.csv2("cea_long.csv")
cancer_patient_final = read.csv2("cancer_patient_final.csv")

cancer_patient_final$cat_age_risk = (cancer_patient_final$CURRENT_AGE_DEID < 65)

cancer_patient_final$PRIOR_MED_TO_MSK = as.factor(cancer_patient_final$PRIOR_MED_TO_MSK)
cancer_patient_final$PRIOR_MED_TO_MSK = relevel(cancer_patient_final$PRIOR_MED_TO_MSK, ref = "No prior medications")

#levels(cancer_patient_final$PRIOR_MED_TO_MSK) 
cancer_patient_final2 = cancer_patient_final %>% filter(PRIOR_MED_TO_MSK!="Unknown") 

kmfit = survfit(Surv(TimeS, Delta)~prior_med,  data = cancer_patient_final2)

dim(cancer_patient_final2)
plot(kmfit)

cea_long2 = cea_long %>%
  filter(id_p %in% cancer_patient_final2$id_p) 

ca153_long2 =
  ca153_long %>%
  filter(id_p %in% cancer_patient_final2$id_p)


cp1 = cea_long2 %>%
  group_by(ID) %>%
  summarise(n_measures = n())

cp2 = ca153_long2 %>%
  group_by(ID) %>%
  summarise(n_measures = n())

summary(cp1$n_measures)
summary(cp2$n_measures)


dim(cea_long2)
dim(ca153_long2)


ctrl = list(
  int.strategy = "ccd", ## "eb" for velocity and ccd to generate posterior samples
  cfg = TRUE,
  cmin = 0,
  
  priorFixed = list(mean = 0, prec = 0.01,
                    mean.intercept = 0, prec.intercept = 0.01),
  priorAssoc = list(mean = 0, prec = 0.01),
  
  priorSRE_ind = list(mean = 0, prec = 0.01),
  
  # Wishart for var-cov of random effects.
  priorRandom = list(r = 5, R = 1), 
  
  # PC prior for frailty. 
  # P(sd frailty > 1) = 0.01
  priorFrailty = list(prior = "pc.prec", param = c(1, 0.01)),
  
  baselineHyper = list(
    alpha = list(prior = "normal", param = c(0, 0.01), initial = -1)
  )
)

## fitting the survival model with frailty
set.seed(127)
survpart_frailty = joint(formSurv = inla.surv(time = TimeS, event = Delta) ~ prior_med  + (1 | id_p), 
                         control=ctrl,
                          basRisk = "dgompertzsurv", 
                          id = "id_p", 
                          dataSurv = cancer_patient_final2)

summary(survpart_frailty, sdcor=T)


## Joint model
## Model with gaussian assumptions for the random effects

formsurv = inla.surv(time = TimeS, event = Delta) ~ prior_med  + (1 | id_p)
form_cea = log_CEA ~ timej + cat_age_risk + (1+timej|id_p)
form_ca153 = log_CA153 ~ timej + HR + (1+timej|id_p)


set.seed(123)
joint_dgompertz = joint(formSurv = formsurv,
                        basRisk = "dgompertzsurv",
                        id = "id_p", dataSurv = cancer_patient_final2,
                        formLong=list(form_cea,form_ca153),
                        family=c("gaussian","gaussian"), dataLong=list(cea_long2,ca153_long2),
                        timeVar="timej", assoc=c("SRE_ind","SRE_ind"),
                        corLong=T,run=T, control=ctrl)

summary(joint_dgompertz, sdcor=T, digits=3)

#saveRDS(joint_dgompertz, file = "3joint_dgompertz.rds")

#joint_dgompertz=readRDS("3joint_dgompertz.rds")

round(joint_dgompertz$summary.fixed,digits=3)
round(joint_dgompertz$summary.hyperpar,digits=3)

# Survival
round(summary(joint_dgompertz, sdcor = TRUE)$SurvEff[[1]][1,], digits = 3)
round(joint_dgompertz$summary.fixed[1,],digits = 3)
round(summary(joint_dgompertz, sdcor = TRUE)$SurvEff[[1]][3,], digits = 3)
round(summary(joint_dgompertz, sdcor = TRUE)$AssocLS, digits = 3)
round(summary(joint_dgompertz, sdcor = TRUE)$ReffListS[[1]],digits = 3)

# Longitudinal
round(summary(joint_dgompertz, sdcor = TRUE)$FixedEff[[1]], digits = 3)
round(summary(joint_dgompertz, sdcor = TRUE)$FixedEff[[2]], digits = 3)

## Var-Cov matrix
round(summary(joint_dgompertz, sdcor = TRUE)$ReffList[[1]], digits = 3)



# Resumo dos resultados
summary(joint_dgompertz, sdcor=TRUE)


## generating samples
set.seed(1234) ## Seed for R
inla.seed = 1234 ## seed for R INLA in C

ns=1000
samples_post = inla.posterior.sample(ns, joint_dgompertz, seed = inla.seed) # For latent

names(joint_dgompertz)
joint_dgompertz$names.fixed
joint_dgompertz$misc$configs$contents

# Defining valid IDs and extracting number of samples
ids_valid1 = 1:nrow(cancer_patient_final2)
ids_valid2 = (nrow(cancer_patient_final2)+1):(2*nrow(cancer_patient_final2))
ids_valid3 = (2*nrow(cancer_patient_final2)+1):(3*nrow(cancer_patient_final2))
ids_valid4 = (3*nrow(cancer_patient_final2)+1):(4*nrow(cancer_patient_final2))


# 1. Extracting scalar parameters across all samples (resulting in vectors of length ns)
alpha_s = sapply(samples_post, function(x) x$hyperpar["alpha parameter for dGompertz-surv[3]"])
mu_s = exp(sapply(samples_post, function(x) x$latent["Intercept_S1:1", 1]))
psi_pm_s = sapply(samples_post, function(x) x$latent["prior_med_S1:1", 1])

## Random effects (latent)
call_int_L1 = paste0("IDIntercept_L1:", ids_valid1)
call_slp_L1 = paste0("IDtimej_L1:", ids_valid2)
call_int_L2 = paste0("IDIntercept_L2:", ids_valid3)
call_slp_L2 = paste0("IDtimej_L2:", ids_valid4)

int_L1_s = sapply(samples_post, function(x) x$latent[call_int_L1, 1])
slp_L1_s = sapply(samples_post, function(x) x$latent[call_slp_L1, 1])
int_L2_s = sapply(samples_post, function(x) x$latent[call_int_L2, 1])
slp_L2_s = sapply(samples_post, function(x) x$latent[call_slp_L2, 1])


dim(int_L1_s)

# Result: A matrix with 1000 rows (samples) and length(times) columns.
# Extracting vector parameters (resulting in matrices of 915 rows x 1000 columns)
# Each column represents one complete posterior sample for all 915 patients

frailty_s = sapply(samples_post, function(x) x$latent[paste0("IDIntercept_S1:", ids_valid1), 1])

gamma_b01_s = sapply(samples_post, function(x) x$hyperpar["Beta for SRE_Intercept_L1_S1"])
gamma_b11_s = sapply(samples_post, function(x) x$hyperpar["Beta for SRE_timej_L1_S1"])
gamma_b02_s = sapply(samples_post, function(x) x$hyperpar["Beta for SRE_Intercept_L2_S1"])
gamma_b12_s = sapply(samples_post, function(x) x$hyperpar["Beta for SRE_timej_L2_S1"])


sre_int_L1_weighted = t(t(int_L1_s) * gamma_b01_s)
sre_slp_L1_weighted = t(t(int_L1_s) * gamma_b11_s)
sre_int_L2_weighted = t(t(int_L1_s) * gamma_b02_s)
sre_slp_L2_weighted = t(t(int_L1_s) * gamma_b12_s)

# Covariate matrix W
W = as.matrix(cancer_patient_final2[, "prior_med"])

# Matrix multiplication: W (915 x 2) %*% psi_mat (2 x 1000) 
#W_eff = W %*% psi_mat
W_eff = W %*% psi_pm_s

inner_sum = frailty_s + W_eff + 
  sre_int_L1_weighted + sre_slp_L1_weighted + 
  sre_int_L2_weighted + sre_slp_L2_weighted

ratio_s = mu_s / alpha_s

inner_exp_scaled = t( t(exp(inner_sum)) * ratio_s )

p_cures_t = exp(inner_exp_scaled)

# Transposing to achieve the final requested dimension (1000 rows x 915 columns)
p_cures = t(p_cures_t)

dim(p_cures) 

head(p_cures)
#p_est=colMeans(p_cures)
p_est=apply(p_cures, MARGIN = 2, FUN = mean)
p_pop_est=median(p_est)

p_quantiles=quantile(p_est, probs = c(0.025,0.975))

kmfit=survfit(Surv(TimeS, Delta)~1, data=cancer_patient_final2)
plot(kmfit)
abline(h=p_pop_est, col="red",lwd=2)
abline(h=quantile(p_est, probs = c(0.025,0.975)), lty=2)


library(ggplot2)
library(broom)

kmfit = survfit(Surv(TimeS, Delta) ~ prior_med, data = cancer_patient_final2)
kmfit_tidy = tidy(kmfit)

text = sprintf("%.3f (%.3f; %.3f)", 
               p_pop_est, 
               p_quantiles[1], p_quantiles[2])

ggplot(kmfit_tidy, aes(x = time, y = estimate, color = strata)) +
  geom_step(linewidth = 1.5) +
  
  scale_color_manual(
    name = "Prior Medication",
    values = c("prior_med=1" = "darkred", "prior_med=0" = "steelblue"),
    labels = c("prior_med=1" = "Yes", "prior_med=0" = "No")
  ) +
  
  scale_y_continuous(
    limits = c(0,1),
    breaks = seq(0, 1, 0.25)) +
  scale_x_continuous(
    limits = c(0,10)) +
  
  labs(
    x = "Time (Years)", 
    y = "Survival probability"
  )  + 
  theme_bw(base_size = 14) +
  theme(
    legend.position = c(0.98, 0.98), 
    legend.justification = c(1, 1), 
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

## Ploting the individual cures
  
idx_group1 = W[,1]==0
color_group = "steelblue"  #"darkred" 
post_pest_group1 = p_est[idx_group1]
ci_pest_group1 = t(apply(p_cures[,idx_group1], MARGIN=2, FUN = quantile, probs = c(0.025,0.975)))
  

df_group1 = data.frame(
  Subject = 1:sum(idx_group1),              
  Probability = post_pest_group1,      
  Lower = ci_pest_group1[,1],   
  Upper = ci_pest_group1[,2]    
)

group_prob_cure = median(df_group1$Probability)
group_prob_cure_ci = quantile(df_group1$Probability, probs=c(0.025,0.975))

label_text = paste0(
  "hat(p) == '", sprintf("%.3f", group_prob_cure), "'", 
  "~'(", sprintf("%.3f", group_prob_cure_ci[1]), "; ", sprintf("%.3f", group_prob_cure_ci[2]), ")'"
)

# graph
ggplot(df_group1, aes(x = Subject, y = Probability)) +
  
  geom_errorbar(aes(ymin = Lower, ymax = Upper), color = "gray55", width = 0) +
  geom_point(color = "black", size = 2.5) +
  annotate("rect", xmin = -Inf, xmax = Inf, 
           ymin = ci_pest_group1[1], ymax = ci_pest_group1[2], 
           fill = color_group, alpha = 0.15) +
  
  geom_hline(yintercept = group_prob_cure_ci[1], color = color_group, linetype = "dashed", linewidth = 0.9, alpha = 0.6) +
  geom_hline(yintercept = group_prob_cure_ci[2], color = color_group, linetype = "dashed", linewidth = 0.9, alpha = 0.6) +
  
  # Horizontal line
  geom_hline(yintercept = group_prob_cure, color = color_group, linewidth = 1.4) +
  
  ylim(c(0,1)) +
  
  # Anotating the man and the credible interval
  annotate(
    "text", 
    x = Inf, y = Inf, 
    label = label_text, 
    parse = TRUE, 
    hjust = 1.05, vjust = 1.5, 
    color = color_group, 
    size = 7 
  ) +
  
  theme_bw(base_size = 14) +
  
  labs(
    y = "Probability of cure"
  ) +
  
  theme(
    legend.position = "top",
    legend.title = element_text(size = 22),
    legend.text  = element_text(size = 22),
    strip.text   = element_text(size = 22),
    panel.grid.minor = element_blank(),
    
    axis.title.x = element_text(size = 22, margin = margin(t = 15)), 
    axis.title.y = element_text(size = 22,  margin = margin(r = 15)), 
    
    axis.text.x = element_text(size = 20, color = "black"), 
    axis.text.y = element_text(size = 20, color = "black"),
    
    axis.ticks = element_line(linewidth = 1, color = "black")                   
  )

## [The distribution of survival probabilities over time and the risk reduction]

seq_tc = seq(0, 40, length=200)

# The %o% operator calculates the outer product, and exp() is vectorized
hyper_samples = inla.hyperpar.sample(ns,result =joint_dgompertz, improve.marginals = T)
hyper_alpha_samples = hyper_samples[,3]

mean(hyper_alpha_samples)
sd(hyper_alpha_samples)
quantile(hyper_alpha_samples, probs = c(0.025, 0.975))

summary(hyper_alpha_samples)

#hyper_alpha_samples=alphas

c_mat = 1-exp(hyper_alpha_samples %o% seq_tc)

# Calculate statistics across the rows (the 100 alphas) for each column (the 200 time points)
c_mean   = apply(c_mat, 2, mean)
c_median = apply(c_mat, 2, median)
c_lower  = apply(c_mat, 2, quantile, probs = 0.025) # 2.5% lower bound
c_upper  = apply(c_mat, 2, quantile, probs = 0.975) # 97.5% upper bound

# Combine everything into a neat data frame
plot_data = data.frame(
  Time   = seq_tc,
  Mean   = c_mean,
  Median = c_median,
  Lower  = c_lower,
  Upper  = c_upper
)

# Calculate your max observed survival time
max_time = max(cancer_patient_final$TimeS)

library(scales) # Load this to easily format percentages

ggplot(plot_data, aes(x = Time)) +

  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3) +
  
  geom_line(aes(y = Mean, color = "Mean"), linewidth = 1) +
  geom_line(aes(y = Median, color = "Median"), linewidth = 1) +
  
  geom_vline(xintercept = max_time, linetype = "dashed", color = "gray50", linewidth = 1) +
  
  scale_color_manual(values = c("Mean" = "black", "Median" = "darkviolet")) +
  
  scale_x_continuous(breaks = seq(0, 40, by = 5)) +

  scale_y_continuous(breaks = seq(0, 1, by = 0.1), labels = scales::percent) +
  
  labs(
    x = "Time (Years)",
    y = "Risk reduction percentage (1-c)%",
    color = NULL 
  ) +
  theme_bw(base_size = 14) +
  
  theme(
    legend.position = c(0.12, 0.9),
    legend.title = element_text(size = 22),
    legend.background = element_rect(fill = NA, color = NA), 
    legend.text  = element_text(size = 18),
    strip.text   = element_text(size = 22),
    panel.grid.minor = element_blank(),
    
    
    axis.title.x = element_text(size = 22, margin = margin(t = 15)), 
    axis.title.y = element_text(size = 22,  margin = margin(r = 15)), 
    
    axis.text.x = element_text(size = 20, color = "black"), 
    axis.text.y = element_text(size = 20, color = "black"),
    
    axis.ticks = element_line(linewidth = 1, color = "black") 
  )


## Compute the relative risk reduction for follow-up times
seq_tc = c(0.5, 1, 3, 5, 7, 10, 15, 20, 30, 40, 50, 60)
minusc_mat = 1-exp(alpha_s %o% seq_tc)

mean(minusc_mat[,1])
sd(minusc_mat[,6])
quantile(minusc_mat[,12], probs = c(0.025,0.5, 0.975))

#alpha_s

sum_df = data.frame(
  Years  = seq_tc,
  Mean   = apply(minusc_mat, 2, mean),
  Median = apply(minusc_mat, 2, median),
  Lower  = apply(minusc_mat, 2, quantile, probs = 0.025),
  Upper  = apply(minusc_mat, 2, quantile, probs = 0.975)
)



ggplot(sum_df, aes(x = Years)) +
  
  geom_errorbar(
    aes(ymin = Lower, ymax = Upper),
    width = 0,
    color = "grey50",
    linewidth = 1
  ) +
  
  geom_line(
    aes(y = Median),
    color = "steelblue",
    alpha=0.5,
    linewidth = 1.2,
    linetype = "solid"
  ) +
  
  # Median points
  geom_point(
    aes(y = Median, shape = "Median"),
    color = "steelblue",
    size = 4.5
  ) +
  
  geom_vline(
    xintercept = max(cancer_patient_final$TimeS),
    linetype = "dashed",
    color = "black",
    linewidth = 1.2
  ) +
  
  scale_shape_manual(
    values = c("Mean" = 17, "Median" = 16)
  ) +
  
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), labels = scales::percent)+
  
  scale_x_continuous(
    breaks = sum_df$Years[sum_df$Years != 0.5],
    labels = function(x) as.integer(x)
  ) +
  
  labs(
    x = "Time (Years)",
    y = "Risk reduction 100(1-c)%",
    shape = NULL
  )  +
  theme_bw(base_size = 14) +
  

  theme(
    legend.position = c(0.85, 0.1),
    legend.title = element_text(size = 22),
    legend.background = element_rect(fill = NA, color = NA), 
    legend.text  = element_text(size = 20),
    strip.text   = element_text(size = 22),
    panel.grid.minor = element_blank(),
    
    
    axis.title.x = element_text(size = 22, margin = margin(t = 15)), 
    axis.title.y = element_text(size = 22,  margin = margin(r = 15)), 
    
    axis.text.x = element_text(size = 20, color = "black"), 
    axis.text.y = element_text(size = 20, color = "black"),
    
    axis.ticks = element_line(linewidth = 1, color = "black") 
  )

arrays_survtc = array(NA,dim=c(ns, nrow(W), length(seq_tc)))

for(s in 1:ns){
  for(i in 1:nrow(W)){
    for(j in 1:length(seq_tc)){
      arrays_survtc[s,i,j] = exp(inner_exp_scaled[i,s]*(1 - exp(seq_tc[j]*alpha_s[s])))
    }
  }
}


dim(arrays_survtc)

mean(apply(arrays_survtc[,,1], 2, mean))
quantile(apply(arrays_survtc[,,1], 2, mean), probs = c(0.025, 0.975))

mean(apply(arrays_survtc[,,5], 2, mean))
quantile(apply(arrays_survtc[,,5], 2, mean), probs = c(0.025, 0.975))

boxplot(apply(arrays_survtc[,,5], 2, mean))


lista_medias = lapply(1:length(seq_tc), function(j) {
  apply(arrays_survtc[,,j], 2, mean)
})


names(lista_medias) = 1:length(seq_tc)

df_box = stack(lista_medias)

colnames(df_box) = c("Media", "Grupo")

df_box$Grupo = factor(df_box$Grupo, levels = 1:length(seq_tc))

my_blues = colorRampPalette(c("#E0F3F8", "#1E88F9"))(length(seq_tc))

ggplot(df_box, aes(x = Grupo, y = Media, fill = Grupo)) +
  geom_boxplot(show.legend = FALSE) +
  labs(
    x = "Time (Years)",
    y = "Survival probability") +
  scale_x_discrete(
    labels = seq_tc
  ) +
  geom_hline(yintercept = p_pop_est, linetype = "dashed", color = "darkred", linewidth = 1.2) +
  
  annotate(
    "text", 
    x = 11,               
    y = Inf,              
    label = paste0("hat(p) == ", sprintf("%.3f", p_pop_est)), 
    vjust = 2,            
    hjust = 0.7,           
    color = "darkred",                
    fontface = "bold",
    size = 8,
    parse = TRUE          
  ) +
  theme_bw(base_size = 14) +
  

  theme(
    legend.position = c(0.12, 0.5),
    legend.title = element_text(size = 22),
    legend.background = element_rect(fill = NA, color = NA), 
    legend.text  = element_text(size = 20),
    strip.text   = element_text(size = 22),
    panel.grid.minor = element_blank(),
    
    axis.title.x = element_text(size = 22, margin = margin(t = 15)), 
    axis.title.y = element_text(size = 22,  margin = margin(r = 15)), 
    
    axis.text.x = element_text(size = 20, color = "black"), 
    axis.text.y = element_text(size = 20, color = "black"),
    
    axis.ticks = element_line(linewidth = 1, color = "black") 
  ) +
  scale_fill_manual(values = my_blues)


## Showing examples of Collolary 1
p_est[W[,1]==0] 
p_est[W[,1]==1]

length(p_est[W[,1]==0]); length(p_est[W[,1]==1])
which(p_est[W[,1]==0] > 0.2) 
#p_est[525]


## viasulaizing the longitudinal biomarkers
id_chosen= 305 
cancer_patient_final2[id_chosen,] %>% select(TimeS, Delta, PRIOR_MED_TO_MSK)

p_est[id_chosen]
id_pchosen = unique(cancer_patient_final2$id_p)[id_chosen]
dadoslong1_idp = cea_long2[cea_long2$id_p == id_pchosen, ]
dadoslong2_idp = ca153_long2[ca153_long2$id_p == id_pchosen,]

plot(dadoslong2_idp$timej, dadoslong2_idp$log_CEA)
plot(dadoslong2_idp$timej, dadoslong2_idp$log_CA153)

last_time_idp = as.numeric(cancer_patient_final2[id_chosen,] %>% select(TimeS))

time_eval = seq(0,60, 0.01)

idp_survtc = matrix(NA,nrow=ns,ncol=length(time_eval))


for(s in 1:ns){
    for(j in 1:length(time_eval)){
      idp_survtc[s,j] = exp(inner_exp_scaled[id_chosen,s]*(1 - exp(time_eval[j]*alpha_s[s])))
    }
}

head(idp_survtc)

surv_mean_idp = colMeans(idp_survtc)

plot(time_eval, surv_mean_idp, ylim=c(0,1), type="l", lwd=2)
abline(h=p_est[id_chosen], col="red")


df1 = data.frame(
  time_eval = time_eval,
  surv_mean_idp = surv_mean_idp
)

ggplot(df1, aes(x = time_eval, y = surv_mean_idp)) +
  geom_line(linewidth = 1.5) +
  geom_hline(
    yintercept = p_est[id_chosen],
    color = "darkred",
    linewidth = 1.2
  ) +
  ylim(0, 1) + 
  annotate(
    "text", 
    x = Inf, y = Inf, 
    label = paste0(
      "hat(p)[", id_chosen, "] == \"",
      sprintf("%.3f", as.numeric(p_est[id_chosen])),
      "\""
    ),
    parse = TRUE, 
    hjust = 1.2, vjust = 1.5,
    color = "darkred", 
    size = 8
  ) +
  geom_vline(
    xintercept = last_time_idp, 
    color = "black", 
    linewidth = 1.2, 
    linetype = "dashed",
    alpha = 0.6
  ) +
  scale_x_continuous(
    breaks = function(x) {
      b = unique(c(pretty(x), last_time_idp))
      b = b[!(b %in% c(0,1))] 
      sort(b)
    },
    labels = function(b) {
      ifelse(
        b == last_time_idp, 
        sprintf("%.2f", b), 
        as.character(b)     
      )
    }
  ) +
  labs(
    x = "Time (Years)", 
    y = bquote(S[.(id_chosen)]~"(Time)")
  ) + 
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 22),
    legend.text  = element_text(size = 22),
    strip.text   = element_text(size = 22),
    panel.grid.minor = element_blank(),
    
    axis.title.x = element_text(size = 22, margin = margin(t = 15)), 
    axis.title.y = element_text(size = 22,  margin = margin(r = 15)), 
    
    axis.text.x = element_text(size = 20, color = "black"), 
    axis.text.y = element_text(size = 20, color = "black"),
    
    axis.ticks = element_line(linewidth = 1, color = "black")
  )

## survival function based on corollary 2
c_eval = seq(1,0, length=200)
idp_survtc2 = matrix(NA,nrow=ns,ncol=length(c_eval))
for(s in 1:ns){
for(j in 1:length(c_eval)){
idp_survtc2[s,j] = p_cures[s,id_chosen]^(1-c_eval[j])
}
}
head(idp_survtc2)
colMeans(idp_survtc2)
surv_mean_idp2 = colMeans(idp_survtc2)

plot(1-c_eval, surv_mean_idp2, ylim=c(0,1), type="l", lwd=2)
abline(h=p_est[id_chosen], col="red")

c_lastime_idp = rep(NA,length(alpha_s))
for(s in 1:length(alpha_s)){
  c_lastime_idp[s] = exp(last_time_idp*alpha_s[s])
}


mean_val = mean(1 - c_lastime_idp)

ci_lower = quantile(1 - c_lastime_idp, probs = 0.0275) 
ci_upper = quantile(1 - c_lastime_idp, probs = 0.975)


df2 = data.frame(
  time_eval = 1 - c_eval,
  surv_mean_idp = surv_mean_idp2
)

p_lab = sprintf("%.3f", as.numeric(p_est[id_chosen]))


ggplot(df2, aes(x = time_eval, y = surv_mean_idp)) +
  geom_line(linewidth = 1.5) +
  geom_hline(
    yintercept = p_est[id_chosen],
    color = "darkred",
    linewidth = 1.2
  ) +
  ylim(0, 1) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.10), labels = scales::percent) +
  
  annotate(
    "text",
    x = Inf, y = Inf,
    label = "- - (1-c)%",
    hjust = 1.12, vjust = 5.2, 
    color = "grey60",
    size = 8
  ) +
  
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0("hat(p)[147] == \"", p_lab, "\""),
    parse = TRUE,
    hjust = 1.2, vjust = 2, 
    color = "darkred",
    size = 8
  ) +

  annotate(
    "rect",
    xmin = ci_lower, 
    xmax = ci_upper,
    ymin = -Inf, 
    ymax = Inf,
    fill = "grey10", 
    alpha = 0.2  
  ) +
  
  geom_vline(
    xintercept = mean_val, 
    color = "black", 
    linewidth = 1.2, 
    linetype = "dashed",
    alpha = 0.6
  ) +
  
  labs(
    x = "Risk reduction percentage (1-c)%",
    y = expression(S[147]==p[147]^(1-c))
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 22),
    legend.text  = element_text(size = 22),
    strip.text   = element_text(size = 22),
    panel.grid.minor = element_blank(),
    
    axis.title.x = element_text(size = 22, margin = margin(t = 15)),
    axis.title.y = element_text(size = 22, margin = margin(r = 15)),
    
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    
    axis.ticks = element_line(linewidth = 1, color = "black")
  )


mean_val = mean(1 - c_lastime_idp)
ci_lower = quantile(1 - c_lastime_idp, probs = 0.0275) 
ci_upper = quantile(1 - c_lastime_idp, probs = 0.975)

df2 = data.frame(
  time_eval = 1 - c_eval,
  surv_mean_idp = surv_mean_idp2
)

y_mean = approx(x = df2$time_eval, y = df2$surv_mean_idp, xout = mean_val)$y
y_lower = approx(x = df2$time_eval, y = df2$surv_mean_idp, xout = ci_lower)$y
y_upper = approx(x = df2$time_eval, y = df2$surv_mean_idp, xout = ci_upper)$y

p_lab = sprintf("%.3f", as.numeric(p_est[id_chosen]))

my_y_breaks = round(sort(unique(c(seq(0, 1, by = 0.25), y_mean, y_lower, y_upper))), 3)

my_y_breaks = my_y_breaks[my_y_breaks != 0.75]

ggplot(df2, aes(x = time_eval, y = surv_mean_idp)) +
  annotate(
    "rect",
    xmin = ci_lower, xmax = ci_upper,
    ymin = -Inf, ymax = Inf,
    fill = "grey10", alpha = 0.2  
  ) +
  
  geom_line(linewidth = 1.5) +
  geom_hline(yintercept = p_est[id_chosen], color = "darkred", linewidth = 1.2) +
  
  geom_segment(aes(x = 0, xend = mean_val, y = y_mean, yend = y_mean), 
               color = "black", linetype = "dashed", linewidth = 1.2, alpha = 0.6) +
  geom_segment(aes(x = 0, xend = ci_lower, y = y_lower, yend = y_lower), 
               color = "grey40", linetype = "dotted", linewidth = 1) +
  geom_segment(aes(x = 0, xend = ci_upper, y = y_upper, yend = y_upper), 
               color = "grey40", linetype = "dotted", linewidth = 1) +
  
  geom_point(aes(x = mean_val, y = y_mean), size = 4, color = "black") +
  geom_point(aes(x = ci_lower, y = y_lower), size = 3, color = "grey40") +
  geom_point(aes(x = ci_upper, y = y_upper), size = 3, color = "grey40") +
  
  geom_vline(xintercept = mean_val, color = "black", linewidth = 1.2, linetype = "dashed", alpha = 0.6) +
  
  scale_y_continuous(limits = c(0, 1), breaks = my_y_breaks) +
  
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.10),
    labels = function(x) sprintf("%.0f", x * 100)
  ) +
  
  annotate(
    "text", x = Inf, y = Inf, 
    label = paste0("hat(p)[", id_chosen, "] == \"", p_lab, "\""), 
    parse = TRUE, hjust = 1.16, vjust = 3, color = "darkred", size = 8
  ) +
  
  annotate(
    "text", x = Inf, y = Inf, 
    label = sprintf("%.3f (%.3f; %.3f)", 100*mean_val, 100*ci_lower, 100*ci_upper), 
    parse = FALSE, hjust = 1.08, vjust = 2.7, 
    color = "grey40", 
    size = 8
  ) +
  
  labs(
    x = "Risk reduction percentage 100(1-c)%",
    y = bquote(S[.(id_chosen)] == p[.(id_chosen)]^(1-c)) 
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 22),
    legend.text  = element_text(size = 22),
    strip.text   = element_text(size = 22),
    panel.grid.minor = element_blank(),
    
    axis.title.x = element_text(size = 22, margin = margin(t = 15)),
    axis.title.y = element_text(size = 22, margin = margin(r = 15)),
    
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    
    axis.ticks = element_line(linewidth = 1, color = "black")
  )


################
## Longitudinal plots
################

# ids_valid1 = 1:nrow(cancer_patient_final2)
# ids_valid2 = (nrow(cancer_patient_final2)+1):(2*nrow(cancer_patient_final2))
# ids_valid3 = (2*nrow(cancer_patient_final2)+1):(3*nrow(cancer_patient_final2))
# ids_valid4 = (3*nrow(cancer_patient_final2)+1):(4*nrow(cancer_patient_final2))

min(ids_valid1);max(ids_valid1)

#joint_dgompertz$misc$configs$contents
#summary(joint_dgompertz, sdcor = TRUE)$ReffList[[1]]

## INLA make copies of the object
## taking the 
call_int_L1 = paste0("IDIntercept_L1:", ids_valid1)
call_slp_L1 = paste0("IDtimej_L1:", ids_valid2)
call_int_L2 = paste0("IDIntercept_L2:", ids_valid3)
call_slp_L2 = paste0("IDtimej_L2:", ids_valid4)

int_L1_s = sapply(samples_post, function(x) x$latent[call_int_L1, 1])
slp_L1_s = sapply(samples_post, function(x) x$latent[call_slp_L1, 1])
int_L2_s = sapply(samples_post, function(x) x$latent[call_int_L2, 1])
slp_L2_s = sapply(samples_post, function(x) x$latent[call_slp_L2, 1])

dim(int_L1_s)

beta01_s = sapply(samples_post, function(x) x$latent["Intercept_L1:1",1])
beta11_s = sapply(samples_post, function(x) x$latent["timej_L1:1",1])
beta21_s = sapply(samples_post, function(x) x$latent["cat_age_risk_L1:1",1])

beta02_s = sapply(samples_post, function(x) x$latent["Intercept_L2:1",1])
beta12_s = sapply(samples_post, function(x) x$latent["timej_L2:1",1])
beta22_s = sapply(samples_post, function(x) x$latent["HR_L2:1",1])


## L1 (CEA)
# sd_int1_amostra = apply(int_L1_s, 2, sd)
# mean(sd_int1_amostra)
# sd(sd_int1_amostra)
# quantile(sd_int1_amostra, c(0.025, 0.5, 0.975))
# 
# sd_slp1_amostra = apply(slp_L1_s, 2, sd)
# mean(sd_slp1_amostra)
# sd(sd_slp1_amostra)
# quantile(sd_slp1_amostra, c(0.025, 0.5, 0.975))

## L2 (CA 15-3)
# summary(joint_dgompertz, sdcor = TRUE)$ReffList[[1]]
# 
# sd_int2_amostra <- apply(int_L2_s, 2, sd)
# mean(sd_int2_amostra)
# sd(sd_int2_amostra)
# quantile(sd_int2_amostra, c(0.025, 0.5, 0.975))
# 
# sd_slp2_amostra <- apply(slp_L2_s, 2, sd)
# mean(sd_slp2_amostra)
# sd(sd_slp2_amostra)
# quantile(sd_slp2_amostra, c(0.025, 0.5, 0.975))


id_chosen = 620 #(prior medication)#305 #(no prior medication) 
id_pchosen = unique(cancer_patient_final2$id_p)[id_chosen]

dadoslong1_idp = cea_long2[cea_long2$id_p == id_pchosen, ]
dadoslong2_idp = ca153_long2[ca153_long2$id_p == id_pchosen, ]

# Correção 1 e 2: Extrair como vetor e pegar apenas o 1º valor (assumindo que a idade/risco é fixa por paciente)
X11 = dadoslong1_idp %>% pull(cat_age_risk) %>% .[1]
X12 = dadoslong2_idp %>% pull(HR) %>% .[1]

# max(cancer_patient_final$TimeS)
timej_seq = seq(0, 10, length = 100)

## Longitudinal 1
m_cea = sapply(timej_seq, function(t) {
  (beta01_s + int_L1_s[id_chosen, ]) + 
    (beta11_s + slp_L1_s[id_chosen, ]) * t + 
    (beta21_s * as.numeric(X11)) 
})


dim(m_cea)

mean_id1 = apply(m_cea, MARGIN=2, FUN = mean)
cilower_id1 = apply(m_cea, MARGIN=2, FUN = quantile, probs = 0.025)
ciupper_id1 = apply(m_cea, MARGIN=2, FUN = quantile, probs = 0.975)

## making graphs
last_time_idp = as.numeric(cancer_patient_final2[id_chosen,] %>% select(TimeS))
last_time_idp = round(last_time_idp,2)


pred_data = data.frame(
  timej_seq = timej_seq,
  mean_id1 = mean_id1,
  cilower_id1 = cilower_id1,
  ciupper_id1 = ciupper_id1
)


# meus_breaks = c(3, 4, last_time_idp, 6, 7, 8, 9, 10)
# meus_labels = c("3", "4",last_time_idp, "6", "7", "8", "9","10")

meus_breaks = c(0, 0.5, 1, last_time_idp,2, 3, 4, 5)
meus_labels = c("0","0.5", "1", last_time_idp, "2", "3", "4", "5")


ggplot(data = dadoslong1_idp, aes(x = timej, y = log_CEA)) +

  geom_point(size=3) +
  
  geom_ribbon(
    data = pred_data,
    aes(
      x = timej_seq,
      ymin = cilower_id1,
      ymax = ciupper_id1
    ),
    inherit.aes = FALSE,
    fill = "steelblue",
    alpha = 0.25
  ) +
  
  geom_line(
    data = pred_data,
    aes(x = timej_seq, y = mean_id1),
    inherit.aes = FALSE,
    color = "steelblue",
    linewidth = 1.2
  ) +
  
  coord_cartesian( ylim = c(0, 1)) +
  
  geom_vline(
    xintercept = last_time_idp,
    linetype = "dashed",
    color = "grey50",
    linewidth = 1.2
  ) +
  
  scale_x_continuous(
    limits = c(0,5),
    breaks = meus_breaks,
    labels = meus_labels
  )  +

  labs(
    x = "Time (Years)", 
    y = bquote("Log CEA"[.(id_chosen)]) 
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    legend.position = "top",
    legend.title = element_text(size = 22),
    legend.text  = element_text(size = 22),
    strip.text   = element_text(size = 22),
    panel.grid.minor = element_blank(),
    
    axis.title.x = element_text(size = 22, margin = margin(t = 15)),
    axis.title.y = element_text(size = 22, margin = margin(r = 15)),
    
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    
    axis.ticks = element_line(linewidth = 1, color = "black")
  )


## Longitudinal 2
m_ca153 = sapply(timej_seq, function(t) {
  (beta02_s + int_L2_s[id_chosen, ]) + 
    (beta12_s + slp_L2_s[id_chosen, ]) * t + 
    (beta22_s * as.numeric(X12)) 
})


dim(m_ca153)

mean_id2 = apply(m_ca153, MARGIN=2, FUN = mean)
cilower_id2 = apply(m_ca153, MARGIN=2, FUN = quantile, probs = 0.025)
ciupper_id2 = apply(m_ca153, MARGIN=2, FUN = quantile, probs = 0.975)


## making graphs
plot(dadoslong2_idp$timej, dadoslong2_idp$log_CA153, xlim = c(0,10))
lines(timej_seq, mean_id2, col = "steelblue")
lines(timej_seq, cilower_id2, col = "steelblue")
lines(timej_seq, ciupper_id2, col = "steelblue")


pred_data2 = data.frame(
  timej_seq = timej_seq,
  mean_id2 = mean_id2,
  cilower_id2 = cilower_id2,
  ciupper_id2 = ciupper_id2
)


ggplot(data = dadoslong2_idp, aes(x = timej, y = log_CA153)) +
  
  geom_point(size=3) +
  
  geom_ribbon(
    data = pred_data2,
    aes(
      x = timej_seq,
      ymin = cilower_id2,
      ymax = ciupper_id2
    ),
    inherit.aes = FALSE,
    fill = "steelblue",
    alpha = 0.25
  ) +
  
  geom_line(
    data = pred_data2,
    aes(x = timej_seq, y = mean_id2),
    inherit.aes = FALSE,
    color = "steelblue",
    linewidth = 1.2
  ) +
  
  scale_x_continuous(
    limits = c(0,5),
    breaks = meus_breaks,
    labels = meus_labels
  )  +
  
  theme_bw(base_size = 14) +
  
  geom_vline(
    xintercept = last_time_idp,
    linetype = "dashed",
    color = "grey50",
    linewidth = 1.2
  ) +
  
  labs(
    x = "Time (Years)", 
    y = bquote("Log CA 15-3"[.(id_chosen)]) 
  ) +
  
  theme(
    legend.position = "top",
    legend.title = element_text(size = 22),
    legend.text  = element_text(size = 22),
    strip.text   = element_text(size = 22),
    panel.grid.minor = element_blank(),
    
    axis.title.x = element_text(size = 22, margin = margin(t = 15)),
    axis.title.y = element_text(size = 22, margin = margin(r = 15)),
    
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    
    axis.ticks = element_line(linewidth = 1, color = "black")
  )




