## -----
## Code of application
## Epileptic data
## A NEW APPROACH FOR JOINT CURE MODELLING USING THE
## DEFECTIVE GOMPERTZ DISTRIBUTION 
## -----

## [Libraries]
## Please install the INLA package in your machine (https://www.r-inla.org/download-install)
library(survival)
library(joineRML)
library(INLA)
library(INLAjoint)

library(dplyr)
library(ggplot2)

## 1. [Load data]
set.seed(1)
data=joineRML::epileptic.qol
data$time = data$time/365.25
data$with.time = data$with.time/365.25

epileptic.l = data[, -c(2,4,9)]
epileptic.s = data[!duplicated(data$id), c(1:3, 9)] # filtering the survival time and censoring

## 2. [Descriptive analysis]

head(epileptic.l)
head(epileptic.s) # CBZ (Carbamazepine) == 0 e LTG (Lamotrigine) ==1

dim(epileptic.l)
dim(epileptic.s)

## longitudinal measures per patient
table(table(epileptic.l$id))
prop.table(table(table(epileptic.l$id)))

## proportion of event and censoring
prop.table(table(epileptic.s$with.status2))


## The kaplan-Meier curve
# kmfit=survfit(Surv(with.time, with.status2)~trt,data = epileptic.s)
# plot(kmfit, xlab = "Time (Years)", ylab = "Survival probability",
#      main = "", ylim = c(0, 1), conf.int=F)
# grid()
# 
# dev.off()

## Analyzing the characteristics of the biomarkers
## for the poisson, mean=variance
mean(epileptic.l$anxiety,na.rm=T);var(epileptic.l$anxiety,na.rm=T)
mean(epileptic.l$aep,na.rm=T);var(epileptic.l$aep,na.rm=T) # breaks the assumption
mean(epileptic.l$depress,na.rm=T);var(epileptic.l$depress,na.rm=T)


disp_intra_trt = epileptic.l %>%
  group_by(trt, id) %>% 
  summarise(
    n_medidas = n(),
    media_aep = mean(depress, na.rm = TRUE),
    var_aep = var(depress, na.rm = TRUE),
    .groups = "drop" # Remove o agrupamento após o sumário
  ) %>%
  filter(n_medidas > 1) %>%
  mutate(indice_dispersao = var_aep / media_aep)


disp_intra_trt_anx = epileptic.l %>%
  group_by(trt, id) %>% 
  summarise(
    n_medidas = n(),
    media_anx = mean(anxiety, na.rm = TRUE),
    var_anx = var(anxiety, na.rm = TRUE),
    .groups = "drop" 
  ) %>%
  filter(n_medidas > 1) %>%
  mutate(
    indice_dispersao = var_anx / media_anx,
    diff_abs = abs(var_anx - media_anx) 
  )


ggplot(disp_intra_trt_anx, aes(x = as.factor(trt), y = diff_abs, fill = as.factor(trt))) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5, color = "darkgray") +
  labs(
    x = "Treatment",
    y = "|Variance - Mean| intra-subject for Anxiety"
  ) +
  scale_fill_manual(values = c("darkred", "steelblue")) +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text   = element_text(size = 22),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 22, margin = margin(t = 15)), 
    axis.title.y = element_text(size = 22,  margin = margin(r = 15)), 
    
    axis.text.x = element_text(size = 20, color = "black"), 
    axis.text.y = element_text(size = 20, color = "black"),
    
    axis.ticks = element_line(linewidth = 1, color = "black")
  )

######

disp_intra_trt_dep = epileptic.l %>%
  group_by(trt, id) %>% 
  summarise(
    n_medidas = n(),
    media_dep = mean(depress, na.rm = TRUE),
    var_dep = var(depress, na.rm = TRUE),
    .groups = "drop" 
  ) %>%
  filter(n_medidas > 1) %>%
  mutate(
    indice_dispersao = var_dep / media_dep,
    diff_abs = abs(var_dep - media_dep)
  )


ggplot(disp_intra_trt_dep, aes(x = as.factor(trt), y = diff_abs, fill = as.factor(trt))) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5, color = "darkgray") +
  labs(
    x = "Treatment",
    y = "|Variance - Mean| intra-subject for Depression"
  ) +
  scale_fill_manual(values = c("darkred", "steelblue")) +
  scale_y_continuous(breaks = seq(0, 120, by = 10)) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text   = element_text(size = 22),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 22, margin = margin(t = 15)), 
    axis.title.y = element_text(size = 22,  margin = margin(r = 15)), 
    
    axis.text.x = element_text(size = 20, color = "black"), 
    axis.text.y = element_text(size = 20, color = "black"),
    
    axis.ticks = element_line(linewidth = 1, color = "black")
  )



# compute stats of biomarkers
desc_biomarker = function(data, biomarker, id_col, value_col) {
  
  x = data[[value_col]]
  id = data[[id_col]]
  
  tibble(
    Biomarker = biomarker,
    n = sum(!is.na(x)),
    N = n_distinct(id[!is.na(x)]),
    Mean = mean(x, na.rm = TRUE),
    SD = sd(x, na.rm = TRUE),
    Median = median(x, na.rm = TRUE),
    Q1 = quantile(x, 0.25, na.rm = TRUE, names = FALSE),
    Q3 = quantile(x, 0.75, na.rm = TRUE, names = FALSE),
    Min = min(x, na.rm = TRUE),
    Max = max(x, na.rm = TRUE),
    P95 = quantile(x, 0.95, na.rm = TRUE, names = FALSE),
    Missing_pct = 100 * mean(is.na(x))
  )
}

# Substitua "id_p" e "value" pelos nomes reais das colunas
stats_anxiety = desc_biomarker(
  epileptic.l, "anxiety", "id", "anxiety"
)

stats_depress = desc_biomarker(
  epileptic.l, "depress", "id", "depress"
)

rbind(stats_anxiety,
      stats_depress)



## Fitting the survival and joint models

ctrl2 = list(
  int.strategy = "ccd", ## "eb" for velocity and ccd to generate posterior samples
  cfg = TRUE,
  cmin = 0,
  
  priorFixed = list(mean = 0, prec = 0.01,
                    mean.intercept = 0, prec.intercept = 0.01),
  priorAssoc = list(mean = 0, prec = 0.01),
  
  priorSRE_ind = list(mean = 0, prec = 0.01),
  
  # Wishart for var-cov of random effects.
  priorRandom = list(r = 5, R = 1), 
  
  cpo = T,dic = T, waic = T,
  #seed = 123,
  baselineHyper = list(
      alpha = list(prior = "normal", param = c(0, 0.01), initial = -1)
    )
)



set.seed(127)
survpart = joint(formSurv = inla.surv(with.time, with.status2) ~ trt, 
                         control=ctrl2,
                         basRisk = "dgompertzsurv", 
                         id = "id", 
                         dataSurv = epileptic.s)

summary(survpart, sdcor=T)

survpart$summary.fixed
survpart$summary.hyperpar## defective model

# survpart$cpo$cpo
# survpart$cpo$pit


#hist(survpart$cpo$cpo)


## [Fitting two longitudinal biomarkers]
## (2L) anxiety & depress 

formsurv = inla.surv(with.time, with.status2) ~ trt
form_anxiety = anxiety ~ time + trt + (1+time|id)
form_depress = depress ~ time + trt + (1+time|id)
#form_aep = aep ~ time + trt + (1+time|id)

head(epileptic.l)

set.seed(123)

set.seed(123)
joint_dgompertz = joint(formSurv = formsurv,
                         basRisk = "dgompertzsurv",
                         id = "id", dataSurv = epileptic.s,
                         formLong=list(form_anxiety, form_depress),
                         family=c("poisson","poisson"), 
                         dataLong=list(epileptic.l,epileptic.l),
                         timeVar="time", assoc=c("SRE_ind","SRE_ind"),
                         corLong=T,run=T, control=ctrl2)

set.seed(123)
joint_dgompertz5 = joint(formSurv = formsurv,
                        basRisk = "dgompertzsurv",
                        id = "id", dataSurv = epileptic.s,
                        formLong=list(form_anxiety, form_depress),
                        family=c("nbinomial","nbinomial"), 
                        dataLong=list(epileptic.l,epileptic.l),
                        timeVar="time", assoc=c("SRE_ind","SRE_ind"),
                        corLong=T,run=T, control=ctrl2)



summary(joint_dgompertz)
summary(joint_dgompertz5)


round(joint_dgompertz$summary.hyperpar, 3)
round(joint_dgompertz$summary.fixed, 3)

round(joint_dgompertz5$summary.hyperpar, 3)
round(joint_dgompert5$summary.fixed, 3)


INLA::inla.priors.used(joint_dgompertz5)
 
joint_dgompertz$waic$waic
joint_dgompertz5$waic$waic

joint_dgompertz$dic$dic
joint_dgompertz5$dic$dic




# hist(joint_dgompertz$cpo$pit)
# hist(joint_dgompertz5$cpo$pit)

round(joint_dgompertz5$summary.hyperpar, 3)
round(joint_dgompertz5$summary.fixed, 3)



## Generating posterior samples from the joint model with  negative binomial distribution for biomarkers
set.seed(123)
ns=1000
samples_post = inla.posterior.sample(ns, joint_dgompertz5, seed = 123)


names(joint_dgompertz5)
joint_dgompertz5$names.fixed
joint_dgompertz5$misc$configs$contents
joint_dgompertz5$marginals.hyperpar$`alpha parameter for dGompertz-surv[4]`

# Defining valid IDs and extracting number of samples
ids_valid1 = 1:nrow(epileptic.s)
ids_valid2 = (nrow(epileptic.s)+1):(2*nrow(epileptic.s))
ids_valid3 = (2*nrow(epileptic.s)+1):(3*nrow(epileptic.s))
ids_valid4 = (3*nrow(epileptic.s)+1):(4*nrow(epileptic.s))


# 1. Extracting scalar parameters across all samples (resulting in vectors of length ns)
alpha_s = sapply(samples_post, function(x) x$hyperpar["alpha parameter for dGompertz-surv[3]"])
mu_s = exp(sapply(samples_post, function(x) x$latent["Intercept_S1:1", 1]))
psi_pm_s = sapply(samples_post, function(x) x$latent["trtLTG_S1:1", 1])

## Random effects (latent)
call_int_L1 = paste0("IDIntercept_L1:", ids_valid1)
call_slp_L1 = paste0("IDtime_L2:", ids_valid2)
call_int_L2 = paste0("IDIntercept_L2:", ids_valid3)
call_slp_L2 = paste0("IDtime_L2:", ids_valid4)


int_L1_s = sapply(samples_post, function(x) x$latent[call_int_L1, 1])
slp_L1_s = sapply(samples_post, function(x) x$latent[call_slp_L1, 1])
int_L2_s = sapply(samples_post, function(x) x$latent[call_int_L2, 1])
slp_L2_s = sapply(samples_post, function(x) x$latent[call_slp_L2, 1])



## gammas for survival
gamma_b01_s = sapply(samples_post, function(x) x$hyperpar["Beta for SRE_Intercept_L1_S1"])
gamma_b11_s = sapply(samples_post, function(x) x$hyperpar["Beta for SRE_time_L1_S1"])
gamma_b02_s = sapply(samples_post, function(x) x$hyperpar["Beta for SRE_Intercept_L2_S1"])
gamma_b12_s = sapply(samples_post, function(x) x$hyperpar["Beta for SRE_time_L2_S1"])


sre_int_L1_weighted = t(t(int_L1_s) * gamma_b01_s)
sre_slp_L1_weighted = t(t(slp_L1_s) * gamma_b11_s)
sre_int_L2_weighted = t(t(int_L2_s) * gamma_b02_s)
sre_slp_L2_weighted = t(t(slp_L2_s) * gamma_b12_s)



# Covariate matrix W
W = as.matrix((epileptic.s$trt=="LTG")*1)

# Matrix multiplication: W (n x 1) %*% psi_mat (1 x 1000) 
W_eff = W %*% psi_pm_s

inner_sum = W_eff + 
  sre_int_L1_weighted + sre_slp_L1_weighted + 
  sre_int_L2_weighted + sre_slp_L2_weighted 

ratio_s = mu_s / alpha_s

inner_exp_scaled = t( t(exp(inner_sum)) * ratio_s )

p_cures_t = exp(inner_exp_scaled)

# Transposing to achieve the final requested dimension (1000 rows x 915 columns)
p_cures = t(p_cures_t)

dim(p_cures) 


#p_est=colMeans(p_cures)
p_est=apply(p_cures, MARGIN = 2, FUN = mean)
p_pop_est=mean(p_est)
p_quantiles=quantile(p_est, probs = c(0.025,0.975))

p_est_LTG=mean(apply(p_cures[,W==1], MARGIN = 2, FUN = mean))
p_est_CBZ=mean(apply(p_cures[,W==0], MARGIN = 2, FUN = mean))


 

## posterior samples from the joint model with  poisson distribution for biomarkers
set.seed(123)
ns=1000
samples_post1 = inla.posterior.sample(ns, joint_dgompertz, seed = 123)


names(joint_dgompertz)
joint_dgompertz$names.fixed
joint_dgompertz$misc$configs$contents
joint_dgompertz$marginals.hyperpar$`alpha parameter for dGompertz-surv[3]`


# 1. Extracting scalar parameters across all samples (resulting in vectors of length ns)
alpha_s1 = sapply(samples_post1, function(x) x$hyperpar["alpha parameter for dGompertz-surv[3]"])
mu_s1 = exp(sapply(samples_post1, function(x) x$latent["Intercept_S1:1", 1]))
psi_pm_s1 = sapply(samples_post1, function(x) x$latent["trtLTG_S1:1", 1])

## Random effects (latent)
int_L1_s1 = sapply(samples_post1, function(x) x$latent[call_int_L1, 1])
slp_L1_s1 = sapply(samples_post1, function(x) x$latent[call_slp_L1, 1])
int_L2_s1 = sapply(samples_post1, function(x) x$latent[call_int_L2, 1])
slp_L2_s1 = sapply(samples_post1, function(x) x$latent[call_slp_L2, 1])



## gammas for survival
gamma_b01_s1 = sapply(samples_post1, function(x) x$hyperpar["Beta for SRE_Intercept_L1_S1"])
gamma_b11_s1 = sapply(samples_post1, function(x) x$hyperpar["Beta for SRE_time_L1_S1"])
gamma_b02_s1 = sapply(samples_post1, function(x) x$hyperpar["Beta for SRE_Intercept_L2_S1"])
gamma_b12_s1 = sapply(samples_post1, function(x) x$hyperpar["Beta for SRE_time_L2_S1"])


sre_int_L1_weighted1 = t(t(int_L1_s1) * gamma_b01_s1)
sre_slp_L1_weighted1 = t(t(slp_L1_s1) * gamma_b11_s1)
sre_int_L2_weighted1 = t(t(int_L2_s1) * gamma_b02_s1)
sre_slp_L2_weighted1 = t(t(slp_L2_s1) * gamma_b12_s1)


# Matrix multiplication: W (n x 1) %*% psi_mat (1 x 1000) 
#W_eff = W %*% psi_mat
W_eff1 = W %*% psi_pm_s1

inner_sum1 = W_eff1 + 
        sre_int_L1_weighted1 + sre_slp_L1_weighted1 + 
        sre_int_L2_weighted1 + sre_slp_L2_weighted1

ratio_s1 = mu_s1 / alpha_s1

inner_exp_scaled1 = t(t(exp(inner_sum1)) * ratio_s1)

p_cures_t1 = exp(inner_exp_scaled1)

# Transposing to achieve the final requested dimension (1000 rows x 915 columns)
p_cures1 = t(p_cures_t1)


#p_est=colMeans(p_cures)
p_est1=apply(p_cures1, MARGIN = 2, FUN = mean)
p_pop_est1=median(p_est1)
p_quantiles1=quantile(p_est1, probs = c(0.025,0.975))

p_est_LTG1=mean(apply(p_cures1[,W==1], MARGIN = 2, FUN = mean))
p_est_CBZ1=mean(apply(p_cures1[,W==0], MARGIN = 2, FUN = mean))



## Kaplan-Meier vs graph
kmfit = survfit(Surv(with.time, with.status2) ~ trt, data = epileptic.s)


p = ggsurvplot(
  fit            = kmfit,
  data           = epileptic.s,
  palette        = c("black", "black"),   
  size           = 1.2,                   
  xlim           = c(0, 9),
  ylim           = c(0, 1),
  break.y.by     = 0.25,
  xlab           = "Time (Years)",
  ylab           = "Survival probability",
  
  legend         = "none",                
  legend.title   = "Treatment",           
  
  ggtheme        = theme_bw(base_size = 20),
  risk.table     = TRUE,                
  risk.table.col = "strata",
  legend.labs    = c("CBZ", "LTG"),     
  
  risk.table.height = 0.3, 
  fontsize          = 8      
)


p$plot = p$plot +
  geom_hline(yintercept = p_est_CBZ1, color = "darkred", linetype = "dashed", linewidth = 1.2) +
  geom_hline(yintercept = p_est_LTG1, color = "steelblue", linetype = "dashed", linewidth = 1.2) +
  
  annotate("text", x = 7, y = p_est_CBZ1 + 0.05, 
           label = paste0("hat(p)[CBZ] == '", sprintf("%.3f", p_est_CBZ1), "'"), 
           parse = TRUE, color = "darkred", size = 7, hjust = 0) +
  
  annotate("text", x = 7, y = p_est_LTG1 + 0.05, 
           label = paste0("hat(p)[LTG] == '", sprintf("%.3f", p_est_LTG1), "'"), 
           parse = TRUE, color = "steelblue", size = 7, hjust = 0) +
  
  theme(
    strip.text       = element_text(size = 22),
    panel.grid.minor = element_blank(),
    
    axis.title.x = element_text(size = 24, margin = margin(t = 15)),
    axis.title.y = element_text(size = 24, margin = margin(r = 15)),
    
    axis.text.x = element_text(size = 24, color = "black"),
    axis.text.y = element_text(size = 24, color = "black"),
    
    axis.ticks = element_line(linewidth = 1.5, color = "black"),
    axis.ticks.length = grid::unit(0.25, "cm")
  )


p$table = p$table + 
  theme(
    axis.title.y = element_text(size = 20, margin = margin(r = 15)), 
    axis.text.x  = element_text(size = 20, color = "black"),         
    axis.text.y  = element_text(size = 20, color = "black"),         
    plot.title   = element_text(size = 22),
    

    axis.ticks        = element_line(linewidth = 1, color = "black"), 
    axis.ticks.length = unit(0.25, "cm")                              
  )

print(p)

p


################ qunatifying the impact of the treatment in the individual hazard

exp(psi_pm_s1)

mean(exp(psi_pm_s1))

quantile(exp(psi_pm_s1), probs = c(0.025, 0.975))

################ Negative Binomial

p2 = ggsurvplot(
  fit            = kmfit,
  data           = epileptic.s,
  palette        = c("black", "black"),   
  size           = 1.2,                   
  xlim           = c(0, 9),
  ylim           = c(0, 1),
  break.y.by     = 0.25,
  xlab           = "Time (Years)",
  ylab           = "Survival probability",
  
  legend         = "none",                
  legend.title   = "Treatment",           
  
  ggtheme        = theme_bw(base_size = 20),
  risk.table     = TRUE,                
  risk.table.col = "strata",
  legend.labs    = c("CBZ", "LTG"),     
  
  risk.table.height = 0.3, 
  
  fontsize          = 8      
)


p2$plot = p2$plot +
  geom_hline(yintercept = p_est_CBZ, color = "darkgreen", linetype = "dashed", linewidth = 1.2) +
  geom_hline(yintercept = p_est_LTG, color = "purple4", linetype = "dashed", linewidth = 1.2) +
  
  annotate("text", x = 7, y = p_est_CBZ + 0.05, 
           label = paste0("hat(p)[CBZ] == '", sprintf("%.3f", p_est_CBZ), "'"), 
           parse = TRUE, color = "darkgreen", size = 7, hjust = 0) +
  
  annotate("text", x = 7, y = p_est_LTG + 0.05, 
           label = paste0("hat(p)[LTG] == '", sprintf("%.3f", p_est_LTG), "'"), 
           parse = TRUE, color = "purple4", size = 7, hjust = 0) +
  
  theme(
    strip.text       = element_text(size = 22),
    panel.grid.minor = element_blank(),
    
    axis.title.x = element_text(size = 24, margin = margin(t = 15)),
    axis.title.y = element_text(size = 24, margin = margin(r = 15)),
    
    axis.text.x = element_text(size = 24, color = "black"),
    axis.text.y = element_text(size = 24, color = "black"),
    
    axis.ticks = element_line(linewidth = 1.5, color = "black"),
    axis.ticks.length = grid::unit(0.25, "cm")
  )


p2$table = p2$table + 
  theme(
    axis.title.y = element_text(size = 20, margin = margin(r = 15)), 
    axis.text.x  = element_text(size = 20, color = "black"),         
    axis.text.y  = element_text(size = 20, color = "black"),         
    plot.title   = element_text(size = 22),
    
    axis.ticks        = element_line(linewidth = 1, color = "black"),
    axis.ticks.length = unit(0.25, "cm")                             
  )

print(p2)

p
p2



## Ploting the individual cures -----------------

idx_group1 = W[,1]==1
color_group = "steelblue"  #"darkred"

# idx_group1 = W[,1]==0
# color_group = "darkred"  #"darkred"


post_pest_group1 = apply(p_cures1[,idx_group1], MARGIN=2, FUN = mean)
ci_pest_group1 = t(apply(p_cures1[,idx_group1], MARGIN=2, FUN = quantile, probs = c(0.025,0.975)))

df_group1 = data.frame(
  Subject = 1:sum(idx_group1),              
  Probability = post_pest_group1,       
  Lower = ci_pest_group1[,1],   
  Upper = ci_pest_group1[,2]    
)

group_prob_cure = mean(df_group1$Probability)
group_prob_cure_ci = quantile(df_group1$Probability, probs=c(0.025,0.975))

label_text = paste0(
  "hat(p) == '", sprintf("%.3f", group_prob_cure), "'"
)

ggplot(df_group1, aes(x = Subject, y = Probability)) +
  
  geom_linerange(aes(ymin = Lower, ymax = Upper), color = "grey80", linewidth = 1) +
  

  geom_point(color = "black", size = 2.5) +
  
  geom_hline(yintercept = group_prob_cure, color = color_group, linewidth = 1.4) +

  ylim(c(0,1)) +
  
  annotate(
    "text", 
    x = Inf, y = Inf, 
    label = label_text, 
    parse = TRUE, 
    hjust = 1.05, vjust = 1.5, 
    color = color_group, 
    size = 8 
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



## Ploting the dgompertz and gompertz -----------------

ctrl3 = list(
  int.strategy = "eb", ## "eb" for velocity and ccd to generate posterior samples
  cfg = TRUE,
  cmin = 0,
  
  priorFixed = list(mean = 0, prec = 0.01,
                    mean.intercept = 0, prec.intercept = 0.01),
  priorAssoc = list(mean = 0, prec = 0.01),
  
  priorSRE_ind = list(mean = 0, prec = 0.01),
  
  # Wishart for var-cov of random effects.
  priorRandom = list(r = 3, R = 1), 
  
  cpo = T,dic = T, waic = T,
  
  baselineHyper = list(
    alpha = list(prior = "loggamma", param = c(0.01, 0.01), initial = 0.1)
  )
)

set.seed(123)
joint_gompertz2 = joint(formSurv = formsurv,
                         basRisk = "gompertzsurv",
                         id = "id", dataSurv = epileptic.s,
                         formLong=list(form_anxiety, form_depress),
                         family=c("poisson","poisson"), 
                         dataLong=list(epileptic.l,epileptic.l),
                         timeVar="time", assoc=c("SRE_ind","SRE_ind"),
                         corLong=T,run=T, control=ctrl_inlajoint)



summary(joint_gompertz2, sdcor=T)
summary(joint_dgompertz, sdcor=T)


joint_dgompertz$waic$waic
joint_gompertz2$waic$waic

joint_dgompertz$dic$dic
joint_gompertz2$waic$waic


## comparing the the survival curves using the defective behavior


## Sampling the model non-defective Gompertz

ns=1000
samples_post_gomp = inla.posterior.sample(ns, joint_gompertz2, seed = 123)


names(joint_gompertz2)
joint_gompertz2$names.fixed
joint_gompertz2$misc$configs$contents
joint_gompertz2$marginals.hyperpar$`alpha parameter for Gompertz-surv[3]`


# Extracting scalar parameters across all samples (resulting in vectors of length ns)
alpha_s2 = sapply(samples_post_gomp, function(x) x$hyperpar["alpha parameter for Gompertz-surv[3]"])
mu_s2 = exp(sapply(samples_post_gomp, function(x) x$latent["Intercept_S1:1", 1]))
psi_pm_s2 = sapply(samples_post_gomp, function(x) x$latent["trtLTG_S1:1", 1])


int_L1_s2 = sapply(samples_post_gomp, function(x) x$latent[call_int_L1, 1])
slp_L1_s2 = sapply(samples_post_gomp, function(x) x$latent[call_slp_L1, 1])
int_L2_s2 = sapply(samples_post_gomp, function(x) x$latent[call_int_L2, 1])
slp_L2_s2 = sapply(samples_post_gomp, function(x) x$latent[call_slp_L2, 1])
# int_L3_s2 = sapply(samples_post_gomp, function(x) x$latent[call_int_L3, 1])
# slp_L3_s2 = sapply(samples_post_gomp, function(x) x$latent[call_slp_L3, 1])


## gammas for survival
gamma_b01_s2 = sapply(samples_post_gomp, function(x) x$hyperpar["Beta for SRE_Intercept_L1_S1"])
gamma_b11_s2 = sapply(samples_post_gomp, function(x) x$hyperpar["Beta for SRE_time_L1_S1"])
gamma_b02_s2 = sapply(samples_post_gomp, function(x) x$hyperpar["Beta for SRE_Intercept_L2_S1"])
gamma_b12_s2 = sapply(samples_post_gomp, function(x) x$hyperpar["Beta for SRE_time_L2_S1"])
# gamma_b03_s2 = sapply(samples_post_gomp, function(x) x$hyperpar["Beta for SRE_Intercept_L2_S1"])
# gamma_b13_s2 = sapply(samples_post_gomp, function(x) x$hyperpar["Beta for SRE_time_L2_S1"])

sre_int_L1_weighted2 = t(t(int_L1_s2) * gamma_b01_s2)
sre_slp_L1_weighted2 = t(t(slp_L1_s2) * gamma_b11_s2)
sre_int_L2_weighted2 = t(t(int_L2_s2) * gamma_b02_s2)
sre_slp_L2_weighted2 = t(t(slp_L2_s2) * gamma_b12_s2)
# sre_int_L3_weighted2 = t(t(int_L3_s2) * gamma_b02_s2)
# sre_slp_L3_weighted2 = t(t(slp_L3_s2) * gamma_b12_s2)


# Matrix multiplication: W (n x 1) %*% psi_mat (1 x 1000) 
#W_eff = W %*% psi_mat
W_eff2 = W %*% psi_pm_s2

inner_sum2 = W_eff2 + 
  sre_int_L1_weighted2 + sre_slp_L1_weighted2 + 
  sre_int_L2_weighted2 + sre_slp_L2_weighted2 
  # sre_int_L3_weighted2 + sre_slp_L3_weighted2


ratio_s2 = mu_s2 / alpha_s2

inner_exp_scaled2 = t( t(exp(inner_sum2)) * ratio_s2)



t_seq = seq(0, 10, length.out = 200)

n_ind  = nrow(inner_sum)   
n_samp = ncol(inner_sum)  
n_time = length(t_seq)     

S_array = array(NA, dim = c(n_ind, n_samp, n_time))
exp_term = t(t(exp(inner_sum)) * mu_s)

S_array2 = array(NA, dim = c(n_ind, n_samp, n_time))
exp_term2 = t(t(exp(inner_sum2)) * mu_s2)


for (k in 1:n_time) {
  t = t_seq[k]
  
  time_term = (1 / alpha_s) * (exp(alpha_s * t) - 1)
  S_array[, , k] = t(exp(-t(exp_term) * time_term))
  
  time_term2 = (1 / alpha_s2) * (exp(alpha_s2 * t) - 1)
  S_array2[, , k] = t(exp(-t(exp_term2) * time_term2))
}

S_mean_indiv = apply(S_array, c(1, 3), mean)
S_mean_pop = apply(S_array, 3, mean)

S_mean_indiv2 = apply(S_array2, c(1, 3), mean)
S_mean_pop2 = apply(S_array2, 3, mean)


# plot(kmfit2, conf.int = F, xlim = c(0, 10), lwd = 2)
# lines(t_seq, S_mean_pop, lwd = 2, col = "darkgreen")
# lines(t_seq, S_mean_pop2, lwd = 2, col = "purple4")

kmfit2 = survfit(Surv(with.time, with.status2) ~ 1, data = epileptic.s)

# 3. Data frames para as curvas adicionais
df_surv_pop = data.frame(
  time = t_seq,
  surv_defgomp = as.numeric(S_mean_pop),
  surv_gomp = as.numeric(S_mean_pop2)
)

df_pop2 = data.frame(
  time = t_seq,
  surv = as.numeric(S_mean_pop2)
)


# 2. Gráfico de Kaplan-Meier
p3 = ggsurvplot(
  fit               = kmfit2,
  data              = epileptic.s,
  palette           = "black",
  size              = 1.2,
  xlim              = c(0, 9),
  ylim              = c(0, 1),
  break.y.by        = 0.25,
  xlab              = "Time (Years)",
  ylab              = "Survival probability \n",
  legend            = "none",
  ggtheme           = theme_bw(base_size = 24),
  risk.table        = TRUE,
  risk.table.col    = "strata",
  risk.table.height = 0.3,
  risk.table.ylab   = "",        # Remove o título "strata" da tabela
  risk.table.y.text = FALSE,     # Remove os rótulos de grupo da lateral
  conf.int          = FALSE,
  fontsize          = 7
)

# Garante a remoção total do eixo Y da tabela de risco no ggplot2
p3$table = p3$table + 
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

# 4. Adicionar as curvas sobre o Kaplan-Meier
p3$plot = p3$plot +
  geom_line(
    data = df_surv_pop,
    aes(x = time, y = surv_defgomp, color = "Defective Gompertz"),
    inherit.aes = FALSE,
    linetype = "solid",
    linewidth = 1.2
  ) +
  geom_line(
    data = df_pop2,
    aes(x = time, y = surv, color = "Gompertz"),
    inherit.aes = FALSE,
    linetype = "solid",
    linewidth = 1.2
  ) +
  scale_color_manual(
    values = c(
      "Defective Gompertz" = "darkgreen",
      "Gompertz"     = "purple4"
    ),
    breaks = c("Defective Gompertz", "Gompertz"),
    na.value = "black"
  ) +
  labs(color = NULL) +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold")
  )

# Imprimir o gráfico + tabela
p3


## Prior sensivity analysis considering alpha


## Fitting the survival and joint models

ctrl2 = list(
  int.strategy = "ccd", ## "eb" for velocity and ccd to generate posterior samples
  cfg = TRUE,
  cmin = 0,
  
  priorFixed = list(mean = 0, prec = 0.01,
                    mean.intercept = 0, prec.intercept = 0.01),
  priorAssoc = list(mean = 0, prec = 0.01),
  
  priorSRE_ind = list(mean = 0, prec = 0.01),
  
  # Wishart for var-cov of random effects.
  priorRandom = list(r = 5, R = 1), 
  
  baselineHyper = list(
    alpha = list(prior = "normal", param = c(0, 0.01), initial = -1)
  )
)


ctrl3 = list(
  int.strategy = "eb", ## "eb" for velocity and ccd to generate posterior samples
  cfg = TRUE,
  cmin = 0, priorFixed = list(mean = 0, prec = 0.01,
                    mean.intercept = 0, prec.intercept = 0.01),
  priorAssoc = list(mean = 0, prec = 0.01),
  priorSRE_ind = list(mean = 0, prec = 0.01),
  # Wishart for var-cov of random effects.
  priorRandom = list(r = 5, R = 1), 
  
  baselineHyper = list(
    alpha = list(prior = "normal", param = c(0, 0.05), initial = -1)
  )
)

ctrl4 = list(
  int.strategy = "eb", ## "eb" for velocity and ccd to generate posterior samples
  cfg = TRUE,
  cmin = 0, priorFixed = list(mean = 0, prec = 0.01,
                              mean.intercept = 0, prec.intercept = 0.01),
  priorAssoc = list(mean = 0, prec = 0.01),
  priorSRE_ind = list(mean = 0, prec = 0.01),
  # Wishart for var-cov of random effects.
  priorRandom = list(r = 5, R = 1), 
  
  baselineHyper = list(
    alpha = list(prior = "normal", param = c(0, 0.5), initial = -1)
  )
)



## [Fitting two longitudinal biomarkers]
## (2L) anxiety & depress 

formsurv = inla.surv(with.time, with.status2) ~ trt
form_anxiety = anxiety ~ time + trt + (1+time|id)
form_depress = depress ~ time + trt + (1+time|id)



joint_dgompertz13 = joint(formSurv = formsurv,
                        basRisk = "dgompertzsurv",
                        id = "id", dataSurv = epileptic.s,
                        formLong=list(form_anxiety, form_depress),
                        family=c("poisson","poisson"), 
                        dataLong=list(epileptic.l,epileptic.l),
                        timeVar="time", assoc=c("SRE_ind","SRE_ind"),
                        corLong=T,run=T, control=ctrl3)

joint_dgompertz14 = joint(formSurv = formsurv,
                          basRisk = "dgompertzsurv",
                          id = "id", dataSurv = epileptic.s,
                          formLong=list(form_anxiety, form_depress),
                          family=c("poisson","poisson"), 
                          dataLong=list(epileptic.l,epileptic.l),
                          timeVar="time", assoc=c("SRE_ind","SRE_ind"),
                          corLong=T,run=T, control=ctrl4)


summary(joint_dgompertz)
summary(joint_dgompertz13)
summary(joint_dgompertz14)


joint_dgompertz$summary.hyperpar[1,]
joint_dgompertz13$summary.hyperpar[1,]
joint_dgompertz14$summary.hyperpar[1,]


joint_dgompertz$summary.hyperpar
joint_dgompertz13$summary.hyperpar
joint_dgompertz14$summary.hyperpar











