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

## 1. [Load data]
set.seed(1)
data=joineRML::epileptic.qol
data$time <- data$time/365.25
data$with.time <- data$with.time/365.25

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
kmfit=survfit(Surv(with.time, with.status2)~trt,data = epileptic.s)
plot(kmfit, xlab = "Time (Years)", ylab = "Survival probability",
     main = "", ylim = c(0, 1), conf.int=F)
grid()

dev.off()

## Analyzing the characteristics of the biomarkers
## for the poisson, mean=variance
mean(epileptic.l$anxiety,na.rm=T);var(epileptic.l$anxiety,na.rm=T)
mean(epileptic.l$aep,na.rm=T);var(epileptic.l$aep,na.rm=T) # breaks the assumption
mean(epileptic.l$depress,na.rm=T);var(epileptic.l$depress,na.rm=T)

library(dplyr)
library(ggplot2)

# Calculando estatísticas intra-sujeito para o biomarcador 'aep' (como exemplo)
# Adicionando 'trt' ao agrupamento
disp_intra_trt <- epileptic.l %>%
  group_by(trt, id) %>% 
  summarise(
    n_medidas = n(),
    media_aep = mean(depress, na.rm = TRUE),
    var_aep = var(depress, na.rm = TRUE),
    .groups = "drop" # Remove o agrupamento após o sumário
  ) %>%
  filter(n_medidas > 1) %>%
  mutate(indice_dispersao = var_aep / media_aep)

# Calculando o resumo numérico do índice de dispersão por tratamento
# resumo_dispersao <- disp_intra_trt %>%
#   group_by(trt) %>%
#   summarise(
#     pacientes = n(),
#     media_indice = mean(indice_dispersao, na.rm = TRUE),
#     mediana_indice = median(indice_dispersao, na.rm = TRUE)
#   )

# print(resumo_dispersao)


disp_intra_trt_anx <- epileptic.l %>%
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
    diff_abs = abs(var_anx - media_anx) # Criando o módulo da diferença
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

disp_intra_trt_dep <- epileptic.l %>%
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
    diff_abs = abs(var_dep - media_dep) # Criando o módulo da diferença
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



## 3. [Inla estimation of the survival] 

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
  
  baselineHyper = list(
    alpha = list(prior = "normal", param = c(0, 0.01), initial = -1)
  )
)


set.seed(127)
survpart = joint(formSurv = inla.surv(with.time, with.status2) ~ trt, 
                         control=ctrl,
                         basRisk = "dgompertzsurv", 
                         id = "id", 
                         dataSurv = epileptic.s)

summary(survpart, sdcor=T)

survpart$summary.fixed
survpart$summary.hyperpar## defective model


## [Fitting two longitudinal biomarkers]
## (2L) anxiety & depress 

formsurv = inla.surv(with.time, with.status2) ~ trt
form_anxiety = anxiety ~ time + trt + (1+time|id)
form_depress = depress ~ time + trt + (1+time|id)

head(epileptic.l)

set.seed(123)
joint_dgompertz = joint(formSurv = formsurv,
                        basRisk = "dgompertzsurv",
                        id = "id", dataSurv = epileptic.s,
                        formLong=list(form_anxiety, form_depress),
                        family=c("poisson","poisson"), 
                        dataLong=list(epileptic.l,epileptic.l),
                        timeVar="time", assoc=c("SRE_ind","SRE_ind"),
                        corLong=T,run=T, control=ctrl)


sumamry(joint_dgompertz)








shared_long3 = joint(formSurv=inla.surv(with.time, with.status2) ~ trt,
                     dataSurv=epileptic.s, basRisk="weibullsurv", NbasRisk=50,
                     formLong=list(anxiety ~ time + trt + (1+time|id),
                                   depress ~ time + trt + (1+time|id)),
                     family=c("poisson","poisson"), dataLong=epileptic.l, id="id",
                     timeVar="time", assoc=c("SRE_ind","SRE_ind"), corLong=T,
                     run=F, control=list(int.strategy="eb", cfg = TRUE,
                                         config=TRUE,
                                         priorSRE_ind=list(mean=0,prec=1),
                                         priorFixed=list(mean.intercept=0,prec.intercept=0.01,
                                                         mean=0, prec=0.01)))

shared_long3$.args$family = c("poisson", "poisson", "dgompertzsurv")
shared_long3$.args$control.family[[3]] = survpart$.args$control.family
shared_long3$.args$control.inla$cmin = 0

init.time=Sys.time()
modLS3 = joint.run(shared_long3)
end.time=Sys.time()
end.time-init.time ## 22 seconds

SDC = summary(modLS3, sdcor=T)$ReffList[[1]];SDC
class(modLS3) = "inla"

## summary of the model
modLS3$summary.fixed
modLS3$summary.hyperpar
SDC


## Refitting the model to generate a sample of the posterior
refit3 = shared_long3
refit3$.args$control.inla$int.strategy = "ccd"
refit3$.args$control.compute$config = TRUE

init.time2=Sys.time()
refit3 = inla.rerun(refit3)
end.time2=Sys.time()

end.time2-init.time2 # 1 hour to refit using the ccd strategy

length(refit3$misc$configs)
length(refit3$misc$configs$config)
#str(refit3$misc$configs$config[[1]])

refit3$summary.fixed
refit3$summary.hyperpar

refit3$size.random

names(refit3)

#save(refit3, file = "refit3.RData")

set.seed(123)
ns=1000
samples_post = inla.posterior.sample(ns, refit3)

names(refit3$summary.random) ## parametros random
## example of the 20-th sample
samples_post[[20]]$hyperpar
samples_post[[20]]$latent

# Create matrizes to kepp the values of the samples
id_trat = which(epileptic.s$trt=="LTG")
id_placebo = which(epileptic.s$trt!="LTG")
id_int_L1 = 4250:(4249+2176)-1   
id_time_L1 =(4249+2176):(4249+2176*2-1)

id_int_L2 = (4249+2176*2):(4249+2176*3-1)  
id_time_L2 = (4249+2176*3):(4249+2176*4-1) 

id_SRE_int_L1 = (4249+2176*4):(4249+2176*5-1)
id_SRE_time_L1 = (4249+2176*5):(4249+2176*6-1)

id_SRE_int_L2 = (4249+2176*6):(4249+2176*7-1)
id_SRE_time_L2 = (4249+2176*7):(4249+2176*8-1)



matp_trat = matrix(data=NA,nrow=length(id_trat),ncol=ns)
matp_placebo = matrix(data=NA,nrow=length(id_placebo),ncol=ns)
vals_alpha = numeric(ns)
vals_sgamma = numeric(ns)

ori_int = 1:nrow(epileptic.s)
ori_time = (nrow(epileptic.s)+1):(2*nrow(epileptic.s))

cop_int = (2*nrow(epileptic.s)+1):(3*nrow(epileptic.s))
cop_time =  (3*nrow(epileptic.s)+1):(4*nrow(epileptic.s))

for(i in 1:ns){
  s_alpha = samples_post[[i]]$hyperpar["alpha parameter for dGompertz-surv[3]"]
  s_gamma0 = tail(samples_post[[i]]$latent,8)[1]
  s_psi = tail(samples_post[[i]]$latent,8)[2]
  
  sre1_int = samples_post[[i]]$latent[id_int_L2][ori_int]*samples_post[[i]]$hyperpar["Beta for SRE_Intercept_L1_S1"]
  sre1_time = samples_post[[i]]$latent[id_int_L2][ori_time]*samples_post[[i]]$hyperpar["Beta for SRE_time_L1_S1"]
  
  sre2_int = samples_post[[i]]$latent[id_int_L2][cop_int]*samples_post[[i]]$hyperpar["Beta for SRE_Intercept_L2_S1"]
  sre2_time = samples_post[[i]]$latent[id_int_L2][cop_time]*samples_post[[i]]$hyperpar["Beta for SRE_time_L2_S1"]
  
  s_p_trat = exp((1/s_alpha)*exp(s_gamma0+s_psi)*exp(sre1_int[id_trat]+sre1_time[id_trat]+sre2_int[id_trat]+sre2_time[id_trat]))
  s_p_placebo = exp((1/s_alpha)*exp(s_gamma0)*exp(sre1_int[id_placebo]+sre1_time[id_placebo]+sre2_int[id_placebo]+sre2_time[id_placebo]))
  
  matp_trat[,i] = s_p_trat
  matp_placebo[,i] = s_p_placebo
  vals_alpha[i] = s_alpha
  vals_sgamma[i] = s_gamma0
}

# 
# sd(samples_post[[s_sample]]$latent[id_int_L2][ori_int])
# sd(samples_post[[s_sample]]$latent[id_int_L2][cop_int])
# 
# sd(samples_post[[s_sample]]$latent[id_int_L1][cop_int])
# sd(samples_post[[s_sample]]$latent[id_int_L1][ori_int])
# 
# sd(samples_post[[s_sample]]$latent[id_time_L2][ori_time])
# sd(samples_post[[s_sample]]$latent[id_time_L2][cop_time])
# 
# sd(samples_post[[s_sample]]$latent[id_time_L1][cop_time])
# sd(samples_post[[s_sample]]$latent[id_time_L1][ori_time])
# 

# cor(samples_post[[s_sample]]$latent[id_int_L2][ori_int],
#     samples_post[[s_sample]]$latent[id_int_L2][cop_int])
# 
# cor(samples_post[[s_sample]]$latent[id_time_L2][ori_time],
#     samples_post[[s_sample]]$latent[id_time_L2][cop_time])
# 
# 
# cor(samples_post[[s_sample]]$latent[id_int_L2][ori_int],
#     samples_post[[s_sample]]$latent[id_time_L2][cop_time])

## The point estimation (posterior means) of individual cures
rowMeans(matp_trat)
rowMeans(matp_placebo)

mean(rowMeans(matp_trat))
mean(rowMeans(matp_placebo))

## c.i para as medias
quantile(rowMeans(matp_trat), probs = c(0.025,0.975))
quantile(rowMeans(matp_placebo), probs = c(0.025,0.975))

# LTG: [0.3102605 ; 0.7594128]
# CBZ: [0.2078553 ; 0.6694665]

## media dos c.i's
colMeans(t(apply(matp_trat, 1, quantile, probs = c(0.025, 0.975))))
colMeans(t(apply(matp_placebo, 1, quantile, probs = c(0.025, 0.975))))

# LTG: [0.4498585 ; 0.6772675]
# CBZ: [0.3288015 ; 0.5698826]

## [Plot of average probability of cure with credibility intervals]
library(ggplot2)

data_cure =  data.frame(
  id = 1:length(id_trat),
  cure_mean = rowMeans(matp_trat),
  cure_lower = t(apply(matp_trat, 1, quantile, probs = c(0.025, 0.975)))[,1],
  cure_upper = t(apply(matp_trat, 1, quantile, probs = c(0.025, 0.975)))[,2]
)

mean_trat = mean(data_cure$cure_mean)

ggplot(data_cure, aes(x = id, y = cure_mean)) +
  geom_errorbar(aes(ymin = cure_lower, ymax = cure_upper), color = "grey", width = 0.2) +
  geom_point(size = 1.1) +
  geom_hline(yintercept = mean_trat, color = "steelblue", linewidth = 0.6) +  # ajuste o valor conforme o seu gráfico
  labs(x = "Subject", y = "Probability of cure") +
  theme_bw()


plot_idcrue = function(mat,colhline){
  credibility = t(apply(mat, 1, quantile, probs = c(0.025, 0.975)))
  
  data_cure =  data.frame(
    id = 1:nrow(mat),
    cure_mean = rowMeans(mat),
    cure_lower = credibility [,1],
    cure_upper = credibility [,2]
  )
  
  mean_cure_mean = mean(data_cure$cure_mean)
    
  a = ggplot(data_cure, aes(x = id, y = cure_mean)) +
      geom_errorbar(aes(ymin = cure_lower, ymax = cure_upper), color = "grey50", width = 0.2) +
      geom_point(size = 1.1) +
      geom_hline(yintercept = mean_cure_mean, color = colhline, linewidth = 0.6) + 
      labs(x = "Subject", y = "Probability of cure") +
      coord_cartesian(ylim = c(0, 1)) +
      theme_bw()
  
  return(a)
}

plot_idcrue(mat=matp_trat, colhline="steelblue")
plot_idcrue(mat=matp_placebo, colhline="darkred")

## [Plot of the survival curve with the average cure fraction per group]
library(survminer)
library(survival)

survplot = ggsurvplot(kmfit,
           data = epileptic.s,
           conf.int = FALSE,         # Intervalo de confiança
           risk.table = FALSE,       # Tabela de risco abaixo do gráfico
           xlab = "Time (Years)", 
           ylab = "Survival Probability",
           palette = c("grey10", "grey10"),
           legend = "none",
           ggtheme = theme_bw(),
           size=0.8)

round(mean(rowMeans(matp_trat)),3)
round(mean(rowMeans(matp_placebo)),3)

survplot$plot + 
  geom_hline(yintercept = mean(rowMeans(matp_trat)), color = "steelblue", linewidth = 0.6) +
  geom_hline(yintercept = mean(rowMeans(matp_placebo)), color = "darkred", linewidth = 0.6) +
  annotate("text", x = 5.45, y = mean(rowMeans(matp_trat)), label = expression(hat(p)[LTG] == "0.570"),
           hjust = 0, vjust = -0.5, color = "steelblue", size = 4) +
  annotate("text", x = 5.45, y = mean(rowMeans(matp_placebo)), label = expression(hat(p)[CBZ] == "0.453"),
           hjust = 0, vjust = -0.5, color = "darkred", size = 4) 
  
