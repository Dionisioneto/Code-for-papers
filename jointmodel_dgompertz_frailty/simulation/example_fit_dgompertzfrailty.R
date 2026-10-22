rm(list = ls())
library(INLA)
library(INLAjoint)
library(MASS)

## 1. Function to generate data from defective Gompertz data
gen.cure.gompertz = function(n, a, b, p, max_cen) { # max_cen parameter for estabilization
  if (length(b) == 1) b <- rep(b, n)
  if (length(p) == 1) p <- rep(p, n)
  
  # rm = 0 cured, rm = 1 event
  rm = rbinom(n=n, size=1, prob=1-p)
  t = rep(NA, n)
  
  for(i in 1:n){
    t[i] = ifelse(rm[i]==0,
                  Inf,
                  log((-(a/b[i]) * log(1 - runif(n=1, min=0, max=1-p[i]))) + 1)*(1/a))
  }
  
  u2 = runif(n=n, 0, max_cen) 
  t2 = pmin(t, u2)
  delta = ifelse(t < u2, 1, 0)
  
  return(cbind(t2, delta))
}

## 2. Function to generate longitudinal and survival data from defective Gompertz data

gen.joint.cure.dgompertz = function(n, alpha, gammas, sd.fr, betas.fixed, sd.e1,
                                    sd.int, sd.slope, rho,
                                    gap, followup, max_cen){
  
  nmesindiv = followup/gap + 1 # number of individual measurement
  cov = sd.int * sd.slope * rho
  Sigma = matrix(c(sd.int^2, cov,
                   cov, sd.slope^2), ncol=2, nrow=2)
  
  MVnorm <- mvrnorm(n = n, mu = rep(0,2), Sigma = Sigma)
  
  b_int = MVnorm[,1]
  b_intY = rep(b_int, each=nmesindiv) # b0
  b_slo = MVnorm[,2]
  b_sloY = rep(b_slo, each=nmesindiv) # b1
  
  binXs = rbinom(n, 1, 0.75) # binary covariate 1

  bin = rep(rbinom(n, 1, 0.5), each=nmesindiv)
  #ctsX = round(runif(n, 0, 1), 5) # continuous covariate
  #cts = rep(ctsX, each=nmesindiv)
  mestime = seq(0, followup, gap) # measurement times
  timej = rep(mestime, n) # time column
  nmesy = nmesindiv * n # number of longi measurements
  
  # linear predictor for longitudinal
  b_0 = betas.fixed[1]; b_1 = betas.fixed[2]; b_2 = betas.fixed[3]#; b_3 = betas.fixed[4]
  
  linPred = (b_0 + b_intY) + (b_1 + b_sloY)*timej + b_2*bin #+ b_3*cts
  
  # count longitudinal outcome
  biomarker = rnorm(nmesy, mean = linPred, sd = sd.e1)
  idY = as.factor(rep(1:n, each=nmesindiv)) # id
  
  # longitudinal dataset
  longDat = data.frame(id=as.integer(idY), timej, biomarker, bin)#, cts)
  
  # survival submodel / frailty
  frailty = rnorm(n=n, mean=0, sd=sd.fr)
  
  eta = cbind(1, b_int, b_slo, binXs) %*% gammas 
  
  mu = exp(eta + frailty)
  
  # Probability of cure (must be negative for defective Gompertz)
  p_par = exp(mu/alpha) 
  
  # Passando o max_cen para a função
  survDat = gen.cure.gompertz(n=n, a=alpha, b=mu, p=p_par, max_cen=max_cen) 
  colnames(survDat) = c("time", "delta")
  survDat = data.frame(id=1:n, survDat, binXs, b_int, b_slo, p_par)
  
  
  # 1. Traz o tempo de sobrevivência/censura para o dataset longitudinal
  longDat = merge(longDat, survDat[, c("id", "time")], by="id")
  
  # 2. Mantém apenas as medidas longitudinais que ocorreram ANTES ou NO MOMENTO do evento/censura
  longDat = longDat[longDat$timej <= longDat$time, ]
  
  # 3. Remove a coluna 'time' do dataset longitudinal (mantém o dataset limpo, pois 'time' já está no survDat)
  longDat$time = NULL
  
  # 4. (Opcional, mas recomendado) Reordena os dados por indivíduo e tempo
  longDat = longDat[order(longDat$id, longDat$timej), ]
  rownames(longDat) = NULL # reseta os índices das linhas
  
  # ---------------------------------------------------------
  
  return(list(survdata=survDat, longdata=longDat))
}


## biomarkers
beta0=-0.26 # intercept
beta1=0.3  # slope
beta2=0.3 # binary covariate
beta3=0    # continuous covariate
sd.e10 = 0.58

betas.fixed0=c(beta0, beta1, beta2, beta3)
sd.int0 = 0.40
sd.slo0 = 0.56
rho0 = -0.46 

alpha0=-0.2 ##defective
gamma0=-1;gamma1=0.7;gamma2=0.5;gamma3=0.5
gammas0 = c(gamma0,gamma1,gamma2,gamma3)
sd.frailty0 = 0.5


parameters = c(alpha0,gammas0,sd.frailty0,betas.fixed0[1:3],sd.e10,sd.int0,sd.slo0,rho0)
names(parameters)=c("alpha","gamma0","gamma1", "gamma2","psi1", "sd.frailty0",
                    "beta0","beta1","beta2",
                    "sd.e10", "sd.int","sd.slo","rho")


data = gen.joint.cure.dgompertz(n=500, alpha=alpha0,
                                gammas=gammas0, sd.fr=sd.frailty0, betas.fixed=betas.fixed0,
                                sd.e1=sd.e10,
                                sd.int=sd.int0, sd.slope=sd.slo0, rho=rho0,
                                gap=(1/12), followup=2.5, max_cen=35) #gap=0.6,followup=1.8

survdata = data$survdata
longdata = data$longdata

head(survdata)
head(longdata,20)

mean(survdata$delta)
table(table(longdata$id))

library(survival)
kmfit=survfit(Surv(time,delta)~1, data = survdata)
plot(kmfit)
abline(h=mean(survdata$p_par), col="red")



ctrl_stan_like = list(
  int.strategy = "ccd", ## ccd to generate posterior samples, eb if you want velocity
  cfg = TRUE,
  cmin = 0,
  
  # prec = 0.01 (variance = 100) can be very high in the log-scale
  # prec = 0.05 guarantes high variance and low-informative prior
  priorFixed = list(mean = 0, prec = 0.01,
                    mean.intercept = 0, prec.intercept = 0.01),
  priorAssoc = list(mean = 0, prec = 0.01),
  
  priorSRE_ind = list(mean = 0, prec = 0.01),
  
  # r (2 random effects + 1) for the degress of freedom
  priorRandom = list(r = 3, R = 1), 
  
  # PC prior para a frailty. 
  # inv gamma(alpha=0.001, beta=0.001) forces the density to low values
  # the pc prior is a better choice and the parameters chosen ensures 
  # P(sd frailty > 1) = 0.05, in the log-scale, it allow for a big variation
  # in the individual risk
  
  priorFrailty = list(prior = "pc.prec", param = c(1, 0.1)),
  
  # A higher uncertainty beacause this parameter palys a big role of identifying cure fraction 
  baselineHyper = list(
    alpha = list(prior = "normal", param = c(0, 0.01), initial = -1)
  )
)


form_y1 = biomarker ~ timej + bin + (1+timej|id) 
formsurv=inla.surv(time, delta)~binXs  + (1 | id) ## frailty 

set.seed(123)
mdjoint = joint(formSurv = formsurv,
                basRisk = "dgompertzsurv", 
                id = "id", dataSurv = survdata,
                formLong=list(form_y1),
                family=c("gaussian"), dataLong=list(longdata),
                timeVar="timej", assoc=c("SRE_ind"), 
                corLong=F,run=T, control=ctrl_stan_like)





summary(mdjoint, sdcor=TRUE)

summary(mdjoint, sdcor = TRUE)$SurvEff[[1]] ## except for gamma0
summary(mdjoint, sdcor = TRUE)$FixedEff
summary(mdjoint, sdcor = TRUE)$AssocLS

## survival and association
alpha_hat = summary(mdjoint, sdcor = TRUE)$SurvEff[[1]][1,]
gamma0_hat = mdjoint$summary.fixed[1,] 
gammas_hat = summary(mdjoint, sdcor = TRUE)$AssocLS
psi_hat = summary(mdjoint, sdcor = TRUE)$SurvEff[[1]][3,]
sd.frailty_hat = summary(mdjoint, sdcor = TRUE)$ReffListS[[1]]

## Fixed effects
betas_hat = summary(mdjoint, sdcor = TRUE)$FixedEff[[1]][1:3,]
sd.e1_hat = summary(mdjoint, sdcor = TRUE)$FixedEff[[1]][4,]

## sd of random effects and correlation

sds_reff_hat = summary(mdjoint, sdcor = TRUE)$ReffList[[1]][1:2,]
rho_hat = summary(mdjoint, sdcor = TRUE)$ReffList[[1]][3,]


#parameters

# alphapha_hat, gamma0_hat, gammas_hat, psi_hat, frailty_hat, betas_hat, sd.e1_hat, sds_hat, rho_hatal

means_hat = c(alpha_hat[,"mean"], gamma0_hat[,"mean"], gammas_hat[,"mean"], psi_hat[,"mean"], sd.frailty_hat[,"mean"], betas_hat[,"mean"], sd.e1_hat[,"mean"],sds_reff_hat[,"mean"], rho_hat["mean"]); names(means_hat)=names(parameters)

length(means_hat); length(parameters)

#(means_hat-parameters) ## Bias
#(means_hat-parameters)^2 ## squared error

sds_hat = c(alpha_hat[,"sd"], gamma0_hat[,"sd"], gammas_hat[,"sd"], psi_hat[,"sd"], sd.frailty_hat[,"sd"], betas_hat[,"sd"], sd.e1_hat[,"sd"], sds_reff_hat[,"sd"], rho_hat["sd"]); names(sds_hat)=names(parameters)

lw_ci_hat =  c(alpha_hat[,"0.025quant"], gamma0_hat[,"0.025quant"], gammas_hat[,"0.025quant"], psi_hat[,"0.025quant"], sd.frailty_hat[,"0.025quant"], betas_hat[,"0.025quant"], sd.e1_hat[,"0.025quant"], sds_reff_hat[,"0.025quant"], rho_hat["0.025quant"]); names(lw_ci_hat)=names(parameters)


up_ci_hat =  c(alpha_hat[,"0.975quant"], gamma0_hat[,"0.975quant"], gammas_hat[,"0.975quant"], psi_hat[,"0.975quant"], sd.frailty_hat[,"0.975quant"], betas_hat[,"0.975quant"], sd.e1_hat[,"0.975quant"], sds_reff_hat[,"0.975quant"], rho_hat["0.975quant"]); names(up_ci_hat)=names(parameters)

#lw_ci_hat; up_ci_hat

bias = (means_hat-parameters) 
squared_error = (means_hat-parameters)^2 
sds_hat
check_ic = (parameters > lw_ci_hat) & (parameters < up_ci_hat)

bias;squared_error;sds_hat;check_ic

bias; check_ic

## Checking the value of general time
alphas_samples = inla.hyperpar.sample(n=500,mdjoint, improve.marginals = T)[,"alpha parameter for dGompertz-surv[2]"]

mean(alphas_samples)

tc_vals = c(0.5, 1, 3, 5, 10, 15, 20, 25, 30, 35, 40)

oneminus_c_samples = sapply(tc_vals, function(minusc) {
  1-exp(minusc*alphas_samples)
})


colnames(oneminus_c_samples) = paste0("c_", tc_vals)

# Posterior mean of t_c for each c
onmeminus_mean = colMeans(oneminus_c_samples)
onmeminus_sd = apply(oneminus_c_samples, MARGIN = 2, FUN = sd)
onmeminus_ci = apply(oneminus_c_samples, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))

oneminus_true = 1-exp(alpha0*tc_vals)

oneminus_bias = (colMeans(oneminus_c_samples) - oneminus_true)
oneminus_sd = apply(oneminus_c_samples, MARGIN = 2, FUN = sd)
oneminus_ci =  apply(oneminus_c_samples, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))

oneminus_check_ci = (oneminus_true >= oneminus_ci[1,]) & (oneminus_true <= oneminus_ci[2,])

oneminus_bias
oneminus_sd
oneminus_check_ci



# if i want to generate posterior samples (not only for hyperparameters)
set.seed(1234) ## Seed for R
inla.seed = 1234 ## seed for R INLA in C

ns=1000
samples_post = inla.posterior.sample(ns, mdjoint, seed = inla.seed) # For latent

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

