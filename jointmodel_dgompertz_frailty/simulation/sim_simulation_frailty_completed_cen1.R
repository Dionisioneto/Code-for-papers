rm(list = ls())
library(INLA)
library(INLAjoint)
library(MASS)

setwd("~/Desktop/PhD/01_Joint_Model_Defective_Frailty/simulation_results/sim3")

## 1. Function to generate data from defective Gompertz data
gen.cure.gompertz = function(n, a, b, p, max_cen) { # max_cen parameter for estabilization
  if (length(b) == 1) b <- rep(b, n)
  if (length(p) == 1) p <- rep(p, n)
  
  # rm = 0 significa curado, rm = 1 significa evento possível
  rm = rbinom(n=n, size=1, prob=1-p)
  t = rep(NA, n)
  
  for(i in 1:n){
    t[i] = ifelse(rm[i]==0,
                  Inf,
                  log((-(a/b[i]) * log(1 - runif(n=1, min=0, max=1-p[i]))) + 1)*(1/a))
  }
  
  # [CORRIGIDO] Geração de censura mais estável em vez de depender de max(t_finite)
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
  
  # Agora utilizando os elementos do argumento gammas
  eta = cbind(1, b_int, b_slo, binXs) %*% gammas 
  
  mu = exp(eta + frailty)
  
  # Probabilidade de cura (Lembrete: alpha deve ser estritamente negativo na dGompertz!)
  p_par = exp(mu/alpha) 
  
  # Passando o max_cen para a função
  survDat = gen.cure.gompertz(n=n, a=alpha, b=mu, p=p_par, max_cen=max_cen) 
  colnames(survDat) = c("time", "delta")
  survDat = data.frame(id=1:n, survDat, binXs, b_int, b_slo, p_par)
  
  # ---------------------------------------------------------
  # CORREÇÃO: Truncar os dados longitudinais pelo tempo de sobrevida
  # ---------------------------------------------------------
  
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

## Control for the priors
ctrl_stan_like = list(
  int.strategy = "ccd",
  cfg = TRUE,
  cmin = 0,
  
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
  priorFrailty = list(prior = "pc.prec", param = c(1, 0.01)),
  # A higher uncertainty beacause this parameter palys a big role of identifying cure fraction 
  baselineHyper = list(
    alpha = list(prior = "normal", param = c(0, 0.01), initial = -1)
  )
)

ctrl_stan_like0 = list(
  int.strategy = "ccd",
  cfg = TRUE,
  cmin = 0,
  
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
  
  # A higher uncertainty beacause this parameter palys a big role of identifying cure fraction 
  baselineHyper = list(
    alpha = list(prior = "normal", param = c(0, 0.01), initial = -1)
  )
)


#-----
## scenario 1
## biomarkers
beta0=-0.26 # intercept
beta1=0.15  # slope
beta2=0.15 # binary covariate
beta3=0    # continuous covariate
sd.e10 = 0.58

betas.fixed0=c(beta0, beta1, beta2, beta3)
sd.int0 = 0.40
sd.slo0 = 0.56
rho0 = -0.46 

## survival
alpha0=-0.1 ##defective
gamma0=-2;gamma1=0.74;gamma2=0.5;gamma3=0.76
gammas0 = c(gamma0,gamma1,gamma2,gamma3)
sd.frailty0 = 0.5
#-----

#-----
## scenario 2
## biomarkers
# beta0=-0.26 # intercept
# beta1=0.3  # slope
# beta2=0.3 # binary covariate
# beta3=0    # continuous covariate
# sd.e10 = 0.58
# 
# betas.fixed0=c(beta0, beta1, beta2, beta3)
# sd.int0 = 0.40
# sd.slo0 = 0.56
# rho0 = -0.46 
# 
# ## survival
# alpha0=-0.2 ##defective
# gamma0=-1;gamma1=0.3;gamma2=0.3;gamma3=0.5
# gammas0 = c(gamma0,gamma1,gamma2,gamma3)
# sd.frailty0 = 0.7
#-----


parameters = c(alpha0,gammas0,sd.frailty0,betas.fixed0[1:3],sd.e10,sd.int0,sd.slo0,rho0)
names(parameters)=c("alpha","gamma0","gamma1", "gamma2","psi1", "sd.frailty0",
                    "beta0","beta1","beta2",
                    "sd.e10", "sd.int","sd.slo","rho")


##  Monte Carlo studies
niter = 200
nsamples = c(1000)

## table for parameters
tc_vals = c(0.5, 1, 3, 5, 10, 15, 20, 25, 30, 35, 40)
## the order for estimates is always (alpha, gamma0, gamma1, gamma2, sd.frailty, beta0, beta1, beta3, sd.int0, sd.slope0, rho)

res.matrix1 = array(data=NA, dim = c(niter, length(parameters)*4, length(nsamples)),
                    dimnames=list(1:niter, c(rep("mean",length(parameters)), rep("sd",length(parameters)),
                                             rep("bias",length(parameters)), rep("cp",length(parameters))),nsamples))

res.matrix10 = array(data=NA, dim=c(niter,(length(parameters)-1)*4,length(nsamples)), 
                     dimnames=list(1:niter,c(rep("mean",(length(parameters)-1)),
                                             rep("sd",(length(parameters)-1)),rep("bias",(length(parameters)-1)),
                                             rep("cp",(length(parameters)-1))),nsamples))

res.matrix2 = array(data=NA, dim = c(niter, length(tc_vals)*4, length(nsamples)),
                    dimnames=list(1:niter, c(rep("mean",length(tc_vals)), rep("sd",length(tc_vals)),
                                             rep("bias",length(tc_vals)), rep("cp",length(tc_vals))),nsamples))

res.matrix20 = array(data=NA, dim = c(niter, length(tc_vals)*4, length(nsamples)),
                     dimnames=list(1:niter, c(rep("mean",length(tc_vals)), rep("sd",length(tc_vals)),
                                              rep("bias",length(tc_vals)), rep("cp",length(tc_vals))),nsamples))


res.matrix3 = lapply(nsamples, function(n) {
  array(NA_real_, dim = c(niter, n)) 
})

names(res.matrix3) = paste0("n_", nsamples)

res.matrix30 = res.matrix3

censoring_values = array(data=NA, dim = c(niter, 1, length(nsamples)),
                         dimnames=list(1:niter, "censoring",nsamples))

waic_values = array(data=NA, dim = c(niter, 2, length(nsamples)),
                    dimnames=list(1:niter, c("mdjoint", "mdjoint0"),nsamples))

save_points = seq(100,1000,100)

# Max trials
max_tentativas = 10 

for(j in 1:(length(nsamples))){
  for (i in 6:niter) {
    
    sucesso = FALSE 
    tentativas = 0
    
    while (!sucesso && tentativas < max_tentativas) {
      tentativas = tentativas + 1
      
      tryCatch({
        
        print(paste("Monte Carlo realization:", i, "| Sample size:", nsamples[j], "| Try:", tentativas))
        
        ## ====================================================================
        ## (1) Generate Data
        ## ====================================================================
        data = gen.joint.cure.dgompertz(n=nsamples[j], alpha=alpha0,
                                        gammas=gammas0, sd.fr=sd.frailty0, betas.fixed=betas.fixed0,
                                        sd.e1=sd.e10,
                                        sd.int=sd.int0, sd.slope=sd.slo0, rho=rho0,
                                        gap=(1/12), followup=2, max_cen=35)
        
        survdata = data$survdata
        longdata = data$longdata
        
        ## ====================================================================
        ## (2) Fit the joint models
        ## ====================================================================
        
        # Modelo 1 (com frailty)
        form_y1 = biomarker ~ timej + bin + (1+timej|id) 
        formsurv = inla.surv(time, delta)~binXs  + (1 | id) ## frailty 
        mdjoint = joint(formSurv = formsurv,
                        basRisk = "dgompertzsurv", 
                        id = "id", dataSurv = survdata,
                        formLong=list(form_y1),
                        family=c("gaussian"), dataLong=list(longdata),
                        timeVar="timej", assoc=c("SRE_ind"), 
                        corLong=F, run=T, control=ctrl_stan_like)
        
        # Modelo 0 (excluding frailty)
        formsurv0 = inla.surv(time, delta)~binXs  ## no frailty
        mdjoint0 = joint(formSurv = formsurv0,
                         basRisk = "dgompertzsurv",
                         id = "id", dataSurv = survdata,
                         formLong=list(form_y1),
                         family=c("gaussian"), dataLong=list(longdata),
                         timeVar="timej", assoc=c("SRE_ind"),
                         corLong=F, run=T, control=ctrl_stan_like0)
        
        ## ====================================================================
        ## (3) Get values and compute biases (Model 1)
        ## ====================================================================
        
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
        
        means_hat = c(alpha_hat[,"mean"], gamma0_hat[,"mean"], gammas_hat[,"mean"], psi_hat[,"mean"], sd.frailty_hat[,"mean"], betas_hat[,"mean"], sd.e1_hat[,"mean"], sds_reff_hat[,"mean"], rho_hat["mean"])
        
        sds_hat = c(alpha_hat[,"sd"], gamma0_hat[,"sd"], gammas_hat[,"sd"], psi_hat[,"sd"], sd.frailty_hat[,"sd"], betas_hat[,"sd"], sd.e1_hat[,"sd"], sds_reff_hat[,"sd"], rho_hat["sd"])
        
        lw_ci_hat =  c(alpha_hat[,"0.025quant"], gamma0_hat[,"0.025quant"], gammas_hat[,"0.025quant"], psi_hat[,"0.025quant"], sd.frailty_hat[,"0.025quant"], betas_hat[,"0.025quant"], sd.e1_hat[,"0.025quant"], sds_reff_hat[,"0.025quant"], rho_hat["0.025quant"])
        
        up_ci_hat =  c(alpha_hat[,"0.975quant"], gamma0_hat[,"0.975quant"], gammas_hat[,"0.975quant"], psi_hat[,"0.975quant"], sd.frailty_hat[,"0.975quant"], betas_hat[,"0.975quant"], sd.e1_hat[,"0.975quant"], sds_reff_hat[,"0.975quant"], rho_hat["0.975quant"])
        
        bias = (means_hat - parameters) 
        check_ic = (parameters > lw_ci_hat) & (parameters < up_ci_hat)
        
        names(bias)=names(parameters); names(check_ic) = names(parameters)
        
        res.matrix1[i, 1:length(parameters), j] = means_hat
        res.matrix1[i, (length(parameters) + 1):(2*length(parameters)), j] = sds_hat
        res.matrix1[i, (2*length(parameters) + 1):(3*length(parameters)), j] = bias
        res.matrix1[i, (3*length(parameters) + 1):(4*length(parameters)), j] = check_ic*1
        censoring_values[i, , j] = 1 - mean(survdata$delta)
        
        ## ====================================================================
        ## (4) Get values and compute biases (Model 0)
        ## ====================================================================
        
        alpha_hat0 = summary(mdjoint0, sdcor = TRUE)$SurvEff[[1]][1,]
        gamma0_hat0 = mdjoint0$summary.fixed[1,] 
        gammas_hat0 = summary(mdjoint0, sdcor = TRUE)$AssocLS
        psi_hat0 = summary(mdjoint0, sdcor = TRUE)$SurvEff[[1]][3,]
        
        betas_hat0 = summary(mdjoint0, sdcor = TRUE)$FixedEff[[1]][1:3,]
        sd.e1_hat0 = summary(mdjoint0, sdcor = TRUE)$FixedEff[[1]][4,]
        
        sds_reff_hat0 = summary(mdjoint0, sdcor = TRUE)$ReffList[[1]][1:2,]
        rho_hat0 = summary(mdjoint0, sdcor = TRUE)$ReffList[[1]][3,]
        
        means_hat0 = c(alpha_hat0[,"mean"], gamma0_hat0[,"mean"], gammas_hat0[,"mean"], psi_hat0[,"mean"], betas_hat0[,"mean"], sd.e1_hat0[,"mean"], sds_reff_hat0[,"mean"], rho_hat0["mean"]) 
        names(means_hat0) = names(parameters)[-6]
        
        sds_hat0 = c(alpha_hat0[,"sd"], gamma0_hat0[,"sd"], gammas_hat0[,"sd"], psi_hat0[,"sd"], betas_hat0[,"sd"], sd.e1_hat0[,"sd"], sds_reff_hat0[,"sd"], rho_hat0["sd"])
        
        lw_ci_hat0 = c(alpha_hat0[,"0.025quant"], gamma0_hat0[,"0.025quant"], gammas_hat0[,"0.025quant"], psi_hat0[,"0.025quant"], betas_hat0[,"0.025quant"], sd.e1_hat0[,"0.025quant"], sds_reff_hat0[,"0.025quant"], rho_hat0["0.025quant"])
        
        up_ci_hat0 = c(alpha_hat0[,"0.975quant"], gamma0_hat0[,"0.975quant"], gammas_hat0[,"0.975quant"], psi_hat0[,"0.975quant"], betas_hat0[,"0.975quant"], sd.e1_hat0[,"0.975quant"], sds_reff_hat0[,"0.975quant"], rho_hat0["0.975quant"])
        
        names(sds_hat0) = names(parameters)[-6]; names(lw_ci_hat0) = names(parameters)[-6]; names(up_ci_hat0) = names(parameters)[-6]; 
        
        bias0 = means_hat0 - parameters[-6]
        squared_error0 = (means_hat0 - parameters[-6])^2
        check_ic0 = (parameters[-6] > lw_ci_hat0) & (parameters[-6] < up_ci_hat0)
        
        #bias0; check_ic0
        #bias; check_ic
        
        len_p0 = length(parameters) - 1 
        res.matrix10[i, 1:len_p0, j] = means_hat0
        res.matrix10[i, (len_p0 + 1):(2*len_p0), j] = sds_hat0
        res.matrix10[i, (2*len_p0 + 1):(3*len_p0), j] = bias0
        res.matrix10[i, (3*len_p0 + 1):(4*len_p0), j] = check_ic0*1
        
        ## ====================================================================
        ## (5) Compute (1-c)
        ## ====================================================================
        
        alphas_samples = inla.hyperpar.sample(n=500, mdjoint, improve.marginals = T)[,"alpha parameter for dGompertz-surv[2]"]
        alphas_samples0 = inla.hyperpar.sample(n=500, mdjoint0, improve.marginals = T)[,"alpha parameter for dGompertz-surv[2]"]
        
        oneminus_c_samples = sapply(tc_vals, function(minusc) { 1-exp(minusc*alphas_samples) })
        oneminus_c_samples0 = sapply(tc_vals, function(minusc) { 1-exp(minusc*alphas_samples0) })
        
        colnames(oneminus_c_samples) = paste0("c_", tc_vals); colnames(oneminus_c_samples0) = paste0("c_", tc_vals)
        
        onmeminus_mean = colMeans(oneminus_c_samples)
        oneminus_sd = apply(oneminus_c_samples, MARGIN = 2, FUN = sd)
        oneminus_ci = apply(oneminus_c_samples, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))
        
        oneminus_true = 1-exp(alpha0*tc_vals)
        oneminus_bias = (onmeminus_mean - oneminus_true)
        oneminus_check_ci = (oneminus_true >= oneminus_ci[1,]) & (oneminus_true <= oneminus_ci[2,])
        
        onmeminus_mean0 = colMeans(oneminus_c_samples0)
        oneminus_sd0 = apply(oneminus_c_samples0, MARGIN = 2, FUN = sd)
        oneminus_ci0 = apply(oneminus_c_samples0, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))
        
        oneminus_bias0 = (onmeminus_mean0 - oneminus_true)
        oneminus_check_ci0 = (oneminus_true >= oneminus_ci0[1,]) & (oneminus_true <= oneminus_ci0[2,])
        
        res.matrix2[i, 1:length(tc_vals), j] = onmeminus_mean
        res.matrix2[i, (length(tc_vals) + 1):(2*length(tc_vals)),j] = oneminus_sd
        res.matrix2[i, (2*length(tc_vals) + 1):(3*length(tc_vals)),j] = oneminus_bias
        res.matrix2[i, (3*length(tc_vals) + 1):(4*length(tc_vals)),j] = oneminus_check_ci*1
        
        res.matrix20[i, 1:length(tc_vals), j] = onmeminus_mean0
        res.matrix20[i, (length(tc_vals) + 1):(2*length(tc_vals)),j] = oneminus_sd0
        res.matrix20[i, (2*length(tc_vals) + 1):(3*length(tc_vals)),j] = oneminus_bias0
        res.matrix20[i, (3*length(tc_vals) + 1):(4*length(tc_vals)),j] = oneminus_check_ci0*1
        
        ## ====================================================================
        ## (6) Computing the cure
        ## ====================================================================
        ns=1000
        
        # Para mdjoint
        samples_post = inla.posterior.sample(ns, mdjoint, seed = 1234)
        alpha_s = sapply(samples_post, function(x) x$hyperpar["alpha parameter for dGompertz-surv[2]"])
        mu_s = exp(sapply(samples_post, function(x) x$latent["Intercept_S1:1", 1]))
        psi_pm_s = sapply(samples_post, function(x) x$latent["binXs_S1:1", 1])
        
        ids_valid1 = 1:nrow(survdata)
        ids_valid2 = (nrow(survdata)+1):(2*nrow(survdata))
        
        call_int_L1 = paste0("IDIntercept_L1:", ids_valid1)
        call_slp_L1 = paste0("IDtimej_L1:", ids_valid2)
        int_L1_s = sapply(samples_post, function(x) x$latent[call_int_L1, 1])
        slp_L1_s = sapply(samples_post, function(x) x$latent[call_slp_L1, 1])
        
        frailty_s = sapply(samples_post, function(x) x$latent[paste0("IDIntercept_S1:", ids_valid1), 1])
        gamma_b01_s = sapply(samples_post, function(x) x$hyperpar["Beta for SRE_Intercept_L1_S1"])
        gamma_b11_s = sapply(samples_post, function(x) x$hyperpar["Beta for SRE_timej_L1_S1"])
        
        sre_int_L1_weighted = t(t(int_L1_s) * gamma_b01_s)
        sre_slp_L1_weighted = t(t(int_L1_s) * gamma_b11_s)
        
        W = as.matrix(survdata[,"binXs"])
        W_eff = W %*% psi_pm_s
        inner_sum = frailty_s + W_eff + sre_int_L1_weighted + sre_slp_L1_weighted 
        
        ratio_s = mu_s / alpha_s
        inner_exp_scaled = t( t(exp(inner_sum)) * ratio_s )
        p_cures = t(exp(inner_exp_scaled))
        bias_cure = (colMeans(p_cures) - survdata$p_par)
        
        # Para mdjoint0
        samples_post0 = inla.posterior.sample(ns, mdjoint0, seed = 1234) 
        alpha_s0 = sapply(samples_post0, function(x) x$hyperpar["alpha parameter for dGompertz-surv[2]"])
        mu_s0 = exp(sapply(samples_post0, function(x) x$latent["Intercept_S1:1", 1]))
        psi_pm_s0 = sapply(samples_post0, function(x) x$latent["binXs_S1:1", 1])
        
        int_L1_s0 = sapply(samples_post0, function(x) x$latent[call_int_L1, 1])
        slp_L1_s0 = sapply(samples_post0, function(x) x$latent[call_slp_L1, 1])
        gamma_b01_s0 = sapply(samples_post0, function(x) x$hyperpar["Beta for SRE_Intercept_L1_S1"])
        gamma_b11_s0 = sapply(samples_post0, function(x) x$hyperpar["Beta for SRE_timej_L1_S1"])
        
        sre_int_L1_weighted0 = t(t(int_L1_s0) * gamma_b01_s0)
        sre_slp_L1_weighted0 = t(t(int_L1_s0) * gamma_b11_s0)
        
        W_eff0 = W %*% psi_pm_s0
        inner_sum0 = W_eff0 + sre_int_L1_weighted0 + sre_slp_L1_weighted0
        
        ratio_s0 = mu_s0 / alpha_s0
        inner_exp_scaled0 = t( t(exp(inner_sum0)) * ratio_s0 )
        p_cures0 = t(exp(inner_exp_scaled0))
        bias_cure0 = (colMeans(p_cures0) - survdata$p_par)
        
        res.matrix3[[j]][i,] = bias_cure
        res.matrix30[[j]][i,] = bias_cure0
        
        waic_values[i,,j] = c(mdjoint$waic$waic, mdjoint0$waic$waic)
        
        # =====================================================================
        # FINAL DO BLOCO TRY: Se chegou aqui, nada deu erro.
        sucesso = TRUE 
        
      }, error = function(e) {
        
        
        message("\n--- ERROR ---")
        message("Iteration: ", i, " | Sample size: ", nsamples[j])
        message("Error message: ", e$message)
        
        if (tentativas < max_tentativas) {
          message("Excluding this dataset, generating a new dataset and trying again (Tentativa ", tentativas + 1, " of ", max_tentativas, ")...\n")
        } else {
          message("Maximum of trials reached. Skiping this iteration.\n")
        }
        
      }) 
    } 
    
    
    ## ====================================================================
    ## (7) Save results (Save the progress, regardless the error in the iteration)
    ## ====================================================================
    if (i %in% save_points) {
      write.csv2(res.matrix1[,,j], file = paste("cen_1_jmres_", "rep_", niter, "_n_", nsamples[j], ".csv", sep = ""))
      write.csv2(res.matrix2[,,j], file = paste("cen_1_jmres_oneminus_c", "rep_", niter, "_n_", nsamples[j], ".csv", sep = ""))
      write.csv2(censoring_values[,,j], file = paste("cen_1_10jmcens_", "rep_", niter, "_n_", nsamples[j], ".csv", sep = ""))
      
      write.csv2(res.matrix10[,,j], file = paste("cen_1_0jmres_", "rep_", niter, "_n_", nsamples[j], ".csv", sep = ""))
      write.csv2(res.matrix20[,,j], file = paste("cen_1_0jmres_oneminus_c", "rep_", niter, "_n_", nsamples[j], ".csv", sep = ""))
      
      write.csv2(res.matrix3[[j]], file = paste("cen_1_3jmres_", "rep_", niter, "_n_", nsamples[j], ".csv", sep = ""))
      write.csv2(res.matrix30[[j]], file = paste("cen_1_30jmres_oneminus_c", "rep_", niter, "_n_", nsamples[j], ".csv", sep = ""))
    }
    
  } 
} 


# colMeans(res.matrix1[10:20,,1])
# colMeans(res.matrix1[21:30,,1])
# colMeans(res.matrix1[31:40,,1])
# colMeans(res.matrix1[41:50,,1])
# colMeans(res.matrix1[51:60,,1])
# colMeans(res.matrix1[61:70,,1])
# 
# colMeans(res.matrix1[51:80,,1])
# 
# colMeans(res.matrix1[17:80,,1])
# 
# colMeans(res.matrix1[1:80,,1])
# 
# 
# res.matrix1[1:51,,1]
# colMeans(res.matrix1[1:51,,1])
# 
# 





















































