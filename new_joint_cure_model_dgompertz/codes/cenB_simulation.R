rm(list = ls())
library(INLA)
library(INLAjoint)
library(MASS)

## 1. Function to generate data from defective Gompertz data
gen.cure.gompertz = function(n, a, b, p, max_cen) { # max_cen parameter 
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

gen.joint.cure.dgompertz = function(n, alpha, gammas, sd.fr, betas.fixed, size.rbinom,
                                    sd.int, sd.slope, rho,
                                    gap, followup, max_cen){
  
  nmesindiv = followup/gap + 1 # number of individual measurement
  cov = sd.int * sd.slope * rho
  Sigma = matrix(c(sd.int^2, cov,
                   cov, sd.slope^2), ncol=2, nrow=2)
  
  MVnorm = mvrnorm(n = n, mu = rep(0,2), Sigma = Sigma)
  
  b_int = MVnorm[,1]
  b_intY = rep(b_int, each=nmesindiv) # b0
  b_slo = MVnorm[,2]
  b_sloY = rep(b_slo, each=nmesindiv) # b1
  
  binXs = rbinom(n, 1, 0.85) # binary covariate 1
  
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
  biomarker = rnbinom(nmesy, size = size.rbinom, mu = exp(linPred))
  #biomarker =pmax(1, pmin(biomarker, 40))
  
  idY = as.factor(rep(1:n, each=nmesindiv)) # id
  
  # longitudinal dataset
  longDat = data.frame(id=as.integer(idY), timej, biomarker, bin)#, cts)
  
  # survival submodel / frailty
  eta = cbind(1, b_int, b_slo, binXs) %*% gammas 
  
  mu = exp(eta)
  
  # Probability of cure (must be negative for defective Gompertz)
  p_par = exp(mu/alpha) 
  
  survDat = gen.cure.gompertz(n=n, a=alpha, b=mu, p=p_par, max_cen=max_cen) 
  colnames(survDat) = c("time", "delta")
  survDat = data.frame(id=1:n, survDat, binXs, b_int, b_slo, p_par)
  
  longDat = merge(longDat, survDat[, c("id", "time")], by="id")
  
  longDat = longDat[longDat$timej <= longDat$time, ]
  longDat$time = NULL
  
  longDat = longDat[order(longDat$id, longDat$timej), ]
  rownames(longDat) = NULL
  
  # ---------------------------------------------------------
  
  return(list(survdata=survDat, longdata=longDat))
}



## biomarkers
beta0=2.5# intercept
beta1=-0.2 # slope
beta2=0.25 # binary covariate

#betas.fixed0=c(beta0, beta1, beta2, beta3)
betas.fixed0=c(beta0, beta1, beta2)
sd.int0 = 0.5
sd.slo0 = 0.5
rho0 = -0.1

alpha0=-0.4 ##defective
gamma0=-0.1;gamma1=0.35;gamma2=0.35;gamma3=-0.5
gammas0 = c(gamma0,gamma1,gamma2,gamma3)
size.rbinom0 = 100

parameters = c(alpha0, gammas0, betas.fixed0, sd.int0, sd.slo0, rho0, size.rbinom0)
names(parameters)=c("alpha","gamma0","gamma1", "gamma2","psi1",
                    "beta0","beta1","beta2",
                    "sd.int","sd.slo","rho", "size.rbinom0")



nsamples = c(100,300,500,1000)
niter = 1000


res.matrix_poisson = array(data=NA, dim = c(niter, length(parameters[-12])*4, length(nsamples)),
                   dimnames=list(1:niter, c(rep("mean",length(parameters[-12])), rep("sd",length(parameters[-12])),
                                            rep("bias",length(parameters[-12])), rep("cp",length(parameters[-12]))),nsamples))

res.matrix_nb = array(data=NA, dim = c(niter, length(parameters)*4, length(nsamples)),
                           dimnames=list(1:niter, c(rep("mean",length(parameters)), rep("sd",length(parameters)),
                                                    rep("bias",length(parameters)), rep("cp",length(parameters))),nsamples))

censoring_values = array(data=NA, dim = c(niter, 1, length(nsamples)),
                            dimnames=list(1:niter, "censoring",nsamples))



ctrl_inlajoint_negbinom1 = list(
  int.strategy = "eb", ## ccd to generate posterior samples, eb if you want velocity
  cfg = TRUE,
  cmin = 0,
  
  priorFixed = list(mean = 0, prec = 0.01,
                    mean.intercept = 0, prec.intercept = 0.01),
  priorAssoc = list(mean = 0, prec = 0.01),
  
  priorSRE_ind = list(mean = 0, prec = 0.01),
  
  # r (2 random effects + 1) for the degress of freedom
  priorRandom = list(r = 3, R = 1),
  
  # A higher uncertainty beacause this parameter palys a big role of identifying cure fraction
  baselineHyper = list(
    alpha = list(prior = "normal", param = c(0, 0.01), initial = 0)
  )
)

save_points = seq(50,1000,by=50)


for(j in 1:length(nsamples)){
  for (i in 1:niter) {

    print(paste("Monte Carlo realization: ", i, ", sample size:", nsamples[j]))
    
    ## Generate data
    data_gen = gen.joint.cure.dgompertz(n=nsamples[j],alpha=alpha0,gammas=gammas0,
                             betas.fixed=betas.fixed0,
                             sd.int=sd.int0,sd.slope=sd.slo0,
                             rho=rho0,size.rbinom=size.rbinom0,
                             max_cen=40,
                             gap=(1/12), followup=2.5)
    
    survdata = data_gen$survdata
    longdata = data_gen$longdata
    
    # Fitting the joint model
    form_y1 = biomarker ~ timej + bin + (1 + timej | id)
    formsurv = inla.surv(time, delta)~binXs  
    
    mdjoint_poisson = joint(formSurv = formsurv,
                            basRisk = "dgompertzsurv", 
                            id = "id", dataSurv = survdata,
                            formLong=list(form_y1),
                            family=c("poisson"), dataLong=list(longdata),
                            timeVar="timej", assoc=c("SRE_ind"), 
                            control=ctrl_inlajoint_negbinom1,
                            corLong=F,run=T)
    
    mdjoint_negbinom = joint(formSurv = formsurv,
                             basRisk = "dgompertzsurv", 
                             id = "id", dataSurv = survdata,
                             formLong=list(form_y1),
                             family=c("nbinomial"), dataLong=list(longdata),
                             timeVar="timej", assoc=c("SRE_ind"), 
                             control=ctrl_inlajoint_negbinom1,
                             corLong=F,run=T)

  
    ### Extracting values from the jint model with Poisson distribution
    alpha_hat = summary(mdjoint_poisson, sdcor = TRUE)$SurvEff[[1]][1,]
    gamma0_hat = mdjoint_poisson$summary.fixed[1,] 
    gammas_hat = summary(mdjoint_poisson, sdcor = TRUE)$AssocLS
    psi_hat = summary(mdjoint_poisson, sdcor = TRUE)$SurvEff[[1]][3,]
    
    ## Fixed effects
    betas_hat = summary(mdjoint_poisson, sdcor = TRUE)$FixedEff[[1]][1:3,]
    
    ## sd of random effects and correlation
    sds_reff_hat = summary(mdjoint_poisson, sdcor = TRUE)$ReffList[[1]][1:2,]
    rho_hat = summary(mdjoint_poisson, sdcor = TRUE)$ReffList[[1]][3,]
    
    #parameters
    means_hat = c(alpha_hat[,"mean"], gamma0_hat[,"mean"], gammas_hat[,"mean"], psi_hat[,"mean"], betas_hat[,"mean"], sds_reff_hat[,"mean"], rho_hat["mean"])
    
    names(means_hat)=names(parameters[-12])
    length(means_hat); length(parameters[-12])
    
    sds_hat = c(alpha_hat[,"sd"], gamma0_hat[,"sd"], gammas_hat[,"sd"], psi_hat[,"sd"], betas_hat[,"sd"], sds_reff_hat[,"sd"], rho_hat["sd"]); names(sds_hat)=names(parameters[-12])
    lw_ci_hat =  c(alpha_hat[,"0.025quant"], gamma0_hat[,"0.025quant"], gammas_hat[,"0.025quant"], psi_hat[,"0.025quant"], betas_hat[,"0.025quant"], sds_reff_hat[,"0.025quant"], rho_hat["0.025quant"]); names(lw_ci_hat)=names(parameters[-12])
    up_ci_hat =  c(alpha_hat[,"0.975quant"], gamma0_hat[,"0.975quant"], gammas_hat[,"0.975quant"], psi_hat[,"0.975quant"], betas_hat[,"0.975quant"], sds_reff_hat[,"0.975quant"], rho_hat["0.975quant"]); names(up_ci_hat)=names(parameters[-12])
    
    
    bias = (means_hat-parameters[-12]) 
    squared_error = (means_hat-parameters[-12])^2 
    check_ic = (parameters[-12] > lw_ci_hat) & (parameters[-12] < up_ci_hat)
    
    
    
    ### Extracting values from the jint model with Negative Binomial distribution
    alpha_hat_nb = summary(mdjoint_negbinom, sdcor = TRUE)$SurvEff[[1]][1,]
    gamma0_hat_nb = mdjoint_negbinom$summary.fixed[1,] 
    gammas_hat_nb = summary(mdjoint_negbinom, sdcor = TRUE)$AssocLS
    psi_hat_nb = summary(mdjoint_negbinom, sdcor = TRUE)$SurvEff[[1]][3,]
    
    ## Fixed effects
    betas_hat_nb = summary(mdjoint_negbinom, sdcor = TRUE)$FixedEff[[1]][2:4,]
    size_hat_nb = summary(mdjoint_negbinom, sdcor = TRUE)$FixedEff[[1]][1,]
    
    ## sd of random effects and correlation
    sds_reff_hat_nb = summary(mdjoint_negbinom, sdcor = TRUE)$ReffList[[1]][1:2,]
    rho_hat_nb = summary(mdjoint_negbinom, sdcor = TRUE)$ReffList[[1]][3,]
    
    #parameters
    means_hat_nb = c(alpha_hat_nb[,"mean"], gamma0_hat_nb[,"mean"], gammas_hat_nb[,"mean"], 
                     psi_hat_nb[,"mean"], betas_hat_nb[,"mean"], sds_reff_hat_nb[,"mean"],
                     rho_hat_nb["mean"], size_hat_nb["mean"])
    
    means_hat_nb=as.numeric(means_hat_nb)
    
    names(means_hat_nb)=names(parameters)
    length(means_hat_nb); length(parameters)
    
    sds_hat_nb = c(alpha_hat_nb[,"sd"], gamma0_hat_nb[,"sd"], 
                   gammas_hat_nb[,"sd"], psi_hat_nb[,"sd"],
                   betas_hat_nb[,"sd"], sds_reff_hat_nb[,"sd"], 
                   rho_hat_nb["sd"], size_hat_nb["sd"]); 
    
    sds_hat_nb = as.numeric(sds_hat_nb)
    names(sds_hat_nb)=names(parameters)
    
    lw_ci_hat_nb =  c(alpha_hat_nb[,"0.025quant"], gamma0_hat_nb[,"0.025quant"], gammas_hat_nb[,"0.025quant"], 
                      psi_hat_nb[,"0.025quant"], betas_hat_nb[,"0.025quant"], sds_reff_hat_nb[,"0.025quant"], 
                      rho_hat_nb["0.025quant"], size_hat_nb["0.025quant"]); names(lw_ci_hat_nb)=names(parameters)
    
    lw_ci_hat_nb=as.numeric(lw_ci_hat_nb)
    names(lw_ci_hat_nb)=names(parameters)
    
    up_ci_hat_nb =  c(alpha_hat_nb[,"0.975quant"], gamma0_hat_nb[,"0.975quant"], gammas_hat_nb[,"0.975quant"], 
                      psi_hat_nb[,"0.975quant"], betas_hat_nb[,"0.975quant"], sds_reff_hat_nb[,"0.975quant"], 
                      rho_hat_nb["0.975quant"], size_hat_nb["0.975quant"]); names(up_ci_hat_nb)=names(parameters)
    
    up_ci_hat_nb=as.numeric(up_ci_hat_nb)
    names(up_ci_hat_nb)=names(parameters)
    
    bias_nb = (as.numeric(means_hat_nb)-parameters) 
    squared_error_nb = (as.numeric(means_hat_nb)-parameters)^2 
    check_ic_nb = (parameters > lw_ci_hat_nb) & (parameters < up_ci_hat_nb)
    
    res_poisson = rbind(
      means_hat,
      squared_error,
      bias,
      check_ic)
    
    res_nb = rbind(
      means_hat_nb,
      squared_error_nb,
      bias_nb,
      check_ic_nb)
    
    # Saving results
    res.matrix_poisson[i, 1:11, j] = res_poisson[1,]
    res.matrix_poisson[i, 12:22, j] = res_poisson[2,]
    res.matrix_poisson[i, 23:33, j] = res_poisson[3,]
    res.matrix_poisson[i, 34:44, j] = res_poisson[4,]
    
    censoring_values[i, , j] = 1 - mean(survdata$delta)
    
    res.matrix_nb[i, 1:12, j] = res_nb[1,]
    res.matrix_nb[i, 13:24, j] = res_nb[2,]
    res.matrix_nb[i, 25:36, j] = res_nb[3,]
    res.matrix_nb[i, 37:48, j] = res_nb[4,]
    
    if (i %in% save_points) {
    write.csv2(res.matrix_poisson[,,j], file = paste("cenB_jm_poisson_", "rep_", niter, "_n_", nsamples[j], ".csv", sep = ""))
    write.csv2(res.matrix_nb[,,j], file = paste("cenB_jm_nb_", "rep_", niter, "_n_", nsamples[j], ".csv", sep = ""))
    write.csv2(censoring_values[,,j], file = paste("cenB_jm_censoring_", "rep_", niter, "_n_", nsamples[j], ".csv", sep = ""))
    }
  }
}



# res.matrix_poisson[1:34,,1]
# 
# colMeans(res.matrix_poisson[1:50,,1])
# colMeans(res.matrix_nb[1:50,,1])
# 
# colMeans(res.matrix_nb[1:500,,1])
# 
# 
# 
# 
# colMeans(res.matrix_poisson[,,2])
# 
# colMeans(res.matrix_nb[,,1])
# colMeans(res.matrix_nb[,,2])




