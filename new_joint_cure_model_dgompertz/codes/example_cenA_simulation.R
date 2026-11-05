## example fitting data to gompertz distribution
rm(list = ls())
library(INLA)
library(INLAjoint)
library(MASS)


## [Function to generate data from defective Gompertz data]
gen.cure.gompertz = function(n,a,b,p=0) {
  if (length(b) == 1) b <- rep(b, n)
  if (length(p) == 1) p <- rep(p, n)
  t=rep(NA,n)
  
  for(i in 1:n){
    t[i]=log((-(a/b[i])*log(1-runif(n=1,min=0,max=1-p[i]))) + 1)*(1/a)
  }
  
  t_finite = ifelse(t==Inf,0,t)
  u2 = runif(n=n,0,20)
  t2 = pmin(t,u2)
  delta = ifelse(t<u2,1,0)
  return(as.data.frame(cbind(time=t2,delta)))
}



## [Function to generate longitudinal and survival data from defective Gompertz data]
gen.joint.cure.dgompertz = function(n, alpha, gammas, betas.fixed,
                                    sd.int, sd.slope, rho,
                                    gap, followup){
  
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
  biomarker = rpois(nmesy, lambda = exp(linPred))
  idY = as.factor(rep(1:n, each=nmesindiv)) # id
  
  # longitudinal dataset
  longDat = data.frame(id=as.integer(idY), timej, biomarker, bin)#, cts)
  

  
  eta = cbind(1, b_int, b_slo, binXs) %*% gammas 
  
  mu = exp(eta)
  
  survDat = gen.cure.gompertz(n=n, a=alpha, b=mu) 
  colnames(survDat) = c("time", "delta")
  survDat = data.frame(id=1:n, survDat, binXs, b_int, b_slo)
  
  # survival time and censoring for the longitudinal dataset
  longDat = merge(longDat, survDat[, c("id", "time")], by="id")
  
  # Only longitudinal measurements which occurs before the moment of event or censoring
  longDat = longDat[longDat$timej <= longDat$time, ]
  
  # Remove the column 'time' of the longitudinal dataset 
  longDat$time = NULL
  
  # Organize individuals and time
  longDat = longDat[order(longDat$id, longDat$timej), ]
  rownames(longDat) = NULL 
  
  # ---------------------------------------------------------
  
  return(list(survdata=survDat, longdata=longDat))
}


beta0=2.5# intercept
beta1=-0.1 # slope
beta2=0.1 # binary covariate

#betas.fixed0=c(beta0, beta1, beta2, beta3)
betas.fixed0=c(beta0, beta1, beta2)
sd.int0 = 0.4
sd.slo0 = 0.4
rho0 = -0.3

alpha0=0.2 ##defective
gamma0=-0.1;gamma1=0.55;gamma2=0.5;gamma3=-0.75
gammas0 = c(gamma0,gamma1,gamma2,gamma3)

parameters = c(alpha0,gammas0,betas.fixed0,sd.int0,sd.slo0,rho0)
names(parameters)=c("alpha","gamma0","gamma1", "gamma2","psi1",
                    "beta0","beta1","beta2",
                    "sd.int","sd.slo","rho")



data_gen = gen.joint.cure.dgompertz(n=500,alpha=alpha0,gammas=gammas0,
                         betas.fixed=betas.fixed0,
                         sd.int=sd.int0,sd.slope=sd.slo0,
                         rho=rho0,
                         gap=0.5*(1/12), followup=4)

survdata = data_gen$survdata
longdata = data_gen$longdata


head(survdata)
head(longdata)

dim(longdata)
table(table(longdata$id))



library(survival)

kmfit = survfit(Surv(time, delta)~1, data=survdata)

plot(kmfit)


ctrl_inlajoint = list(
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

form_y1 =biomarker ~ timej + bin + (1 + timej | id)
formsurv=inla.surv(time, delta)~binXs  

mdjoint = joint(formSurv = formsurv,
                basRisk = "dgompertzsurv", 
                id = "id", dataSurv = survdata,
                formLong=list(form_y1),
                family=c("poisson"), dataLong=list(longdata),
                timeVar="timej", assoc=c("SRE_ind"), 
                control=ctrl_inlajoint,
                corLong=F,run=T)

summary(mdjoint, sdcor=TRUE)


## recuperando os valores dos parametros
alpha_hat = summary(mdjoint, sdcor = TRUE)$SurvEff[[1]][1,]
gamma0_hat = mdjoint$summary.fixed[1,] 
gammas_hat = summary(mdjoint, sdcor = TRUE)$AssocLS
psi_hat = summary(mdjoint, sdcor = TRUE)$SurvEff[[1]][3,]

## Fixed effects
betas_hat = summary(mdjoint, sdcor = TRUE)$FixedEff[[1]][1:3,]

## sd of random effects and correlation
sds_reff_hat = summary(mdjoint, sdcor = TRUE)$ReffList[[1]][1:2,]
rho_hat = summary(mdjoint, sdcor = TRUE)$ReffList[[1]][3,]


means_hat = c(alpha_hat[,"mean"], gamma0_hat[,"mean"], gammas_hat[,"mean"], psi_hat[,"mean"], betas_hat[,"mean"], sds_reff_hat[,"mean"], rho_hat["mean"])
names(means_hat)=names(parameters)
length(means_hat); length(parameters)


sds_hat = c(alpha_hat[,"sd"], gamma0_hat[,"sd"], gammas_hat[,"sd"], psi_hat[,"sd"], betas_hat[,"sd"], sds_reff_hat[,"sd"], rho_hat["sd"]); names(sds_hat)=names(parameters)
lw_ci_hat =  c(alpha_hat[,"0.025quant"], gamma0_hat[,"0.025quant"], gammas_hat[,"0.025quant"], psi_hat[,"0.025quant"], betas_hat[,"0.025quant"], sds_reff_hat[,"0.025quant"], rho_hat["0.025quant"]); names(lw_ci_hat)=names(parameters)
up_ci_hat =  c(alpha_hat[,"0.975quant"], gamma0_hat[,"0.975quant"], gammas_hat[,"0.975quant"], psi_hat[,"0.975quant"], betas_hat[,"0.975quant"], sds_reff_hat[,"0.975quant"], rho_hat["0.975quant"]); names(up_ci_hat)=names(parameters)

bias = (means_hat-parameters) 
squared_error = (means_hat-parameters)^2 

sds_hat
check_ic = (parameters > lw_ci_hat) & (parameters < up_ci_hat)


bias
check_ic









