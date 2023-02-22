# Goal: Does the efficiency factor work on RR and RD as well?
# Experiment: Use Jacobian to generate variance of RR and RD estimates
# Compare CI by Jacobian to CI by nonparametric bootstrapping 

rm(list=ls())
setwd('/Users/liyunfan/Documents/Simulation_logistic/trteffect-simulation-functions')
library(boot)


#Number of patients in a dataset
N=500
#Generate M datasets
M=50
#Within each dataset, use B bootstrap samples
num_boot=1000

#This function generates probabilities predicted by the unadjusted model
#Input: dataset, treatment arm assignment='arm'
#Output: predicted probabilities, in the control and treatment arm respectively
lrprd_unadj.fun = function(data) {
  m = glm(y ~ as.factor(arm), data=data, family='binomial')
  arm0 = unname(predict(m, newdata=data.frame(arm=0), type='response'))
  arm1 = unname(predict(m, newdata=data.frame(arm=1), type='response'))
  return(list('control'=arm0, 'treatment'=arm1))
}
#This function generates the marginal probabilities prodicted by the adjusted model
#Input: dataset, treatment arm assignment='arm', a SINGLE covariate='c'
#Output: marginal predicted probabilities, in the control and treatment arm respectively
lrprd_adj.fun = function(data) {
  N = dim(data)[1]
  m = glm(y ~ as.factor(arm)+c, data=data, family='binomial')
  arm0 = predict(m, newdata=data.frame(arm=rep(0,N),c=c), type='response')
  arm1 = predict(m, newdata=data.frame(arm=rep(1,N),c=c), type='response')
  return(list('control'=mean(arm0), 'treatment'=mean(arm1)))
}
#These functions generate estimands from marginal probabilities
OR.fun = function(control,treatment) {(treatment/(1-treatment))/(control/(1-control))}
RR.fun = function(control,treatment) {treatment/control}
RD.fun = function(control,treatment) {treatment-control}
#This function generates estimands from both the unadjusted and adjusted models
#For the purpose of bootstrapping
#Input: dataset, treatment arm assignment='arm', a SINGLE covariate='c'
#Output: predicted OR, RR and RD
lrest.fun = function(data, idx) {
  df = data[idx,]
  tmp_unadj = lrprd_unadj.fun(df)
  tmp_adj = lrprd_adj.fun(df)
  OR_unadj = OR.fun(tmp_unadj$control, tmp_unadj$treatment)
  RR_unadj = RR.fun(tmp_unadj$control, tmp_unadj$treatment)
  RD_unadj = RD.fun(tmp_unadj$control, tmp_unadj$treatment)
  OR_adj = OR.fun(tmp_adj$control, tmp_adj$treatment)
  RR_adj = RR.fun(tmp_adj$control, tmp_adj$treatment)
  RD_adj = RD.fun(tmp_adj$control, tmp_adj$treatment)
  return(cbind(OR_unadj,RR_unadj,RD_unadj,OR_adj,RR_adj,RD_adj))
}


########################
#Scenario 1b: Correctly specified model, w/ treatment effect
beta0 = 1
beta1 = 0.75
beta2 = 1
#Save the results
est_original = matrix(NA,nrow=M,ncol=6)
true_cntl = rep(NA,length=M); true_trt = rep(NA,length=M)
true_RR = rep(NA,length=M)
true_RD = rep(NA,length=M)
RR_unadj_ci = matrix(NA,nrow=M,ncol=2)
RD_unadj_ci = matrix(NA,nrow=M,ncol=2)
RR_adj_ci = matrix(NA,nrow=M,ncol=2)
RD_adj_ci = matrix(NA,nrow=M,ncol=2)

RD_est = rep(NA,M)
RD_est_unadj = rep(NA,M)
RD_var = rep(NA,M)
RD_var_unadj = rep(NA,M)

RR_est = rep(NA,M)
RR_est_unadj = rep(NA,M)
lgRR_var = rep(NA,M)
lgRR_var_unadj = rep(NA,M)

for (i in 1:M) {
  if (i %% (M/10) == 0) {print(paste('i =',i))}
  set.seed(5050+i)
  #Generate covariate
  c = rnorm(N,0,1.5)
  #Treatment arm
  arm = c(rep(0,N/2),rep(1,N/2))
  #Generate the dataset
  logit = beta0+beta1*arm+beta2*c
  p = exp(logit)/(1+exp(logit))
  y = rbinom(N,1,p)
  mydata = data.frame(arm,c,y)
  
  logit0 = beta0+beta2*c    #true logit if assigned to control arm
  p0 = exp(logit0)/(1+exp(logit0))
  logit1 = logit0+beta1    #true logit if assigned to treatment arm
  p1 = exp(logit1)/(1+exp(logit1))
  true_cntl[i] = mean(p0)    #true marginal probability in control arm
  true_trt[i] = mean(p1)    #true marginal probability in treatment arm
  true_RR[i] = RR.fun(true_cntl[i],true_trt[i])
  true_RD[i] = RD.fun(true_cntl[i],true_trt[i])
  
  #Fit the adjusted model
  m1 = glm(y ~ as.factor(arm), data=mydata, family='binomial')
  m2 = glm(y ~ as.factor(arm)+c, data=mydata, family='binomial')

  #Calculate the Jacobian and Wald statistic for risk difference (RD) -adjusted model
  tmp1 = exp(cbind(rep(1,N),rep(1,N),c)%*%m2$coefficients)
  tmp0 = exp(cbind(rep(1,N),rep(0,N),c)%*%m2$coefficients)
  denom1_m2 = (1+tmp1)^2    #this denominator will be used a lot
  denom0_m2 = (1+tmp0)^2    #beta1 doesn't appear in this denominator
  Jm2_RD_beta0 = mean(tmp1/denom1_m2)-mean(tmp0/denom0_m2)
  Jm2_RD_beta1 = mean(tmp1/denom1_m2)
  Jm2_RD_beta2 = mean(c*tmp1/denom1_m2)-mean(c*tmp0/denom0_m2)
  Jm2_RD = c(Jm2_RD_beta0, Jm2_RD_beta1, Jm2_RD_beta2)    #this is the Jacobian for RD
  RD_var[i] = Jm2_RD%*%vcov(m2)%*%c(Jm2_RD)    #this is the variance of RD
  RD_est[i] = mean(tmp1/(1+tmp1))-mean(tmp0/(1+tmp0))    #this is the estimate of RD

  #Calculate the Jacobian and Wald statistic for log relative risk (RR) -adjusted model
  Jm2_RR_beta0 = mean(tmp1/denom1_m2)/mean(tmp1/(1+tmp1))-mean(tmp0/denom0_m2)/(mean(tmp0/(1+tmp0)))
  Jm2_RR_beta1 = mean(tmp1/denom1_m2)/mean(tmp1/(1+tmp1))
  Jm2_RR_beta2 = mean(c*tmp1/denom1_m2)/mean(tmp1/(1+tmp1))-mean(c*tmp0/denom0_m2)/mean(tmp0/(1+tmp0))
  Jm2_RR = c(Jm2_RR_beta0, Jm2_RR_beta1, Jm2_RR_beta2)    #this is the Jacobian for RR
  lgRR_var[i] = Jm2_RR%*%vcov(m2)%*%c(Jm2_RR)    #this is the variance of RR
  RR_est[i] = exp(log(mean(tmp1/(1+tmp1)))-log(mean(tmp0/(1+tmp0))))    #this is the estimate of RR

  #Calculate the Jacobian and Wald statistic for risk difference (RD) -unadjusted model
  tmp1 = exp(sum(m1$coefficients))
  tmp0 = exp(m1$coefficients[1])
  denom1_m1 = (1+tmp1)^2    #this denominator will be used a lot
  denom0_m1 = (1+tmp0)^2    #beta1 doesn't appear in this denominator
  Jm1_RD_beta0 = tmp1/denom1_m1-tmp0/denom0_m1
  Jm1_RD_beta1 = tmp1/denom1_m1
  Jm1_RD = c(Jm1_RD_beta0, Jm1_RD_beta1)    #this is the Jacobian for RD
  RD_var_unadj[i] = Jm1_RD%*%vcov(m1)%*%c(Jm1_RD)    #this is the variance of RD
  RD_est_unadj[i] = mean(tmp1/(1+tmp1))-mean(tmp0/(1+tmp0))    #this is the estimate of RD

  #Calculate the Jacobian and Wald statistic for log risk difference (RR) -unadjusted model
  Jm1_RR_beta0 = (tmp1/denom1_m1)/(tmp1/(1+tmp1))-(tmp0/denom0_m1)/(tmp0/(1+tmp0))
  Jm1_RR_beta1 = (tmp1/denom1_m1)/(tmp1/(1+tmp1))
  Jm1_RR = c(Jm1_RR_beta0, Jm1_RR_beta1)    #this is the Jacobian for log RR
  lgRR_var_unadj[i] = Jm1_RR%*%vcov(m1)%*%c(Jm1_RR)    #this is the variance of log RR
  RR_est_unadj[i] = exp(log(mean(tmp1/(1+tmp1)))-log(mean(tmp0/(1+tmp0))))    #this is the estimate of RD

  #Bootstrap the data and get estimates
  set.seed(9988+i)
  boot.out = boot(mydata, lrest.fun, R=num_boot)
  est_original[i,] = boot.out$t0    #these are the estimates using the entire sample
  RR_unadj_ci[i,] = quantile(boot.out$t[,2], c(.025,.975))
  RD_unadj_ci[i,] = quantile(boot.out$t[,3], c(.025,.975))
  RR_adj_ci[i,] = quantile(boot.out$t[,5], c(.025,.975))
  RD_adj_ci[i,] = quantile(boot.out$t[,6], c(.025,.975))
}

#Diagnostics: Does RD perform well?
#Compare the Wald statistics
mean(RD_est*RD_var^(-1/2))    #average Wald statistic of RD, adjusted model
mean(RD_est_unadj*RD_var_unadj^(-1/2))    #average Wald statistic of RD, unadjusted model

#Compare estimated RD to the truth
mean(RD_est)
mean(RD_est_unadj)
summary(true_RD)

#Compare Wald statistic generated CI of R to bootstrap generated CI
mean(RD_est-1.96*sqrt(RD_var))
mean(RD_adj_ci[,1])
mean(RD_est+1.96*sqrt(RD_var))
mean(RD_adj_ci[,2])

mean(RD_est_unadj-1.96*sqrt(RD_var_unadj))
mean(RD_unadj_ci[,1])
mean(RD_est_unadj+1.96*sqrt(RD_var_unadj))
mean(RD_unadj_ci[,2])

#Diagnostics: Does RR perform well?
#Compare the Wald statistics
mean(log(RR_est)*lgRR_var^(-1/2))    #average Wald statistic of log RR, adjusted model
mean(log(RR_est_unadj)*lgRR_var_unadj^(-1/2))    #average Wald statistic of log RR, unadjusted model

#Compare estimated RD to the truth
mean(RR_est)
mean(RR_est_unadj)
summary(true_RR)

#Compare Wald statistic generated CI of RD to bootstrap generated CI
mean(exp(log(RR_est)-1.96*sqrt(lgRR_var)))
mean(RR_adj_ci[,1])
mean(exp(log(RR_est)+1.96*sqrt(lgRR_var)))
mean(RR_adj_ci[,2])

mean(exp(log(RR_est_unadj)-1.96*sqrt(lgRR_var_unadj)))
mean(RR_unadj_ci[,1])
mean(exp(log(RR_est_unadj)+1.96*sqrt(lgRR_var_unadj)))
mean(RR_unadj_ci[,2])




