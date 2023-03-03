
# Study the treatment effect estimation by unadjusted and adjusted LR models
# Use g-computation to get marginal treatment effect
# a larger treatment effect; model is correctly specified

rm(list=ls())
setwd('/Users/liyunfan/Documents/Simulation_logistic')
library(boot)

load('/Users/liyunfan/Documents/Simulation_logistic/trteffect-simulation-results/logistic_trt_largeff.RData')

#Number of patients in a dataset
N=500
#Generate M datasets
M=500
#Within each dataset, use B bootstrap samples
num_boot=5000

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
#Scenario c: Correctly specified model, w/ a larger treatment effect
beta0 = 1
beta1 = 0.85
beta2 = 1
#Save the results
est_original = matrix(NA,nrow=M,ncol=6)
true_cntl = rep(NA,length=M); true_trt = rep(NA,length=M)
true_OR = rep(NA,length=M)
true_RR = rep(NA,length=M)
true_RD = rep(NA,length=M)
OR_unadj_ci = matrix(NA,nrow=M,ncol=2)
RR_unadj_ci = matrix(NA,nrow=M,ncol=2)
RD_unadj_ci = matrix(NA,nrow=M,ncol=2)
OR_adj_ci = matrix(NA,nrow=M,ncol=2)
RR_adj_ci = matrix(NA,nrow=M,ncol=2)
RD_adj_ci = matrix(NA,nrow=M,ncol=2)
p0_var = rep(NA, M)
p0_mean = rep(NA, M)


for (i in 1:M) {
  if (i %% (M/10) == 0) {print(paste('i =',i))}
  set.seed(7635+i)
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
  logit1 = beta0+beta1+beta2*c    #true logit if assigned to treatment arm
  p1 = exp(logit1)/(1+exp(logit1))
  true_cntl[i] = mean(p0)    #true marginal probability in control arm
  true_trt[i] = mean(p1)    #true marginal probability in treatment arm
  true_OR[i] = OR.fun(true_cntl[i],true_trt[i])
  true_RR[i] = RR.fun(true_cntl[i],true_trt[i])
  true_RD[i] = RD.fun(true_cntl[i],true_trt[i])
  
  #Bootstrap the data and get estimates
  set.seed(5475+i)
  boot.out = boot(mydata, lrest.fun, R=num_boot)
  est_original[i,] = boot.out$t0    #these are the estimates using the entire sample
  #Now get CIs from bootstrapping samples
  OR_unadj_ci[i,] = quantile(boot.out$t[,1], c(.025,.975))
  RR_unadj_ci[i,] = quantile(boot.out$t[,2], c(.025,.975))
  RD_unadj_ci[i,] = quantile(boot.out$t[,3], c(.025,.975))
  OR_adj_ci[i,] = quantile(boot.out$t[,4], c(.025,.975))
  RR_adj_ci[i,] = quantile(boot.out$t[,5], c(.025,.975))
  RD_adj_ci[i,] = quantile(boot.out$t[,6], c(.025,.975))
  
  #Mean and variance of p0 to calculate the expected efficiency factor
  p0_var[i] = var(p0)
  p0_mean[i] = mean(p0)
}

#Summary of true estimands
round(c(mean(true_OR),sd(true_OR)),2)
round(c(mean(true_RR),sd(true_RR)),2)
round(c(mean(true_RD),sd(true_RD)),3)
# #Summary of original estimates
# colnames(est_original) = c('OR_unadj','RR_unadj','RD_unadj','OR_adj','RR_adj','RD_adj')
# round(apply(est_original,2,FUN=quantile,probs=0.25),3)
# round(apply(est_original,2,FUN=median),3)
# round(apply(est_original,2,FUN=quantile,probs=0.75),3)
#Deviance from the truth
true_OR_m = mean(true_OR)
true_RR_m = mean(true_RR)
true_RD_m = mean(true_RD)
round(c(mean(true_OR_m-est_original[,1]),sd(true_OR_m-est_original[,1])),3)    #unadjusted OR
round(c(mean(true_RR_m-est_original[,2]),sd(true_RR_m-est_original[,2])),3)    #unadjusted RR
round(c(mean(true_RD_m-est_original[,3]),sd(true_RD_m-est_original[,3])),4)    #unadjusted RD
round(c(mean(true_OR_m-est_original[,4]),sd(true_OR_m-est_original[,4])),3)    #adjusted OR
round(c(mean(true_RR_m-est_original[,5]),sd(true_RR_m-est_original[,5])),3)    #adjusted RR
round(c(mean(true_RD_m-est_original[,6]),sd(true_RD_m-est_original[,6])),4)    #adjusted RD
#Average length of CIs with st.d.
round(c(mean(OR_unadj_ci[,2]-OR_unadj_ci[,1]), sd(OR_unadj_ci[,2]-OR_unadj_ci[,1])),2)
round(c(mean(RR_unadj_ci[,2]-RR_unadj_ci[,1]), sd(RR_unadj_ci[,2]-RR_unadj_ci[,1])),2)
round(c(mean(RD_unadj_ci[,2]-RD_unadj_ci[,1]), sd(RD_unadj_ci[,2]-RD_unadj_ci[,1])),3)
round(c(mean(OR_adj_ci[,2]-OR_adj_ci[,1]), sd(OR_adj_ci[,2]-OR_adj_ci[,1])),2)
round(c(mean(RR_adj_ci[,2]-RR_adj_ci[,1]), sd(RR_adj_ci[,2]-RR_adj_ci[,1])),2)
round(c(mean(RD_adj_ci[,2]-RD_adj_ci[,1]), sd(RD_adj_ci[,2]-RD_adj_ci[,1])),3)
#Ratio of CIs
round(mean((OR_adj_ci[,2]-OR_adj_ci[,1])/(OR_unadj_ci[,2]-OR_unadj_ci[,1])),2)
round(sd((OR_adj_ci[,2]-OR_adj_ci[,1])/(OR_unadj_ci[,2]-OR_unadj_ci[,1])),2)
round(mean((RR_adj_ci[,2]-RR_adj_ci[,1])/(RR_unadj_ci[,2]-RR_unadj_ci[,1])),2)
round(sd((RR_adj_ci[,2]-RR_adj_ci[,1])/(RR_unadj_ci[,2]-RR_unadj_ci[,1])),2)
round(mean((RD_adj_ci[,2]-RD_adj_ci[,1])/(RD_unadj_ci[,2]-RD_unadj_ci[,1])),2)
round(sd((RD_adj_ci[,2]-RD_adj_ci[,1])/(RD_unadj_ci[,2]-RD_unadj_ci[,1])),2)

#The expected efficiency factor
round(mean(sqrt(1-p0_var/(p0_mean*(1-p0_mean)))),2)
round(sqrt(var(sqrt(1-p0_var/(p0_mean*(1-p0_mean))))),2)

save(true_OR, true_RR, true_RD, est_original, OR_unadj_ci, RR_unadj_ci, RD_unadj_ci,
     OR_adj_ci, RR_adj_ci, RD_adj_ci, p0_var, p0_mean, file = "logistic_trt_largeff.RData")

