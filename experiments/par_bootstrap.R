rm(list=ls())
library(boot)

set.seed(7635)

N=500
beta0 = 1
beta1 = 0.85
beta2 = 1
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
true_cntl = mean(p0)    #true marginal probability in control arm
true_trt = mean(p1)    #true marginal probability in treatment arm
true_OR = OR.fun(true_cntl,true_trt)
true_RR = RR.fun(true_cntl,true_trt)
true_RD = RD.fun(true_cntl,true_trt)


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


m = glm(y ~ as.factor(arm)+c, data=mydata, family='binomial')
tmp = predict(m, newdata=mydata, type='response')
num_boot = 10000
boot_results = matrix(NA,nrow=num_boot,ncol=6)
for (i in 1:num_boot) {
  y_b = rbinom(n=500,1,tmp)
  data_b = data.frame(arm,c,y_b)
  colnames(data_b)[3] = 'y'
  boot_results[i,] = lrest.fun(data_b)
}
c(quantile(boot_results[,4],0.025), quantile(boot_results[,4],0.975))
c(quantile(boot_results[,5],0.025), quantile(boot_results[,5],0.975))
c(quantile(boot_results[,6],0.025), quantile(boot_results[,6],0.975))

bootstrap = boot(mydata, lrest.fun, R=num_boot)
boot.ci(bootstrap, type='perc', index=4)$percent[,c(4,5)]
boot.ci(bootstrap, type='perc', index=5)$percent[,c(4,5)]
boot.ci(bootstrap, type='perc', index=6)$percent[,c(4,5)]


