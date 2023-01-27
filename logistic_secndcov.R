
# Study the effect of covariate adjustment on type I error rate and power
# when model 1 already has one covariate; model 2 adds a correlated covariate
# model 2 is correctly specified

rm(list=ls())
setwd('/Users/liyunfan/Documents/Simulation_logistic')

library(MASS)


N=500
set.seed(10004)
#Repeat M times
M=100000

#Scenario: Second covariate in the model
#a: w/o treatment effect
# beta0 = 1
beta1 = 0
# beta2 = 1
beta3 = 1
#Save results
m1_est = rep(NA, M); m2_est = rep(NA, M)
m1_std = rep(NA, M); m2_std = rep(NA, M)
m1_z = rep(NA, M); m2_z = rep(NA, M)
m1_reject = rep(NA, M); m2_reject = rep(NA, M)
p0_var = rep(NA, M)
p0_mean = rep(NA, M)
p0omitc_var = rep(NA,M)
p0omitc_mean = rep(NA,M)

for (i in 1:M) {
  #Generate covariate
  tmp = mvrnorm(n=N, mu=c(0,0), Sigma=matrix(c(2.25,1,1,2.25),ncol=2))
  c = tmp[,1]        #c and u are bivariate normal distributed
  u = tmp[,2]        #u is already in the analysis model; c will be added
  #Treatment arm
  arm = c(rep(0,N/2),rep(1,N/2))
  #Generate probabilities
  logit = beta1*arm+c+beta3*u        #data generating model contains beta3 and u
  p = exp(logit)/(1+exp(logit))
  logit0 = c+beta3*u
  p0 = exp(logit0)/(1+exp(logit0))
  logit0_omitc = beta3*u
  p0_omitc = exp(logit0_omitc)/(1+exp(logit0_omitc))
  #Generate responses
  y = rbinom(N,1,p)
  mydata = data.frame(arm,c,u,y)
  
  #Fit unadjusted model
  m1 = glm(y ~ as.factor(arm)+u, data=mydata, family='binomial')    #c omitted
  m1_est[i] = summary(m1)$coefficients[2,1]
  m1_std[i] = summary(m1)$coefficients[2,2]
  m1_z[i] = summary(m1)$coefficients[2,3]
  m1_reject[i] = summary(m1)$coefficients[2,4]<0.05
  #Fit adjusted model
  m2 = glm(y ~ as.factor(arm)+u+c, data=mydata, family='binomial')
  m2_est[i] = summary(m2)$coefficients[2,1]
  m2_std[i] = summary(m2)$coefficients[2,2]
  m2_z[i] = summary(m2)$coefficients[2,3]
  m2_reject[i] = summary(m2)$coefficients[2,4]<0.05
  p0_var[i] = var(p0)
  p0_mean[i] = mean(p0)
  p0omitc_var[i] = var(p0_omitc)
  p0omitc_mean[i] = mean(p0_omitc)
}


#Type I error
mean(m1_z>qnorm(0.975) | m1_z<qnorm(0.025))
t1err_m1 = sum(m1_reject)/M
mean(m2_z>qnorm(0.975) | m2_z<qnorm(0.025))
t1err_m2 = sum(m2_reject)/M


########################
#b: w/ treatment effect
# beta0 = 1
beta1 = 0.75
# beta2 = 1
beta3 = 1
#Save results
m1_est = rep(NA, M); m2_est = rep(NA, M)
m1_std = rep(NA, M); m2_std = rep(NA, M)
m1_z = rep(NA, M); m2_z = rep(NA, M)
m1_reject = rep(NA, M); m2_reject = rep(NA, M)
p0_var = rep(NA, M)
p0_mean = rep(NA, M)
p0omitc_var = rep(NA,M)
p0omitc_mean = rep(NA,M)

for (i in 1:M) {
  #Generate covariate
  tmp = mvrnorm(n=N, mu=c(0,0), Sigma=matrix(c(2.25,1,1,2.25),ncol=2))
  c = tmp[,1]        #c and u are bivariate normal distributed
  u = tmp[,2]        #u is already in the analysis model
  #Treatment arm
  arm = c(rep(0,N/2),rep(1,N/2))
  #Generate probabilities
  logit = beta1*arm+c+beta3*u
  p = exp(logit)/(1+exp(logit))
  logit0 = c+beta3*u
  p0 = exp(logit0)/(1+exp(logit0))
  logit0_omitc = beta3*u
  p0_omitc = exp(logit0_omitc)/(1+exp(logit0_omitc))
  #Generate responses
  y = rbinom(N,1,p)
  mydata = data.frame(arm,c,u,y)
  
  #Fit unadjusted model
  m1 = glm(y ~ as.factor(arm)+u, data=mydata, family='binomial')
  m1_est[i] = summary(m1)$coefficients[2,1]
  m1_std[i] = summary(m1)$coefficients[2,2]
  m1_z[i] = summary(m1)$coefficients[2,3]
  m1_reject[i] = summary(m1)$coefficients[2,4]<0.05
  #Fit adjusted model
  m2 = glm(y ~ as.factor(arm)+c+u, data=mydata, family='binomial')
  m2_est[i] = summary(m2)$coefficients[2,1]
  m2_std[i] = summary(m2)$coefficients[2,2]
  m2_z[i] = summary(m2)$coefficients[2,3]
  m2_reject[i] = summary(m2)$coefficients[2,4]<0.05
  p0_var[i] = var(p0)
  p0_mean[i] = mean(p0)
  p0omitc_var[i] = var(p0_omitc)
  p0omitc_mean[i] = mean(p0_omitc)
}


#Variance and 'bias' are inflated by the same factor
bias_factor = mean(m1_est/m2_est)
bias_factor_std = sqrt(var(m1_est/m2_est))
var_factor = mean(m1_std^2/m2_std^2)
var_factor_std = sqrt(var(m1_std^2/m2_std^2))
#The theoretical efficiency factor
eff_fact_m2_unadj = sqrt(1-p0_var/(p0_mean*(1-p0_mean)))    #this factor compares m2 to unadjusted
eff_fact_m1_unadj = sqrt(1-p0omitc_var/(p0omitc_mean*(1-p0omitc_mean)))    #this factor compares m1 to unadjusted
eff_factor = mean(eff_fact_m2_unadj/eff_fact_m1_unadj)    #now look at m2 compared to m1
eff_factor_std = sd(eff_fact_m2_unadj/eff_fact_m1_unadj)
#The actual efficiency factor (ratio of z scores)
aeff_factor = mean(m1_z/m2_z)
aeff_factor_std = sqrt(var(m1_z/m2_z))

#Power
mean(m1_z>qnorm(0.975) | m1_z<qnorm(0.025))
power_m1 = sum(m1_reject)/M
mean(m2_z>qnorm(0.975) | m2_z<qnorm(0.025))
power_m2 = sum(m2_reject)/M

cbind(round(t1err_m1,3)*100, round(t1err_m2,3)*100)
cbind(round(power_m1,3)*100, round(power_m2,3)*100,
      round(bias_factor,2), round(var_factor,2),
      round(mean(p0_mean),2), round(mean(p0_var),2), round(eff_factor,2), round(aeff_factor,2))
cbind(bias_factor_std, var_factor_std, sqrt(var(p0_mean)), eff_factor_std, aeff_factor_std)



