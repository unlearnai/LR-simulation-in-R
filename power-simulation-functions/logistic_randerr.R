
# Study the effect of covariate adjustment on type I error rate and power
# model 1 is unadjusted
# model 2 uses observed covariates with random errors


rm(list=ls())
setwd('/Users/liyunfan/Documents/Simulation_logistic')


N=500
set.seed(10005)
#Repeat M times
M=100000

#Scenariob: covariates with random error
#a: w/o treatment effect
# beta0 = 0
beta1 = 0
# beta2 = 1
#Save the results
m1_est = rep(NA, M); m2_est = rep(NA, M)
m1_std = rep(NA, M); m2_std = rep(NA, M)
m1_z = rep(NA, M); m2_z = rep(NA, M)
m1_reject = rep(NA, M); m2_reject = rep(NA, M)
p0_var = rep(NA, M); p0_err_var = rep(NA,M)
p0_mean = rep(NA, M); p0_err_mean = rep(NA,M)

for (i in 1:M) {
  #Generate covariate
  c = 1+rnorm(N,0,1.5)
  #Treatment arm
  arm = c(rep(0,N/2),rep(1,N/2))
  #Generate probabilities
  logit = beta1*arm+c        #here beta0=0 and beta2=1, so c is the logit in control arm
  p = exp(logit)/(1+exp(logit))
  logit0 = c
  p0 = exp(logit0)/(1+exp(logit0))
  c_err = c+rnorm(N,0,1)            #covariate with random error
  logit0_err = c_err
  p0_err = exp(logit0_err)/(1+exp(logit0_err))
  #Generate responses
  y = rbinom(N,1,p)
  mydata = data.frame(arm,c,c_err,y)
  
  #Fit unadjusted model
  m1 = glm(y ~ as.factor(arm), data=mydata, family='binomial')
  m1_est[i] = summary(m1)$coefficients[2,1]
  m1_std[i] = summary(m1)$coefficients[2,2]
  m1_z[i] = summary(m1)$coefficients[2,3]
  m1_reject[i] = summary(m1)$coefficients[2,4]<0.05
  #Fit adjusted model
  m2 = glm(y ~ as.factor(arm)+c_err, data=mydata, family='binomial')    #using the covariate/logit with error
  m2_est[i] = summary(m2)$coefficients[2,1]
  m2_std[i] = summary(m2)$coefficients[2,2]
  m2_z[i] = summary(m2)$coefficients[2,3]
  m2_reject[i] = summary(m2)$coefficients[2,4]<0.05
  p0_var[i] = var(p0)
  p0_mean[i] = mean(p0)
  p0_err_var[i] = var(p0_err)
  p0_err_mean[i] = mean(p0_err)
}

#Type I error
mean(m1_z>qnorm(0.975) | m1_z<qnorm(0.025))
t1err_m1 = sum(m1_reject)/M
mean(m2_z>qnorm(0.975) | m2_z<qnorm(0.025))
t1err_m2 = sum(m2_reject)/M



########################
#b: w/ treatment effect
# beta0 = 0
beta1 = 0.75
# beta2 = 1
#Save the results
m1_est = rep(NA, M); m2_est = rep(NA, M)
m1_std = rep(NA, M); m2_std = rep(NA, M)
m1_z = rep(NA, M); m2_z = rep(NA, M)
m1_reject = rep(NA, M); m2_reject = rep(NA, M)
p0_var = rep(NA, M); p0_err_var = rep(NA,M)
p0_mean = rep(NA, M); p0_err_mean = rep(NA,M)
cov_cor = rep(NA,M)

for (i in 1:M) {
  #Generate covariate
  c = 1+rnorm(N,0,1.5)
  #Treatment arm
  arm = c(rep(0,N/2),rep(1,N/2))
  #Generate probabilities
  logit = beta1*arm+c        #here beta0=0 and beta2=1, so c is the logit in control arm
  p = exp(logit)/(1+exp(logit))
  logit0 = c
  p0 = exp(logit0)/(1+exp(logit0))
  c_err = c+rnorm(N,0,1)            #covariate with random error
  logit0_err = c_err
  p0_err = exp(logit0_err)/(1+exp(logit0_err))
  cov_cor[i] = cor(p0, p0_err)
  #Generate responses
  y = rbinom(N,1,p)
  mydata = data.frame(arm,c,c_err,y)
  
  #Fit unadjusted model
  m1 = glm(y ~ as.factor(arm), data=mydata, family='binomial')
  m1_est[i] = summary(m1)$coefficients[2,1]
  m1_std[i] = summary(m1)$coefficients[2,2]
  m1_z[i] = summary(m1)$coefficients[2,3]
  m1_reject[i] = summary(m1)$coefficients[2,4]<0.05
  #Fit adjusted model
  m2 = glm(y ~ as.factor(arm)+c_err, data=mydata, family='binomial')    #using the covariate/logit with error
  m2_est[i] = summary(m2)$coefficients[2,1]
  m2_std[i] = summary(m2)$coefficients[2,2]
  m2_z[i] = summary(m2)$coefficients[2,3]
  m2_reject[i] = summary(m2)$coefficients[2,4]<0.05
  p0_var[i] = var(p0)
  p0_mean[i] = mean(p0)
  p0_err_var[i] = var(p0_err)
  p0_err_mean[i] = mean(p0_err)
}

#Variance and 'bias' are inflated by the same factor
bias_factor = mean(m1_est/m2_est)
bias_factor_std = sqrt(var(m1_est/m2_est))
var_factor = mean(m1_std^2/m2_std^2)
var_factor_std = sqrt(var(m1_std^2/m2_std^2))
#The theoretical efficiency factor
eff_factor_err = mean(sqrt(1-p0_err_var/(p0_err_mean*(1-p0_err_mean))))
eff_factor_err_std = sqrt(var(sqrt(1-p0_err_var/(p0_err_mean*(1-p0_err_mean)))))
eff_factor = mean(sqrt(1-p0_var/(p0_mean*(1-p0_mean))))
eff_factor_std = sqrt(var(sqrt(1-p0_var/(p0_mean*(1-p0_mean)))))
#Adjust the expectation by correlation between p0 and p0_err
tmp = p0_var/(p0_mean*(1-p0_mean))
eff_factor_adj = mean(sqrt(1-tmp*cov_cor^2))    #expected efficiency factor adjusted by correlation
eff_factor_adj_std = sqrt(var(sqrt(1-tmp*cov_cor^2)))

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
      round(mean(p0_mean),2), round(mean(p0_var),2), round(eff_factor_adj,2), round(aeff_factor,2))
cbind(bias_factor_std, var_factor_std, sqrt(var(p0_mean)), eff_factor_adj_std, aeff_factor_std)

save(list=c('m1_est','m2_est','m1_std','m2_std','m1_z','m2_z',
            'm1_reject','m2_reject','t1err_m1','t1err_m2',
            'p0_var','p0_err_var','p0_mean','p0_err_mean','cov_cor'), file='logistic_randerr.rdata')


