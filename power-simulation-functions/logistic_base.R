
# Study the effect of covariate adjustment on type I error rate and power
# this is the base case
# model 2 is correctly specified

rm(list=ls())
setwd('/Users/liyunfan/Documents/Simulation_logistic')

N=500
set.seed(9999)
#Repeat M times
M=100000

#Scenario: Efficiency gain is not related to treatment effect
#a: w/o treatment effect
beta0 = 1
beta1 = 0
beta2 = 1
#Save the results
m1_est = rep(NA, M); m2_est = rep(NA, M)
m1_std = rep(NA, M); m2_std = rep(NA, M)
m1_z = rep(NA, M); m2_z = rep(NA, M)
m1_reject = rep(NA, M); m2_reject = rep(NA, M)
p0_var = rep(NA, M)
p0_mean = rep(NA, M)

for (i in 1:M) {
    #Generate covariate
    c = rnorm(N,0,1.5)
    #Treatment arm
    arm = c(rep(0,N/2),rep(1,N/2))
    #Generate probabilities
    logit = beta0+beta1*arm+beta2*c
    p = exp(logit)/(1+exp(logit))
    logit0 = beta0+beta2*c
    p0 = exp(logit0)/(1+exp(logit0))    #p0 are the true logit, if subj in control arm
    #Generate responses
    y = rbinom(N,1,p)
    mydata = data.frame(arm,c,y)

    #Fit unadjusted model
    m1 = glm(y ~ as.factor(arm), data=mydata, family='binomial')
    m1_est[i] = summary(m1)$coefficients[2,1]
    m1_std[i] = summary(m1)$coefficients[2,2]
    m1_z[i] = summary(m1)$coefficients[2,3]
    m1_reject[i] = summary(m1)$coefficients[2,4]<0.05
    #Fit adjusted model, using the true logit
    m2 = glm(y ~ as.factor(arm)+c, data=mydata, family='binomial')
    m2_est[i] = summary(m2)$coefficients[2,1]
    m2_std[i] = summary(m2)$coefficients[2,2]
    m2_z[i] = summary(m2)$coefficients[2,3]
    m2_reject[i] = summary(m2)$coefficients[2,4]<0.05
    p0_var[i] = var(p0)
    p0_mean[i] = mean(p0)
    }

#Type I error
mean(m1_z>qnorm(0.975) | m1_z<qnorm(0.025))
t1err_m1 = sum(m1_reject)/M
mean(m2_z>qnorm(0.975) | m2_z<qnorm(0.025))
t1err_m2 = sum(m2_reject)/M

########################
#b: w/ treatment effect
beta0 = 1
beta1 = 0.75
beta2 = 1
#Save the results
m1_est = rep(NA, M); m2_est = rep(NA, M)
m1_std = rep(NA, M); m2_std = rep(NA, M)
m1_z = rep(NA, M); m2_z = rep(NA, M)
m1_reject = rep(NA, M); m2_reject = rep(NA, M)
p0_var = rep(NA, M)
p0_mean = rep(NA, M)

for (i in 1:M) {
    #Generate covariate
    c = rnorm(N,0,1.5)
    #Treatment arm
    arm = c(rep(0,N/2),rep(1,N/2))
    #Generate probabilities
    logit = beta0+beta1*arm+beta2*c
    p = exp(logit)/(1+exp(logit))
    logit0 = beta0+beta2*c
    p0 = exp(logit0)/(1+exp(logit0))
    #Generate responses
    y = rbinom(N,1,p)
    mydata = data.frame(arm,c,y)

    #Fit unadjusted model
    m1 = glm(y ~ as.factor(arm), data=mydata, family='binomial')
    m1_est[i] = summary(m1)$coefficients[2,1]
    m1_std[i] = summary(m1)$coefficients[2,2]
    m1_z[i] = summary(m1)$coefficients[2,3]
    m1_reject[i] = summary(m1)$coefficients[2,4]<0.05
    #Fit adjusted model
    m2 = glm(y ~ as.factor(arm)+c, data=mydata, family='binomial')
    m2_est[i] = summary(m2)$coefficients[2,1]
    m2_std[i] = summary(m2)$coefficients[2,2]
    m2_z[i] = summary(m2)$coefficients[2,3]
    m2_reject[i] = summary(m2)$coefficients[2,4]<0.05
    p0_var[i] = var(p0)
    p0_mean[i] = mean(p0)
    }

#Variance and 'bias' are inflated by the same factor
bias_factor1 = mean(m1_est/m2_est)
bias_factor1_std = sqrt(var(m1_est/m2_est))
var_factor1 = mean(m1_std^2/m2_std^2)
var_factor1_std = sqrt(var(m1_std^2/m2_std^2))
#The theoretical efficiency factor
eff_factor1 = mean(sqrt(1-p0_var/(p0_mean*(1-p0_mean))))
eff_factor1_std = sqrt(var(sqrt(1-p0_var/(p0_mean*(1-p0_mean)))))
p0b_mean = mean(p0_mean)
p0b_mean_std = sqrt(var(p0_mean))
p0b_var = mean(p0_var)
#The actual efficiency factor (ratio of z scores)
aeff_factor1 = mean(m1_z/m2_z)
aeff_factor1_std = sqrt(var(m1_z/m2_z))

#Power
mean(m1_z>qnorm(0.975) | m1_z<qnorm(0.025))
power1_m1 = sum(m1_reject)/M
mean(m2_z>qnorm(0.975) | m2_z<qnorm(0.025))
power1_m2 = sum(m2_reject)/M


#################################
#c: w/ a larger treatment effect
beta0 = 1
beta1 = 0.85
beta2 = 1
#Save the results
m1_est = rep(NA, M); m2_est = rep(NA, M)
m1_std = rep(NA, M); m2_std = rep(NA, M)
m1_z = rep(NA, M); m2_z = rep(NA, M)
m1_reject = rep(NA, M); m2_reject = rep(NA, M)
p0_var = rep(NA, M)
p0_mean = rep(NA, M)

for (i in 1:M) {
  #Generate covariate
  c = rnorm(N,0,1.5)
  #Treatment arm
  arm = c(rep(0,N/2),rep(1,N/2))
  #Generate probabilities
  logit = beta0+beta1*arm+beta2*c
  p = exp(logit)/(1+exp(logit))
  logit0 = beta0+beta2*c
  p0 = exp(logit0)/(1+exp(logit0))
  #Generate responses
  y = rbinom(N,1,p)
  mydata = data.frame(arm,c,y)
  
  #Fit unadjusted model
  m1 = glm(y ~ as.factor(arm), data=mydata, family='binomial')
  m1_est[i] = summary(m1)$coefficients[2,1]
  m1_std[i] = summary(m1)$coefficients[2,2]
  m1_z[i] = summary(m1)$coefficients[2,3]
  m1_reject[i] = summary(m1)$coefficients[2,4]<0.05
  #Fit adjusted model
  m2 = glm(y ~ as.factor(arm)+c, data=mydata, family='binomial')
  m2_est[i] = summary(m2)$coefficients[2,1]
  m2_std[i] = summary(m2)$coefficients[2,2]
  m2_z[i] = summary(m2)$coefficients[2,3]
  m2_reject[i] = summary(m2)$coefficients[2,4]<0.05
  p0_var[i] = var(p0)
  p0_mean[i] = mean(p0)
}

#Variance and 'bias' are inflated by the same factor
bias_factor2 = mean(m1_est/m2_est)
bias_factor2_std = sqrt(var(m1_est/m2_est))
var_factor2 = mean(m1_std^2/m2_std^2)
var_factor2_std = sqrt(var(m1_std^2/m2_std^2))
#The theoretical efficiency factor
eff_factor2 = mean(sqrt(1-p0_var/(p0_mean*(1-p0_mean))))
eff_factor2_std = sqrt(var(sqrt(1-p0_var/(p0_mean*(1-p0_mean)))))
p0c_mean = mean(p0_mean)
p0c_mean_std = sqrt(var(p0_mean))
p0c_var = mean(p0_var)
#The actual efficiency factor (ratio of z scores)
aeff_factor2 = mean(m1_z/m2_z)
aeff_factor2_std = sqrt(var(m1_z/m2_z))


#Power
mean(m1_z>qnorm(0.975) | m1_z<qnorm(0.025))
power2_m1 = sum(m1_reject)/M
mean(m2_z>qnorm(0.975) | m2_z<qnorm(0.025))
power2_m2 = sum(m2_reject)/M


cbind(round(t1err_m1,3)*100, round(t1err_m2,3)*100)
cbind(round(power1_m1,3)*100, round(power1_m2,3)*100,
      round(bias_factor1,2), round(var_factor1,2),
      round(p0b_mean,2), round(p0b_var,2), round(eff_factor1,2), round(aeff_factor1,2))
cbind(bias_factor1_std, var_factor1_std, p0b_mean_std, eff_factor1_std, aeff_factor1_std)
cbind(round(power2_m1,3)*100, round(power2_m2,3)*100,
      round(bias_factor2,2), round(var_factor2,2),
      round(p0c_mean,2), round(p0c_var,2), round(eff_factor2,2), round(aeff_factor2,2))
cbind(bias_factor2_std, var_factor2_std, p0c_mean_std, eff_factor2_std, aeff_factor2_std)


