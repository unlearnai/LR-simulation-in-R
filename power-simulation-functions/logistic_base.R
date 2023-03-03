
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
m1_reject = rep(NA, M); m2_reject = rep(NA, M)
m1_RD_reject = rep(NA, M); m2_RD_reject = rep(NA, M)
m1_lgRR_reject = rep(NA, M); m2_lgRR_reject = rep(NA, M)

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
    m1_reject[i] = summary(m1)$coefficients[2,4]<0.05
    #Fit adjusted model, using the true logit
    m2 = glm(y ~ as.factor(arm)+c, data=mydata, family='binomial')
    m2_reject[i] = summary(m2)$coefficients[2,4]<0.05

    #Calculate the Jacobian and Wald statistic for risk difference (RD) -adjusted model
    tmp1 = exp(cbind(rep(1,N),rep(1,N),c)%*%m2$coefficients)
    tmp0 = exp(cbind(rep(1,N),rep(0,N),c)%*%m2$coefficients)
    denom1_m2 = (1+tmp1)^2    #this denominator will be used a lot
    denom0_m2 = (1+tmp0)^2    #beta1 doesn't appear in this denominator
    Jm2_RD_beta0 = mean(tmp1/denom1_m2)-mean(tmp0/denom0_m2)
    Jm2_RD_beta1 = mean(tmp1/denom1_m2)
    Jm2_RD_beta2 = mean(c*tmp1/denom1_m2)-mean(c*tmp0/denom0_m2)
    Jm2_RD = c(Jm2_RD_beta0, Jm2_RD_beta1, Jm2_RD_beta2)    #this is the Jacobian for RD
    m2_RD_std = sqrt(Jm2_RD%*%vcov(m2)%*%c(Jm2_RD))    #this is the std of RD estimate
    m2_RD_est = mean(tmp1/(1+tmp1))-mean(tmp0/(1+tmp0))    #this is the estimate of RD
    m2_RD_z = m2_RD_est/m2_RD_std
    m2_RD_reject[i] = m2_RD_z>qnorm(0.975) | m2_RD_z<qnorm(0.025)
    
    #Calculate the Jacobian and Wald statistic for log relative risk (RR) -adjusted model
    Jm2_RR_beta0 = mean(tmp1/denom1_m2)/mean(tmp1/(1+tmp1))-mean(tmp0/denom0_m2)/(mean(tmp0/(1+tmp0)))
    Jm2_RR_beta1 = mean(tmp1/denom1_m2)/mean(tmp1/(1+tmp1))
    Jm2_RR_beta2 = mean(c*tmp1/denom1_m2)/mean(tmp1/(1+tmp1))-mean(c*tmp0/denom0_m2)/mean(tmp0/(1+tmp0))
    Jm2_RR = c(Jm2_RR_beta0, Jm2_RR_beta1, Jm2_RR_beta2)    #this is the Jacobian for RR
    m2_lgRR_std = sqrt(Jm2_RR%*%vcov(m2)%*%c(Jm2_RR))    #this is the std of log RR estimate
    m2_lgRR_est = log(mean(tmp1/(1+tmp1)))-log(mean(tmp0/(1+tmp0)))    #this is the estimate of RR
    m2_lgRR_z = m2_lgRR_est/m2_lgRR_std
    m2_lgRR_reject[i] = m2_lgRR_z>qnorm(0.975) | m2_lgRR_z<qnorm(0.025)
    
    #Calculate the Jacobian and Wald statistic for risk difference (RD) -unadjusted model
    tmp1 = exp(sum(m1$coefficients))
    tmp0 = exp(m1$coefficients[1])
    denom1_m1 = (1+tmp1)^2    #this denominator will be used a lot
    denom0_m1 = (1+tmp0)^2    #beta1 doesn't appear in this denominator
    Jm1_RD_beta0 = tmp1/denom1_m1-tmp0/denom0_m1
    Jm1_RD_beta1 = tmp1/denom1_m1
    Jm1_RD = c(Jm1_RD_beta0, Jm1_RD_beta1)    #this is the Jacobian for RD
    m1_RD_std = sqrt(Jm1_RD%*%vcov(m1)%*%c(Jm1_RD))    #this is the std of RD estimate
    m1_RD_est = mean(tmp1/(1+tmp1))-mean(tmp0/(1+tmp0))    #this is the estimate of RD
    m1_RD_z = m1_RD_est/m1_RD_std
    m1_RD_reject[i] = m1_RD_z>qnorm(0.975) | m1_RD_z<qnorm(0.025)
    
    #Calculate the Jacobian and Wald statistic for log risk difference (RR) -unadjusted model
    Jm1_RR_beta0 = (tmp1/denom1_m1)/(tmp1/(1+tmp1))-(tmp0/denom0_m1)/(tmp0/(1+tmp0))
    Jm1_RR_beta1 = (tmp1/denom1_m1)/(tmp1/(1+tmp1))
    Jm1_RR = c(Jm1_RR_beta0, Jm1_RR_beta1)    #this is the Jacobian for log RR
    m1_lgRR_std = sqrt(Jm1_RR%*%vcov(m1)%*%c(Jm1_RR))    #this is the std of log RR estimate
    m1_lgRR_est = log(mean(tmp1/(1+tmp1)))-log(mean(tmp0/(1+tmp0)))    #this is the estimate of RR
    m1_lgRR_z = m1_lgRR_est/m1_lgRR_std
    m1_lgRR_reject[i] = m1_lgRR_z>qnorm(0.975) | m1_lgRR_z<qnorm(0.025)
}

#Report the type I error
round(cbind(mean(m1_reject), mean(m2_reject), mean(m1_RD_reject), mean(m2_RD_reject), 
      mean(m1_lgRR_reject), mean(m2_lgRR_reject))*100,2)

########################
#b: w/ treatment effect
beta0 = 1
beta1 = 0.75
beta2 = 1
#Save the results
p0_var = rep(NA, M)
p0_mean = rep(NA, M)
m1_est = rep(NA, M); m2_est = rep(NA, M)
m1_std = rep(NA, M); m2_std = rep(NA, M)
m1_z = rep(NA, M); m2_z = rep(NA, M)
m1_reject_tmp = rep(NA, M); m2_reject_tmp = rep(NA, M)
m1_RD_est = rep(NA, M); m2_RD_est = rep(NA, M)
m1_RD_std = rep(NA, M); m2_RD_std = rep(NA, M)
m1_RD_z = rep(NA, M); m2_RD_z = rep(NA, M)
m1_lgRR_est = rep(NA, M); m2_lgRR_est = rep(NA, M)
m1_lgRR_std = rep(NA, M); m2_lgRR_std = rep(NA, M)
m1_lgRR_z = rep(NA, M); m2_lgRR_z = rep(NA, M)

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
    #Mean and variance of p0
    p0_var[i] = var(p0)
    p0_mean[i] = mean(p0)
    
    #Fit unadjusted model
    m1 = glm(y ~ as.factor(arm), data=mydata, family='binomial')
    m1_est[i] = summary(m1)$coefficients[2,1]
    m1_std[i] = summary(m1)$coefficients[2,2]
    m1_z[i] = summary(m1)$coefficients[2,3]
    m1_reject_tmp[i] = summary(m1)$coefficients[2,4]<0.05
    #Fit adjusted model
    m2 = glm(y ~ as.factor(arm)+c, data=mydata, family='binomial')
    m2_est[i] = summary(m2)$coefficients[2,1]
    m2_std[i] = summary(m2)$coefficients[2,2]
    m2_z[i] = summary(m2)$coefficients[2,3]
    m2_reject_tmp[i] = summary(m2)$coefficients[2,4]<0.05
    
    #Calculate the Jacobian and Wald statistic for risk difference (RD) -adjusted model
    tmp1 = exp(cbind(rep(1,N),rep(1,N),c)%*%m2$coefficients)
    tmp0 = exp(cbind(rep(1,N),rep(0,N),c)%*%m2$coefficients)
    denom1_m2 = (1+tmp1)^2    #this denominator will be used a lot
    denom0_m2 = (1+tmp0)^2    #beta1 doesn't appear in this denominator
    Jm2_RD_beta0 = mean(tmp1/denom1_m2)-mean(tmp0/denom0_m2)
    Jm2_RD_beta1 = mean(tmp1/denom1_m2)
    Jm2_RD_beta2 = mean(c*tmp1/denom1_m2)-mean(c*tmp0/denom0_m2)
    Jm2_RD = c(Jm2_RD_beta0, Jm2_RD_beta1, Jm2_RD_beta2)    #this is the Jacobian for RD
    m2_RD_std[i] = sqrt(Jm2_RD%*%vcov(m2)%*%c(Jm2_RD))    #this is the std of RD estimate
    m2_RD_est[i] = mean(tmp1/(1+tmp1))-mean(tmp0/(1+tmp0))    #this is the estimate of RD
    m2_RD_z[i] = m2_RD_est[i]/m2_RD_std[i]
    
    #Calculate the Jacobian and Wald statistic for log relative risk (RR) -adjusted model
    Jm2_RR_beta0 = mean(tmp1/denom1_m2)/mean(tmp1/(1+tmp1))-mean(tmp0/denom0_m2)/(mean(tmp0/(1+tmp0)))
    Jm2_RR_beta1 = mean(tmp1/denom1_m2)/mean(tmp1/(1+tmp1))
    Jm2_RR_beta2 = mean(c*tmp1/denom1_m2)/mean(tmp1/(1+tmp1))-mean(c*tmp0/denom0_m2)/mean(tmp0/(1+tmp0))
    Jm2_RR = c(Jm2_RR_beta0, Jm2_RR_beta1, Jm2_RR_beta2)    #this is the Jacobian for RR
    m2_lgRR_std[i] = sqrt(Jm2_RR%*%vcov(m2)%*%c(Jm2_RR))    #this is the std of log RR estimate
    m2_lgRR_est[i] = log(mean(tmp1/(1+tmp1)))-log(mean(tmp0/(1+tmp0)))    #this is the estimate of RR
    m2_lgRR_z[i] = m2_lgRR_est[i]/m2_lgRR_std[i]
    
    #Calculate the Jacobian and Wald statistic for risk difference (RD) -unadjusted model
    tmp1 = exp(sum(m1$coefficients))
    tmp0 = exp(m1$coefficients[1])
    denom1_m1 = (1+tmp1)^2    #this denominator will be used a lot
    denom0_m1 = (1+tmp0)^2    #beta1 doesn't appear in this denominator
    Jm1_RD_beta0 = tmp1/denom1_m1-tmp0/denom0_m1
    Jm1_RD_beta1 = tmp1/denom1_m1
    Jm1_RD = c(Jm1_RD_beta0, Jm1_RD_beta1)    #this is the Jacobian for RD
    m1_RD_std[i] = sqrt(Jm1_RD%*%vcov(m1)%*%c(Jm1_RD))    #this is the std of RD estimate
    m1_RD_est[i] = mean(tmp1/(1+tmp1))-mean(tmp0/(1+tmp0))    #this is the estimate of RD
    m1_RD_z[i] = m1_RD_est[i]/m1_RD_std[i]
    
    #Calculate the Jacobian and Wald statistic for log risk difference (RR) -unadjusted model
    Jm1_RR_beta0 = (tmp1/denom1_m1)/(tmp1/(1+tmp1))-(tmp0/denom0_m1)/(tmp0/(1+tmp0))
    Jm1_RR_beta1 = (tmp1/denom1_m1)/(tmp1/(1+tmp1))
    Jm1_RR = c(Jm1_RR_beta0, Jm1_RR_beta1)    #this is the Jacobian for log RR
    m1_lgRR_std[i] = sqrt(Jm1_RR%*%vcov(m1)%*%c(Jm1_RR))    #this is the std of log RR estimate
    m1_lgRR_est[i] = log(mean(tmp1/(1+tmp1)))-log(mean(tmp0/(1+tmp0)))    #this is the estimate of RR
    m1_lgRR_z[i] = m1_lgRR_est[i]/m1_lgRR_std[i]
}

#Variance and 'bias' are inflated by the same factor
bias_factor = mean(m1_est/m2_est)
bias_factor_std = sqrt(var(m1_est/m2_est))
var_factor = mean(m1_std^2/m2_std^2)
var_factor_std = sqrt(var(m1_std^2/m2_std^2))
#The theoretical efficiency factor
eff_factor = mean(sqrt(1-p0_var/(p0_mean*(1-p0_mean))))
eff_factor_std = sqrt(var(sqrt(1-p0_var/(p0_mean*(1-p0_mean)))))
E_p0 = mean(p0_mean)
E_p0_std = sqrt(var(p0_mean))
Var_p0 = mean(p0_var)

#############Properties of the model###############
#Report the bias factor, ARE, mean and variance of p0
round(cbind(bias_factor, var_factor, eff_factor, E_p0, Var_p0),2)
#Report the std of the above
round(cbind(bias_factor_std, var_factor_std, eff_factor_std, E_p0_std),2)



#############Inference on beta_1^hat###############
#Power of test on \beta_1^hat
power_m1 = sum(m1_reject_tmp)/M
power_m2 = sum(m2_reject_tmp)/M
#The actual efficiency factor (ratio of z scores)
z_ratio = mean(m1_z/m2_z)
z_ratio_std = sqrt(var(m1_z/m2_z))

#############Inference on RD###############
#Power of test on RD
power_RD_m1 = sum(m1_RD_z>qnorm(0.975) | m1_RD_z<qnorm(0.025))/M
power_RD_m2 = sum(m2_RD_z>qnorm(0.975) | m2_RD_z<qnorm(0.025))/M
#The actual efficiency factor (ratio of z scores)
z_RD_ratio = mean(m1_RD_z/m2_RD_z)
z_RD_ratio_std = sqrt(var(m1_RD_z/m2_RD_z))

#############Inference on log RR###############
#Power of test on log RR
power_lgRR_m1 = sum(m1_lgRR_z>qnorm(0.975) | m1_lgRR_z<qnorm(0.025))/M
power_lgRR_m2 = sum(m2_lgRR_z>qnorm(0.975) | m2_lgRR_z<qnorm(0.025))/M
#The actual efficiency factor (ratio of z scores)
z_lgRR_ratio = mean(m1_lgRR_z/m2_lgRR_z)
z_lgRR_ratio_std = sqrt(var(m1_lgRR_z/m2_lgRR_z))

#Report the mean and variance of p0, expected efficiency factor, and ratio of z-scores and power
#in all three tests
round(cbind(E_p0, Var_p0, eff_factor, z_ratio, power_m1*100, power_m2*100, z_RD_ratio, power_RD_m1*100, power_RD_m2*100,
            z_lgRR_ratio, power_lgRR_m1*100, power_lgRR_m2*100),2)
#Report the std of the above
round(cbind(E_p0_std, eff_factor_std, z_ratio_std, z_RD_ratio_std, z_lgRR_ratio_std),2)


save(list=c('m1_est','m2_est','m1_std','m2_std','m1_z','m2_z',
            'm1_reject','m2_reject','m1_RD_reject','m2_RD_reject','m1_lgRR_reject','m2_lgRR_reject',
            'p0_var','p0_mean',
            'm1_RD_est','m2_RD_est','m1_RD_std','m2_RD_std','m1_RD_z','m2_RD_z',
            'm1_lgRR_est','m2_lgRR_est','m1_lgRR_std','m2_lgRR_std',
            'm1_lgRR_z','m2_lgRR_z'), file='logistic_base.rdata')



