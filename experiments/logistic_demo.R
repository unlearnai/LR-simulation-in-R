
# A demo of sample size reduction in logistic regression
# Balanced design, equal reduction in each arm

rm(list=ls())
setwd('/Users/liyunfan/Documents/Simulation_logistic')

load("logistic_demo.RData")

N=2000
#Scenario 1: With treatment effect; equal size arms
beta0 = 0.5
beta1 = 0.5
beta2 = 2
#Repeat M times
set.seed(9999)
M=10000
m1_est = rep(NA,ncol=M); m2_est = rep(NA,ncol=M)
m1_std = rep(NA,ncol=M); m2_std = rep(NA,ncol=M)
m1_z = rep(NA,ncol=M); m2_z = rep(NA,ncol=M)
# m1_coverage = rep(NA,M); m2_coverage = rep(NA,M)
m1_reject = rep(NA,ncol=M); m2_reject = rep(NA,ncol=M)
p0_var = rep(NA,ncol=M)
p0_mean = rep(NA,ncol=M)
inf_factor = rep(NA,ncol=M)
m1_trteff = rep(NA,ncol=M); m2_trteff = rep(NA,ncol=M)
for (i in 1:M) {
    #Generate covariate
    c = rnorm(N,0,1)
    #Treatment arm
    arm = c(rep(0,N/2),rep(1,N/2))
    #Generate probabilities
    logit = beta0+beta1*arm+beta2*c
    p = exp(logit)/(1+exp(logit))
    logit0 = beta0+beta2*c
    p0 = exp(logit0)/(1+exp(logit0))
    #Generate responses
    y = rbinom(N,1,p)
    mydata = data.frame(arm,c,y, p0)
    #Inflation factor by covariate adjustment
    inf_factor[i] = 1-var(p0)/(mean(p0)*(1-mean(p0)))
    p0_var[i] = var(p0)
    p0_mean[i] = mean(p0)
    
    #Fit unadjusted model
    m1 = glm(y ~ as.factor(arm), data=mydata, family='binomial')
    m1_est[i] = summary(m1)$coefficients[2,1]
    m1_std[i] = summary(m1)$coefficients[2,2]
    m1_z[i] = summary(m1)$coefficients[2,3]
    # m1_coverage[i] = confint(m1)[2,1]<beta1 & confint(m1)[2,2]>beta1
    m1_reject[i] = summary(m1)$coefficients[2,4]<0.05
    #Calculate treatment effect in unadjusted model
    mydata0 = cbind(rep(0,N),mydata[,-1])
    colnames(mydata0)[1] = 'arm'
    mydata1 = cbind(rep(1,N),mydata[,-1])
    colnames(mydata1)[1] = 'arm'
    m1_trteff[i] = mean(predict(m1, mydata1, type='response'))-mean(predict(m1, mydata0, type='response'))
    
    #Take subset of data and fit adjusted model
    sub_ind = sample(1:N, ceiling(N*inf_factor[i]), replace=F)
    sub_data = mydata[sub_ind,]
    m2 = glm(y ~ as.factor(arm)+p0, data=sub_data, family='binomial')
    m2_est[i] = summary(m2)$coefficients[2,1]
    m2_std[i] = summary(m2)$coefficients[2,2]
    m2_z[i] = summary(m2)$coefficients[2,3]
    # m2_coverage[i] = confint(m2)[2,1]<beta1 & confint(m2)[2,2]>beta1
    m2_reject[i] = summary(m2)$coefficients[2,4]<0.05
    #Calculate marginal treatment effect in adjusted model
    sub_data0 = cbind(rep(0,length(sub_ind)),sub_data[,-1])
    colnames(sub_data0)[1] = 'arm'
    sub_data1 = cbind(rep(1,length(sub_ind)),sub_data[,-1])
    colnames(sub_data1)[1] = 'arm'
    m2_trteff[i] = mean(predict(m2, sub_data1, type='response'))-mean(predict(m2, sub_data0, type='response'))
}

summary(m1_z)
summary(m2_z)
sum(m1_reject)
sum(m2_reject)
mean(inf_factor)
summary(m1_trteff)
summary(m2_trteff)

# save.image(file = "logistic_demo.RData")
