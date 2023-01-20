rm(list=ls())
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

m = glm(y ~ as.factor(arm)+c, data=mydata, family='binomial')
tmp = predict(m, newdata=mydata, type='response')
y_b = rbinom(n=500,1,tmp)
data_b = data.frame(arm,c,y_b)
m_b = glm(y ~ as.factor(arm)+c, data=data_b, family='binomial')
m_b$coefficients
m$coefficients

