#Goal: To estimate the correlation between true and estimated underlying probabilities, to
#adjust for over confidence in the efficiency factor
#Experiment: To group participants into bins, calculate the observed probability within each bin
#from binary outcomes, and calculate the correlation between the observed bin probabilities and
#estimated bin probabiilities

rm(list=ls())
library(MASS)
library(boot)
library(dplyr)

N=500

u = 1+rnorm(N,0,1.5)
#Assume control arm data only; no treatment arm and no treatment effect
#Generate probabilities
logit0 = u        #u is the true underlying logit probability in controlled participants
p0 = exp(logit0)/(1+exp(logit0))
c_err = u+rnorm(N,0,1)            #logit with random error
logit0_err = c_err
p0_err = exp(logit0_err)/(1+exp(logit0_err))
#Generate responses
y = rbinom(N,1,p0)
mydata = data.frame(u,c_err,p0,p0_err,y)

#This is the probability we want to estimate, without knowing p0
cor(p0,p0_err)

num_bins = 50
num_rep = 10000
est_cor = rep(NA,num_rep)

for (i in 1:num_rep) {
  #Group the participants into bins
  bin_idx = sample(1:num_bins,N,replace=T)
  bin_mydata = data.frame(mydata,bin_idx)
  #Calculate the probabilities within each bin
  bin_probabilities = bin_mydata %>% group_by(bin_idx) %>%
    summarise_at(vars(c(y,p0,p0_err)), list(name=mean))
  # cor(bin_probabilities$p0_name,bin_probabilities$p0_err_name)
  est_cor[i] = cor(bin_probabilities$y_name,bin_probabilities$p0_err_name)
}
summary(est_cor)
