
# Generate a plot showing the relationship between power gain and the efficiency factor
# in Wald test
# This can be used in linear regression and logistic regression, by changing
# the definition of the efficiency factor

library(ggplot2)
library(tidyverse)

rm(list=ls())
setwd('/Users/liyunfan/Documents/Simulation_logistic')

# Function of Wald test power 
wald_pow = function(x) (pnorm(qnorm(0.025)+x) + pnorm(qnorm(0.025)-x))
seq_w = seq(0,5,by=0.01)
seq_wald_pow = lapply(X=seq_w, FUN=wald_pow)

#Now do a finer plot to zoom in the possible area
seq_f = seq(0.8,1,by=0.05)[5:1]
seq_adj_pow = matrix(, nrow=length(seq_f), ncol=length(seq_w))
for (i in 1:length(seq_f)) {
  f = seq_f[i]
  seq_adj_w = seq_w/f
  seq_adj_pow[i,] = unlist(lapply(X=seq_adj_w, FUN=wald_pow))
}

# Plot power gain by the efficiency factor
seq_adj_pow = t(seq_adj_pow)
colnames(seq_adj_pow) = sprintf("%.2f",seq_f)
all = data.frame(seq_w, seq_adj_pow)
colnames(all)[2:6] = sprintf("%.2f",seq_f)
pow.gathered = all %>%
  as_data_frame() %>%
  gather(key = 'efficiency', value='value',
         -seq_w)
power_curves = ggplot(pow.gathered, aes(x=seq_w, y=value)) +
  geom_line(aes(linetype=efficiency), size=1) +
  scale_linetype_manual(values=c('dotted','dashed','dotdash','twodash','solid')) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0.75, 1, by=0.025), limits=c(0.75,1)) +
  scale_x_continuous(breaks = seq(2.5, 4.5, by=0.25), limits=c(2.75,4.5)) +
  scale_color_viridis_d() +
  xlab('Unadjusted Wald statistic') + ylab('Power') +
  theme(text=element_text(size=15)) +
  labs(linetype=expression(f[EFF]))
ggsave('power_curves.eps', plot=power_curves, device='eps', dpi=800)


