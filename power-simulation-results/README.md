**Results of PROCOVA-LR Simulation Studies on Power**

Results of the simulation studies for six different settings are provided. These settings are stored in the ‘power-simulation-functions’ folder.
For each setting, the simulation study results are stored in the “logistic_settingname.RData" file.
Each logistic_settingname.RData file contains the following variables:

m1_est: estimated beta coefficient () in the unadjusted model
m1_reject: whether the null hypothesis is rejected (=0) in the unadjusted model
m1_std: standard error of the beta coefficient estimate in the unadjusted model
m1_z: z-score for the hypothesis testing in the unadjusted model
m2_est: estimated beta coefficient () in the adjusted model
m2_reject: whether the null hypothesis is rejected (=0) in the adjusted model
m2_std: standard error of the beta coefficient estimate in the adjusted model
m2_z: z-score for the hypothesis testing in the unadjusted model
p0_mean: mean of the underlying true probabilities
p0_var: variance of the underlying true probabilities
t1err_m1: type I error in the unadjusted model
t1err_m2: type I error in the adjusted model
In the settings with misspecified models (‘Random Error’, ‘Shifted’, ‘Omitted Covariate’), the following variables are also included:
cov_cor: correlation between true and available underlying probabilities to adjust the expected efficiency factor
p0_err_mean: mean of the available underlying probabilities
p0_err_var: variance of the available underlying probabilities
