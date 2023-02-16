Functions for PROCOVA-LR Simulation Studies on Power

Six different settings for the simulation studies are considered. The settings are ‘Base’, ‘Large Variance’, ‘Higher Prevalence’, ‘Random Error’, ‘Shifted’, and ‘Omitted Covariate’. Each is implemented in R via a corresponding .R file.

In the first set of simulation settings (‘Base’, ‘Large Variance’, and ‘Higher Prevalence’), the covariate-adjusted model is correctly specified. The data generating model contains a covariate, the covariate-adjusted analysis model is a logistic regression model specified exactly as the data generating model, and the unadjusted model omits the covariate from the analysis. The scenarios are summarized below, where each summary contains the true data generating model, X denotes treatment assignment, and U denotes a known covariate.

In the first setting, we consider three situations, each simulated data set consists of 500 subjects, with 1:1 randomization for the two arms. The three situations for this setting are:

 \mathrm{Scenario \ 1a:} \ \mathrm{Pr}(Y_i^{\mathrm{obs}} = 1 \mid X_i, U_i) = \frac{\mathrm{exp} \left ( 1 + U_i \right )}{1 + \mathrm{exp} \left ( 1 + U_i \right )} \; \mathrm{where} \; U_i \sim \mathrm{Normal}(0,1.5^2) \ \mathrm{independently},

\mathrm{Scenario \ 1b:} \ \mathrm{Pr}(Y_i^{\mathrm{obs}} = 1 \mid X_i, U_i) = \frac{\mathrm{exp} \left ( 1 + 0.75X + U_i \right )}{1 + \mathrm{exp} \left ( 1 + 0.75X_i + U_i \right )} \; \mathrm{where} \; U_i \sim \mathrm{Normal}(0,1.5^2) \ \mathrm{independently}, 

\mathrm{Scenario \ 1c:} \ \mathrm{Pr}(Y_i = 1 \mid X_i, U_i) = \frac{\mathrm{exp} \left ( 1 + 0.85X_i + U_i \right )}{1 + \mathrm{exp} \left ( 1 + 0.85X_i + U_i \right )} \; \mathrm{where} \; U_i \sim \mathrm{Normal}(0,1.5^2)\ \mathrm{independently}.
