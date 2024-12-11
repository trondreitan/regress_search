# regress_search
Regression search algorithms in R, based on information criteria, a use-specified criteria or hypothesis testing. Either full dredge or step-wise.

Regression models allowed for:
A. Linear model (using 'lm') or mixed linear models (using 'lmm' in the 'lme4' package). Specify with family="Normal" 
B) General linear models (using 'glm'), mixed general linear models (so with random effects, using 'glmm' in the 'lme4' package or alternatively 'glmmadmb' in the 'glmmADMB' package or 'glmmTMB' in the 'glmmTMB' package). Specify with the option family="binomial", "poisson" or any other allowed in the 'glm' method.
c) Survival analysis (using 'clogit' in the 'survival' package). Specify with the option family="clogit".
d) Ordinal regression (using 'clm' or 'clmm' in the in the 'ordinal' package). Specify with the option family="ordinal" or family="ordinal_cauchyit".

The search algorithms per default only returns the top model, but more can be returned if specified.

Other standard options: 
a) data - a data frame to be sent to the regression mdethods.
b) respnse - The name of the response variable (should exist in the 'data' data frame).
c) covs - An array of covariate names (should exist in the 'data' data frame).
d) family - Described above.
e) IC - Information criterion to use (sefault "BIC"). Options are "BIC", "AIC" or "AICc".
f) top (default 1) - If larger than 1, returns the top models rather than just one, as a list object. 

More specialized options:
g) offset (if applicable)
h) use.glmmADMB/use.glmmTMB - For using alternative packages for mixed general linear models.i
i) talkative - Lots of output information about what models are being focused on, while running.
j) do.return.status. Returns the model specification (1 for covariates that are used, 0 for uniused) in addition to the best model (as a list).
k) strata - Strate treatment in the ordinal framwork.
l) threshold - Threshold option if using ordinal regression (using 'clm').
m) check.est - Checks if NA values are found in the coefficient estimates, in which case IC=NA.
n) check.se - Checks if NA values are found in the coefficient standard errors, in which case IC=NA.o
o) max.complex - Maximal complexity (if set).

Functions:
1) regress.ic.dredge: Goes through all possible models, so the number of models goes as 2^#covs, where #covs=number of covariates. Will take extraordinary long time for #covs>10. (This also depends on whether randome ffects are used).
2) regress.functioncriterion.dredge: Same as regress.ic.dredge except with a  user-specified criteria. (Was used for training set / test set evaluation.)
3) regress.ic.search: Step-wise traversal of models, typically starting at the null model (no covariates used). Starting model can instead be specified with the option 'start.status' (an array of ones and zeros specifying which covariates are being used in the stating model). Each step tries adding one unused covariate to the previous model, removing one covariate from the previous model and replacing one used with one unused covariate. The best model from all these models is then the starting point for the next step. Contineues until no better model is found.
4) regress.functioncriterion.search: Same as regress.ic.search except with a user-specified criteria.
5) regress.ic.search.2step:  Two-step-wise traversal of models, typically starting at the null model (no covariates used). Starting model can instead be specified with the option 'start.status' (an array of ones and zeros specifying which covariates are being used in the stating model). Each step tries adding one or two unused covariate to the previous model, removing one one or two covariate from the previous model and replacing one or two used with one or two unused covariate. The best model from all these models is then the starting point for the next step. Contineues until no better model is found.
6) regress.random.ic.search: Random regresison model search. For each iteration, tries a given number (num.tries.per.iteration=10000 per default) of randomly picked models. If there is one deemed better than the one from the previous iteration, the function continues, otherwise it stops.
7) regress.random.functioncriterion.search: Same as regress.random.ic.search except with a user-specified criteria.
8) regress.mcmc.ic.search: Uses MCMC to jump between models (1-5 replacements per propoes new model), so aas to allows for leaping between models that differ a bit. Uses the information criteria as 'model probability' in the form exp(-IC/2).
9) regress.mcmc.functioncriterion.search: Same as regress.mcmc.ic.search except with a user-specified criteria.
10) regress.hypo.search: Hypothesis testing based search. Starts at the null model and works itself up (choosing the covariate that is "most" significant, i.e. lowest p-value) until no more statistically significant covariates are found. If bonferroni=TRUE, uses Bonferroni-correction in each step.


