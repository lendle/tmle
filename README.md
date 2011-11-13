tmle implements targeted maximum likelihood estimation,
first described in van der Laan and Rubin, 2006 (Targeted
Maximum Likelihood Learning, The International Journal of
biostatistics, 2(1), 1006.  This implementation calculates the
adjusted marginal difference in mean outcome associated with a
binary point treatment, for continuous or binary outcomes.
Relative risk and odds ratio estimates are also reported for
binary outcomes. Missingness in the outcome is allowed, but not
in treatment assignment or baseline covariate values. Effect
estimation stratified by a binary mediating variable is also
available. The population mean is calculated when there is
missingness, and no variation in the treatment assignment.  An
ID argument can be used to identify repeated measures. Default
settings call SuperLearner to estimate the Q and g portions of
the likelihood, unless values or a user-supplied regression
function are passed in as arguments.