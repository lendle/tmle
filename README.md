Frome DESCRIPTION:
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


From comment in tmle.R:
Targeted Maximum Likelihood Estimation
for non-parametric estimation of the marginal effect of a binary point
treatment, adjusting for treatment (g) and missingness (g.Delta) mechanisms
Parameters include: additive treatment effect: E_W[E(Y|A=1,W) - E(Y|A=0,W)]
and, for binary outcomes, relative risk(RR) and odds ratio(RR)
  mu1 = E_W[E(Y|A=1,W)], mu0 = E_W[E(Y|A=0,W)]
  RR = mu1/mu0
  OR =  mu1/(1-mu1)/(mu0 / (1-mu0))
Controlled direct effect estimation available for optional binary Z intermediate variable,
(P(DeltaY|Z,A,W, Delta=1)P(Delta|Z,A,W)P(Z|A,W)P(A|W)P(W))
EY1 parameter estimated when there is missinginess and no treatment assignment
author: Susan Gruber, sgruber@berkeley.edu
date:   September 25, 2010
revised: December 20, 2010

Models or estimated values for Q, g, g.Z, g.Delta can be user-supplied or 
estimated using super learner.  Cross-validated inital Q can be obtained
by specifying cvQinit=TRUE
