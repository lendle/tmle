#----------------------- calcParameters -----------------------
# Calculate parameter estimates, CIs and p-values
# arguments:
#   Y - outcome
#   A - binary treatment indicator
#   I.Z - Indicator that Z=z (needed for CDE estimation)
#   Delta - missingness indicator
# 	g1W - values of P(Delta=1|Z,A,W)*P(Z=z|A,W)*P(A=1|W)
# 	g0W - values of P(Delta=1|Z,A,W)*P(Z=z|A,W)*P(A=0|W)
#   Delta =- missingness indicator 
#   Q - nx3 matrix(Q(A,W), Q(1,W), Q(0,W))
#   mu1 - targeted estimate of EY1
#   mu0 - targeted estimate of EY0
#   id - subject id
#   family - "gaussian" or "binomial" (or "poisson") 
# returns:
#   list: (EY1, ATE, RR, OR) 
#	 EY1 - population mean outcome (NULL unless no variation in A)
#	 ATE - additive effect: psi, CI, pvalue
#	 RR - relative risk (NULL if family != binomial): psi, CI, pvalue, log.psi, var.log.psi
#	 OR - odds ratio (NULL if family != binomial): psi, CI, pvalue, log.psi, var.log.psi
#--------------------------------------------------------------

##' Calculate Parameter Estimates (calcParameters)
##' 
##' An internal function called by the \code{tmle} function to calculate the
##' population mean effect when there is missingness in the data, but no
##' treatment assignment.  When observations are in treatment and control
##' groups, estimates the additive treatment effect, and if the outcome is
##' binary, also the relative risk and odds ratio parameters.  P-values and 95%
##' confidence intervals are also calculated (on the log scale for RR and OR).
##' 
##' 
##' @param Y continuous or binary outcome variable
##' @param A binary treatment indicator, \code{1} - treatment, \code{0} -
##' control
##' @param I.Z Indicator Z=z, needed for CDE estimation
##' @param Delta indicator of missing outcome. \code{1} - observed, \code{0} -
##' missing
##' @param g1W censoring mechanism estimates, \eqn{P(A=1|W)*P(Delta=1|A,W)}
##' @param g0W censoring mechanism estimates, \eqn{P(A=0|W)*P(Delta=1|A,W)}
##' @param Q a 3-column matrix \code{(Q(A,W), Q(1,W), Q(0,W))}
##' @param mu1 targeted estimate of \eqn{E(Y|A=1,W)}
##' @param mu0 targeted estimate of \eqn{E(Y|A=0,W)}
##' @param id subject identifier
##' @param family family specification for regressions, generally
##' \sQuote{gaussian} for continuous oucomes, \sQuote{binomial} for binary
##' outcomes
##' @return \item{EY1}{Population mean outcome estimate, variance, p-value, 95%
##' confidence interval (missingness only, no treatment assingment), or
##' \code{NULL}} \item{ATE}{additive treatment effect estimate, variance,
##' p-value, 95% confidence interval, or \code{NULL}} \item{RR}{relative risk
##' estimate, p-value, 95% confidence interval, log(RR), variance(log(RR)), or
##' \code{NULL}} \item{OR}{odds ratio estimate, p-value, 95% confidence
##' interval, log(OR), variance(log(OR)), or \code{NULL}}
##' @author Susan Gruber
##' @seealso \code{\link{tmle}}, \code{\link{estimateQ}},
##' \code{\link{estimateG}}
##' @export
calcParameters <- function(Y,A, I.Z, Delta, g1W, g0W, Q, mu1, mu0, id, family){
	n.id <- length(unique(id))
	Y[is.na(Y)] <- 0
	EY1 <- ATE <- RR <- OR <- NULL
	if(length(unique(A))==1){
		pDelta1 <- A*g1W+ (1-A)*g0W   # one of these will be zero, the other pDelta1
		IC.EY1 <- Delta/pDelta1*(Y-Q[,"QAW"]) + Q[,"Q1W"] - mu1
		if(n.id < length(id)){
			IC.EY1 <- as.vector(by(IC.EY1, id, mean))
		}
		EY1$psi <- mu1
		EY1$var.psi <- var(IC.EY1)/n.id
		EY1$CI <- EY1$psi + c(-1.96,1.96)*sqrt(EY1$var.psi)
		EY1$pvalue <- 2*pnorm(-abs(EY1$psi/sqrt(EY1$var.psi)))
	} else {
		IC.ATE<- I.Z*(A/g1W - (1-A)/g0W)*Delta*(Y-Q[,"QAW"]) + Q[,"Q1W"] - Q[,"Q0W"] - (mu1-mu0)
		if(n.id < length(id)){
			IC.ATE <- as.vector(by(IC.ATE, id, mean))
		}
		ATE$psi <- mu1-mu0
		ATE$var.psi <- var(IC.ATE)/n.id
		ATE$CI <- ATE$psi + c(-1.96,1.96)*sqrt(ATE$var.psi)
		ATE$pvalue <- 2*pnorm(-abs(ATE$psi/sqrt(ATE$var.psi)))

		if(family=="binomial"){		
			IC.logRR <- 1/mu1*(I.Z*(A/g1W)*Delta*(Y-Q[,"QAW"]) + Q[,"Q1W"] - mu1) - 
					1/mu0*(I.Z*(1-A)/g0W*Delta*(Y-Q[,"QAW"])+ Q[,"Q0W"] - mu0)
			if(n.id < length(id)){
				IC.logRR <- as.vector(by(IC.logRR, id, mean))
			}
			var.psi.logRR <- var(IC.logRR)/n.id
			RR$psi <- mu1/mu0
			RR$CI  <- exp(log(RR$psi) + c(-1.96,1.96)*sqrt(var.psi.logRR))
			RR$pvalue <- 2*pnorm(-abs(log(RR$psi)/sqrt(var.psi.logRR)))
			RR$log.psi <- log(RR$psi)
			RR$var.log.psi <- var.psi.logRR

			IC.logOR <- 1/(mu1*(1-mu1)) * (I.Z*A/g1W*Delta*(Y-Q[,"QAW"]) + Q[,"Q1W"]) - 
					1/(mu0*(1-mu0)) * (I.Z*(1-A)/g0W*Delta * (Y-Q[,"QAW"]) + Q[,"Q0W"])
			if(n.id < length(id)){
				IC.logOR <- as.vector(by(IC.logOR, id, mean))
			}
			var.psi.logOR <- var(IC.logOR)/n.id
			OR$psi <-  mu1/(1-mu1)/(mu0 / (1-mu0))
			OR$CI  <- exp(log(OR$psi) + c(-1.96,1.96)*sqrt(var.psi.logOR))
			OR$pvalue <- 2*pnorm(-abs(log(OR$psi)/sqrt(var.psi.logOR)))
			OR$log.psi <- log(OR$psi)
			OR$var.log.psi <- var.psi.logOR
		}
	}
	return(list(EY1=EY1, ATE=ATE, RR=RR, OR=OR))
}
