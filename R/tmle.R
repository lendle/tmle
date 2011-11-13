##' Summarization of the results of a call to the tmle routine
##' 
##' These functions are all \link{methods} for class \code{tmle},
##' \code{tmle.list}, \code{summary.tmle}, \code{summary.tmle.list} objects
##' 
##' \code{print.tmle} prints the estimate, variance, p-value, and 95%
##' confidence interval only.  \code{print.summary.tmle}, called indirectly by
##' entering the command \kbd{summary(result)} (where \code{result} has class
##' \code{tmle}), outputs additional information.  Controlled direct effect
##' estimates have class \code{tmle.list}, a list of two objects of class
##' \code{tmle}.  The first item corresponds to \eqn{Z=0}, the second to
##' \eqn{Z=1}
##' 
##' @aliases summary.tmle print.summary.tmle print.tmle summary.tmle.list
##' print.summary.tmle.list print.tmle.list
##' @param object an object of class \code{tmle} or \code{tmle.list}.
##' @param x an object of class \code{tmle} or \code{tmle.list} for summary
##' functions, class \code{summary.tmle} or \code{summary.tmle.list} for print
##' functions.
##' @param \dots currently ignored.
##' @return \item{estimates}{list of parameter estimates, pvalues, and 95%
##' confidence intervals} \item{Qmodel}{working model used to obtain initial
##' estimate of \code{Q} portion of the likelihood, if \code{glm} used}
##' \item{Qterms}{terms in the model for \code{Q}} \item{Qcoef}{coefficient of
##' each term in model for \code{Q}} \item{gmodel}{model used to estimate
##' treatment mechanism \code{g}} \item{gterms}{terms in the treatment
##' mechanism model} \item{gcoef}{coefficient of each term in model for
##' treatment mechanism} \item{gtype}{description of estimation procedure for
##' treatment mechanism, e.g. "glm"} \item{g.Zmodel}{model used to estimate
##' intermediate variable assignment mechanism \code{g.Z}}
##' \item{g.Zterms}{terms in the intermediate mechanism model}
##' \item{g.Zcoef}{coefficient of each term in model for intermediate
##' mechanism} \item{g.Ztype}{description of estimation procedure for
##' intermediate variable} \item{g.Deltamodel}{model used to estimate
##' missingness mechanism \code{g.Delta}} \item{g.Deltaterms}{terms in the
##' missingness mechanism model} \item{g.Deltacoef}{coefficient of each term in
##' model for missingness mechanism} \item{g.Deltatype}{description of
##' estimation procedure for missingness}
##' @author Susan Gruber
##' @seealso \code{\link{tmle}}
##' @examples
##' 
##' # generate data
##'   n <- 500
##'   W <- matrix(rnorm(n*3), ncol=3)
##'   A <- rbinom(n,1, 1/(1+exp(-(.1*W[,1] - .1*W[,2] + .5*W[,3]))))
##'   Y <- A + 2*W[,1] + W[,3] + W[,2]^2 + rnorm(n)
##'   colnames(W) <- paste("W",1:3, sep="")
##' 
##'   result <- tmle(Y,A,W, Qform="Y~A+W1", g1W=rep(.5, n))
##'   summary(result)
##'
##' @export
summary.tmle <- function(object,...) {
  if(identical(class(object), "tmle")){
		Qmodel <- Qcoef <- gmodel <- gcoef <- g.Deltamodel <- g.Deltacoef <- NULL
		g.Zmodel <- g.Zcoef <- NULL
		Qterms <- gterms <- g.Deltaterms <- g.Zterms <- ""
		if(!is.null(object$Qinit)){
		  if (!is.null(object$Qinit$coef)) {
			Qcoef <- object$Qinit$coef	
			if(class(Qcoef) == "matrix"){
				Qterms <- colnames(Qcoef)
			} else {
				Qterms <- names(Qcoef)
			}
			Qmodel <- paste("Y ~ 1")
			if(length(Qterms) > 1) {
				Qmodel <- paste("Y ~ ", paste(Qterms, collapse =" + "))
			} 
		  } else {
		  	Qmodel <- object$Qinit$type
		  }
		}
		if(!is.null(object$g)){
		if (!is.null(object$g$coef)) {
			gbd <- object$g$bound
			gcoef <- object$g$coef
			if(class(gcoef) == "matrix"){
				gterms <- colnames(gcoef)
			} else {
				gterms <- names(gcoef)
			}			
			gmodel <- paste("A ~ 1")
			if(length(gterms) > 1) {
				gmodel <- paste("A ~ ", paste(gterms, collapse =" + "))
			} 
		}}
		if(!is.null(object$g.Z)){
		if (!is.null(object$g.Z$coef)) {
			g.Zcoef <- object$g.Z$coef
			if(class(g.Zcoef) == "matrix"){
				g.Zterms <- colnames(g.Zcoef)
			} else {
				g.Zterms <- names(g.Zcoef)
			}			
			g.Zmodel <- paste("Z ~ 1")
			if(length(g.Zterms) > 1) {
				g.Zmodel <- paste("Z ~", paste(g.Zterms, collapse =" + "))
			} 
		}}
		if(!is.null(object$g.Delta)){
		if (!is.null(object$g.Delta$coef)) {
			g.Deltacoef <- object$g.Delta$coef
			if(class(g.Deltacoef) == "matrix"){
				g.Deltaterms <- colnames(g.Deltacoef)
			} else {
				g.Deltaterms <- names(g.Deltacoef)
			}			
			g.Deltamodel <- paste("Delta ~ 1")
			if(length(g.Deltaterms) > 1) {
				g.Deltamodel <- paste("Delta ~", paste(g.Deltaterms, collapse =" + "))
			} 
		}}
		summary.tmle <- list(estimates=object$estimates,
						Qmodel=Qmodel, Qterms=Qterms, Qcoef=Qcoef, Qtype=object$Qinit$type,
						gbd=gbd, gmodel=gmodel, gterms=gterms, gcoef=gcoef,gtype=object$g$type, 
						g.Zmodel=g.Zmodel, g.Zterms=g.Zterms, g.Zcoef=g.Zcoef, g.Ztype=object$g.Z$type,
						g.Deltamodel=g.Deltamodel, g.Deltaterms=g.Deltaterms, g.Deltacoef=g.Deltacoef, g.Deltatype=object$g.Delta$type )
		class(summary.tmle) <- "summary.tmle"
	} else {
		stop("object must have class 'tmle'")
		summary.tmle <- NULL
	}
	return(summary.tmle)
}

#-------------print.summary.tmle------------------
# print tmle summary object
#-------------------------------------------------
print.summary.tmle <- function(x,...) {
  if(identical(class(x), "summary.tmle")){
  	cat(" Initial estimation of Q\n")
  	cat("\t Procedure:", x$Qtype)
  	if(!(is.na(x$Qcoef[1]))){
  		cat("\n\t Model:\n\t\t", x$Qmodel)
  		cat("\n\n\t Coefficients: \n")
  		terms <- sprintf("%15s", x$Qterms)
  		extra <- ifelse(x$Qcoef >= 0, "  ", " ")
  		for(i in 1:length(x$Qcoef)){
  			cat("\t", terms[i], extra[i], x$Qcoef[i], "\n")
  		}
  	}
  	cat("\n Estimation of g (treatment mechanism)\n")
  	cat("\t Procedure:", x$gtype,"\n")
  	if(!(is.na(x$gcoef[1]))){	
 		cat("\t Model:\n\t\t", x$gmodel, "\n")
  		cat("\n\t Coefficients: \n")
  		terms <- sprintf("%15s", x$gterms)
  		extra <- ifelse(x$gcoef >= 0, "  ", " ")
  		for(i in 1:length(x$gcoef)){
  			cat("\t", terms[i], extra[i], x$gcoef[i], "\n")
  	}}
  	cat("\n Estimation of g.Z (intermediate variable assignment mechanism)\n")
  	cat("\t Procedure:", x$g.Ztype, "\n")
  	if(!(is.na(x$g.Zcoef[1]))){		
  		cat("\t Model:\n\t\t",x$g.Zmodel, "\n")
  		cat("\n\t Coefficients: \n")
  		terms <- sprintf("%15s", x$g.Zterms)
  		extra <- ifelse(x$g.Zcoef >= 0, "  ", " ")
  		for(i in 1:length(x$g.Zcoef)){
  			cat("\t", terms[i], extra[i], x$g.Zcoef[i], "\n")
  	}}
  	cat("\n Estimation of g.Delta (missingness mechanism)\n")
  	cat("\t Procedure:", x$g.Deltatype, "\n")
  	if(!(is.na(x$g.Deltacoef[1]))){		
  		cat("\t Model:\n\t\t",x$g.Deltamodel, "\n")
  		cat("\n\t Coefficients: \n")
  		terms <- sprintf("%15s", x$g.Deltaterms)
  		extra <- ifelse(x$g.Deltacoef >= 0, "  ", " ")
  		for(i in 1:length(x$g.Deltacoef)){
  			cat("\t", terms[i], extra[i], x$g.Deltacoef[i], "\n")
  	}}
  	cat("\n Bounds on g: (", x$gbd,")\n")
	if(!is.null(x$estimates$EY1)){
			cat("\n Population Mean")
			cat("\n   Parameter Estimate: ", signif(x$estimates$EY1$psi,5))
			cat("\n   Estimated Variance: ", signif(x$estimates$EY1$var.psi,5))
			cat("\n              p-value: ", ifelse(x$estimates$EY1$pvalue <= 2*10^-16, "<2e-16",signif(x$estimates$EY1$pvalue,5)))
			cat("\n    95% Conf Interval:",paste("(", signif(x$estimates$EY1$CI[1],5), ", ", signif(x$estimates$EY1$CI[2],5), ")", sep=""),"\n") 
		}

	if(!is.null(x$estimates$ATE)){
			cat("\n Additive Effect")
			cat("\n   Parameter Estimate: ", signif(x$estimates$ATE$psi,5))
			cat("\n   Estimated Variance: ", signif(x$estimates$ATE$var.psi,5))
			cat("\n              p-value: ", ifelse(x$estimates$ATE$pvalue <= 2*10^-16, "<2e-16",signif(x$estimates$ATE$pvalue,5)))
			cat("\n    95% Conf Interval:",paste("(", signif(x$estimates$ATE$CI[1],5), ", ", signif(x$estimates$ATE$CI[2],5), ")", sep=""),"\n") 
		}
		if(!is.null(x$estimates$RR)){
			cat("\n Relative Risk")
			cat("\n   Parameter Estimate: ", signif(x$estimates$RR$psi,5))
			cat("\n              p-value: ", ifelse(x$estimates$RR$pvalue <= 2*10^-16, "<2e-16",signif(x$estimates$RR$pvalue,5)))
			cat("\n    95% Conf Interval:",paste("(", signif(x$estimates$RR$CI[1],5), ", ", signif(x$estimates$RR$CI[2],5), ")", sep=""),"\n") 
			cat("\n              log(RR): ", signif(x$estimates$RR$log.psi,5))
			cat("\n    variance(log(RR)): ", signif(x$estimates$RR$var.log.psi,5), "\n")

		}
		if(!is.null(x$estimates$OR)){
			cat("\n Odds Ratio")
			cat("\n   Parameter Estimate: ", signif(x$estimates$OR$psi,5))
			cat("\n              p-value: ", ifelse(x$estimates$OR$pvalue <= 2*10^-16, "<2e-16",signif(x$estimates$OR$pvalue,5)))
			cat("\n    95% Conf Interval:",paste("(", signif(x$estimates$OR$CI[1],5), ", ", signif(x$estimates$OR$CI[2],5), ")", sep=""),"\n") 
			cat("\n              log(OR): ", signif(x$estimates$OR$log.psi,5))
			cat("\n    variance(log(OR)): ", signif(x$estimates$OR$var.log.psi,5),"\n")

		}
 
  } else {
  	 stop("Error calling print.summary.tmle. 'x' needs to have class 'summary.tmle'\n")
  }
}

#-------------print.tmle------------------
# print object returned by tmle function
#-----------------------------------------
print.tmle <- function(x,...) {
	if(identical(class(x), "tmle")){
		if(!is.null(x$estimates$EY1)){
			cat(" Population Mean")
			cat("\n   Parameter Estimate: ", signif(x$estimates$EY1$psi,5))
			cat("\n   Estimated Variance: ", signif(x$estimates$EY1$var.psi,5))
			cat("\n              p-value: ", ifelse(x$estimates$EY1$pvalue <= 2*10^-16, "<2e-16",signif(x$estimates$EY1$pvalue,5)))
			cat("\n    95% Conf Interval:",paste("(", signif(x$estimates$EY1$CI[1],5), ", ", signif(x$estimates$EY1$CI[2],5), ")", sep=""),"\n") 
		}

		if(!is.null(x$estimates$ATE)){
			cat(" Additive Effect")
			cat("\n   Parameter Estimate: ", signif(x$estimates$ATE$psi,5))
			cat("\n   Estimated Variance: ", signif(x$estimates$ATE$var.psi,5))
			cat("\n              p-value: ", ifelse(x$estimates$ATE$pvalue <= 2*10^-16, "<2e-16",signif(x$estimates$ATE$pvalue,5)))
			cat("\n    95% Conf Interval:",paste("(", signif(x$estimates$ATE$CI[1],5), ", ", signif(x$estimates$ATE$CI[2],5), ")", sep=""),"\n") 
		}
		if(!is.null(x$estimates$RR)){
			cat("\n Relative Risk")
			cat("\n   Parameter Estimate: ", signif(x$estimates$RR$psi,5))
			cat("\n              p-value: ", ifelse(x$estimates$RR$pvalue <= 2*10^-16, "<2e-16",signif(x$estimates$RR$pvalue,5)))
			cat("\n    95% Conf Interval:",paste("(", signif(x$estimates$RR$CI[1],5), ", ", signif(x$estimates$RR$CI[2],5), ")", sep=""),"\n") 
			cat("\n              log(RR): ", signif(x$estimates$RR$log.psi,5))
			cat("\n    variance(log(RR)): ", signif(x$estimates$RR$var.log.psi,5),"\n")

		}
		if(!is.null(x$estimates$OR)){
			cat("\n Odds Ratio")
			cat("\n   Parameter Estimate: ", signif(x$estimates$OR$psi,5))
			cat("\n              p-value: ", ifelse(x$estimates$OR$pvalue <= 2*10^-16, "<2e-16",signif(x$estimates$OR$pvalue,5)))
			cat("\n    95% Conf Interval:",paste("(", signif(x$estimates$OR$CI[1],5), ", ", signif(x$estimates$OR$CI[2],5), ")", sep=""),"\n") 
			cat("\n              log(OR): ", signif(x$estimates$OR$log.psi,5))
			cat("\n    variance(log(OR)): ", signif(x$estimates$OR$var.log.psi,5),"\n")

		}
	} else {
 		stop("Error calling print.tmle. 'x' needs to have class 'tmle'\n")
 }}

#-------------print.tmle.list------------------
# print object returned by tmle 
# when there is a controlled direct effect
#-----------------------------------------
print.tmle.list <- function(x,...) {
	cat("Controlled Direct Effect\n")
	cat("           ----- Z = 0 -----\n")
	print(x[[1]])
	cat("\n           ----- Z = 1 -----\n")
	print(x[[2]])

}

#-------------summary.tmle.list------------------
# create tmle summary object for controlled direct 
# effect
#---------------------------------------	 
summary.tmle.list <- function(object,...) {
	summary.tmle.list <- list(Z0=summary(object[[1]]), Z1=summary(object[[2]]))
	class(summary.tmle.list) <- "summary.tmle.list"
	return(summary.tmle.list)
}

#-------------print.summary.tmle.list------------------
# create tmle summary object for controlled direct 
# effect
#---------------------------------------	 
print.summary.tmle.list <- function(x,...) {
	cat("Controlled Direct Effect\n")
	cat("           ----- Z = 0 -----\n")
	print(x[[1]])
	cat("\n           ----- Z = 1 -----\n")
	print(x[[2]])
}








#-------------------------------tmle----------------------------------------
# estimate marginal treatment effect for binary point treatment
# accounting for missing outcomes.
# EY1 parameter if no variation in A
# arguments:
# Y - outcome
# A - binary treatment indicator, 1-treatment, 0 - control
# W - vector, matrix or dataframe containing baseline covariates
# Z - optional binary intermediate between A and Y - if specified
# 	  we'll compute "controlled direct effect" to get parameter estimate at
#    each value of Z 
# Delta - indicator of missingness.  1 - observed, 0 - missing
# Q - E(Y|Z,A,W), optional nx2 matrix [E(Y|A=0,W), E(Y|A=1,W)] (if CDE estimation, Z=0)
# Q.Z1 - optional nx2 matrix [E(Y|A=0,W), E(Y|A=1,W)] (if CDE estimation, when Z=1)
# g1W - optional values for P(A=1|W)
# gform - optional glm regression formula 
# gbound - one value for symmetrical bounds on g1W, or a vector containing upper and lower bounds
# pZ1 - optional values for P(Z=1|A,W)
# g.Zform - optional glm regression formula
# pDelta1 - optional values for P(Delta=1|Z,A,W)
# g.Deltaform - optional glm regression formula  
# Q.SL.library- optional Super Learner library for estimation of Q, 
# cvQinit - if TRUE obtain cross-validated initial Q
# g.SL.library - optional library for estimation of g and g.Delta
# family - family specification for regression models, defaults to gaussian
# fluctuation - "logistic" (default) or "linear" (for targeting step)
# alpha - bound on predicted probabilities for Q (0.005, 0.995 default)
# id - optional subject identifier
# verbose - flag for controlling printing of messages
#-------------------------------------------------------------------------------






##' Targeted Maximum Likelihood Estimation
##' 
##' Carries out targeted maximum likelihood estimation of marginal additive
##' treatment effect of a binary point treatment on an outcome.  For binary
##' outcomes risk ratio and odds ratio estimates are also reported. The
##' \code{tmle} function is minimally called with arguments \code{(Y,A,W)},
##' where \code{Y} is a continuous or binary outcome variable, \code{A} is a
##' binary treatment variable, (\code{A=1} for treatment, \code{A=0} for
##' control), and \code{W} is a matrix or dataframe of baseline covariates. The
##' population mean outcome is calculated when there is no variation in
##' \code{A}. If values of binary mediating variable \code{Z} are supplied,
##' estimates are returned at each level of \code{Z}. Missingness in the
##' outcome is accounted for in the estimation procedure if missingness
##' indicator \code{Delta} is 0 for some observations.  Repeated measures can
##' be identified using the \code{id} argument.
##' 
##' \code{gbounds} defaults to (0.025, 0.975).  To specify symmetrical
##' truncation levels, only one value has to be provided. The upper bound is
##' set to 1 for EY1 parameter estimated when there is missingness only, no
##' treatment assignment.
##' 
##' Controlled direct effects are estimated when binary covariate \code{Z} is
##' non-null. The tmle function returns an object of class \code{tmle.list}, a
##' list of two items of class \code{tmle}.  The first corresponds to estimates
##' obtained when \code{Z} is fixed at \eqn{0}, the second correspondes to
##' estimates obtained when \code{Z} is fixed at \eqn{1}.
##' 
##' \code{Q.SL.library} defaults to (\sQuote{SL.glm}, \sQuote{SL.step},
##' \sQuote{SL.DSA.2}) when family=\sQuote{gaussian}, (\sQuote{SL.glm},
##' \sQuote{SL.step}, \sQuote{SL.knn}, \sQuote{SL.DSA.2}) when family =
##' \sQuote{binomial}
##' 
##' \code{g.SL.library} Defaults to (\sQuote{SL.glm}, \sQuote{SL.step},
##' \sQuote{SL.knn}, \sQuote{SL.DSA.2})
##' 
##' See \code{SuperLearner} help files for further information.
##' 
##' @param Y continuous or binary outcome variable
##' @param A binary treatment indicator, \code{1} - treatment, \code{0} -
##' control
##' @param W vector, matrix, or dataframe containing baseline covariates
##' @param Z optional binary indicator for intermediate covariate for conrolled
##' direct effect estimation
##' @param Delta indicator of missing outcome or treatment assignment.
##' \code{1} - observed, \code{0} - missing
##' @param Q \eqn{E(Y|A,W)}, optional \eqn{nx2} matrix of initial values for
##' \eqn{Q} portion of the likelihood, \eqn{(E(Y|A=0,W), E(Y|A=1,W))}
##' @param Q.Z1 \eqn{E(Y|Z=1,A,W)}, optional \eqn{nx2} matrix of initial values
##' for \eqn{Q} portion of the likelihood, \eqn{(E(Y|Z=1,A=0,W),
##' E(Y|Z=1,A=1,W))}. (When specified, values for \eqn{E(Y|Z=0,A=0,W),
##' E(Y|Z=0,A=1,W)} are passed in using the \code{Q} argument
##' @param Qform optional regression formula for estimation of \eqn{E(Y|A, W)},
##' suitable for call to \code{glm}
##' @param Qbounds vector of upper and lower bounds on \code{Y} and predicted
##' values for initial \code{Q}
##' @param Q.SL.library optional vector of prediction algorithms to use for
##' \code{SuperLearner} estimation of initial \code{Q}
##' @param cvQinit logical, if \code{TRUE}, estimates cross-validated predicted
##' values using discrete super learning, default=\code{FALSE}
##' @param g1W optional vector of conditional treatment assingment
##' probabilities, \eqn{P(A=1|W)}
##' @param gform optional regression formula of the form \code{A~W}, if
##' specified this overrides the call to \code{SuperLearner}
##' @param gbound value between (0,1) for truncation of predicted
##' probabilities. See \code{Details} section for more information
##' @param pZ1 optional\eqn{nx2} matrix of conditional probabilities
##' \eqn{P(Z=1|A=0,W), P(Z=1|A=1,W)}
##' @param g.Zform optionalregression formula of the form \code{Z~A+W}, if
##' specified this overrides the call to \code{SuperLearner}
##' @param pDelta1 optional matrix of conditional probabilities for missingness
##' mechanism, \eqn{nx2} when \code{Z} is \code{NULL} \eqn{P(Delta=1|A=0,W),
##' P(Delta=1|A=1,W)}.  \eqn{nx4} otherwise, \eqn{P(Delta=1|Z=0,A=0,W),
##' P(Delta=1|Z=0,A=1,W), P(Delta=1|Z=1,A=0,W), P(Delta=1|Z=1,A=1,W)}
##' @param g.Deltaform optional regression formula of the form
##' \code{Delta~A+W}, if specified this overrides the call to
##' \code{SuperLearner}
##' @param g.SL.library optional vector of prediction algorithms to use for
##' \code{SuperLearner} estimation of \code{g1W} or \code{pDelta1}
##' @param family family specification for working regression models, generally
##' \sQuote{gaussian} for continuous outcomes (default), \sQuote{binomial} for
##' binary outcomes
##' @param fluctuation \sQuote{logistic} (default), or \sQuote{linear}
##' @param alpha used to keep predicted initial values bounded away from (0,1)
##' for logistic fluctuation
##' @param id optional subject identifier
##' @param verbose status messages printed if set to \code{TRUE}
##' (default=\code{FALSE})
##' @return \item{estimates}{list with elements EY1 (population mean), ATE
##' (additive treatment effect), RR (relative risk), OR (odds ratio). Each
##' element in the estimates of these is itself a list containing \itemize{
##' \item psi - parameter estimate \item pvalue - two-sided p-value \item CI -
##' 95% confidence interval \item var.psi - Influence-curve based variance of
##' estimate (ATE parameter only) \item log.psi - Parameter estimate on log
##' scale (RR and OR parameters) \item var.log.psi - Influence-curve based
##' variance of estimate on log scale (RR and OR parameters) }}
##' \item{Qinit}{initial estimate of \code{Q}. \code{Qinit$coef} are the
##' coefficients for a \code{glm} model for \code{Q}, if applicable.
##' \code{Qinit$Q} is an \eqn{nx2} matrix, where \code{n} is the number of
##' observations.  Columns contain predicted values for \code{Q(0,W),Q(1,W)}
##' using the initial fit.  \code{Qinit$type} is method for estimating
##' \code{Q}} \item{Qstar}{targeted estimate of \code{Q}, an \eqn{nx2} matrix
##' with predicted values for \code{Q(0,W), Q(1,W)} using the updated fit}
##' \item{g}{treatment mechanism estimate. A list with three items:
##' \code{g$g1W} contains estimates of \eqn{P(A=1|W)} for each observation,
##' \code{g$coef} the coefficients for the model for \eqn{g} when \code{glm}
##' used, \code{g$type} estimation procedure} \item{g.Z}{intermediate covariate
##' assignment estimate (when applicable). A list with three items:
##' \code{g.Z$g1W} an \eqn{nx2} matrix containing values of \eqn{P(Z=1|A=1,W),
##' P(Z=1|A=0,W)} for each observation, \code{g.Z$coef} the coefficients for
##' the model for \eqn{g} when \code{glm} used, \code{g.Z$type} estimation
##' procedure} \item{g.Delta}{missingness mechanism estimate. A list with three
##' items: \code{g.Delta$g1W} an \eqn{nx4} matrix containing values of
##' \eqn{P(Delta=1|Z,A,W)} for each observation, with (Z=0,A=0), (Z=0,A=1),
##' (Z=1,A=0),(Z=1,A=1). (When \code{Z} is \code{NULL}, columns 3 and 4 are
##' duplicates of 1 and 2.) \code{g.Delta$coef} the coefficients for the model
##' for \eqn{g} when \code{glm} used, \code{g.Delta$type} estimation procedure}
##' @author Susan Gruber \email{sgruber@@berkeley.edu}, in collaboration with
##' Mark van der Laan.
##' @seealso \code{\link{summary.tmle}}, \code{\link{estimateQ}},
##' \code{\link{estimateG}}, \code{\link{calcParameters}}
##' @references 1. Gruber, S. and van der Laan, M.J. (2009), Targeted Maximum
##' Likelihood Estimation: A Gentle Introduction. \emph{U.C. Berkeley Division
##' of Biostatistics Working Paper Series}.  Working Paper 252.
##' \url{http://www.bepress.com/ucbbiostat/paper252}
##' 
##' 2. Gruber, S. and van der Laan, M.J.  (2010), A Targeted Maximum Likelihood
##' Estimator of a Causal Effect on a Bounded Continuous Outcome. \emph{The
##' International Journal of Biostatistics}, 6(1), 2010.
##' 
##' 3. van der Laan, M.J. and Rubin, D. (2006), Targeted Maximum Likelihood
##' Learning. \emph{The International Journal of Biostatistics}, 2(1).
##' \url{http://www.bepress.com/ijb/vol2/iss1/11/}
##' 
##' 4. van der Laan, M.J., Rose, S., and Gruber,S., editors, (2009) Readings in
##' Targeted Maximum Likelihood Estimation . \emph{U.C. Berkeley Division of
##' Biostatistics Working Paper Series}.  Working Paper 254.
##' \url{http://www.bepress.com/ucbbiostat/paper254}
##' @examples
##' 
##' library(tmle)
##' 
##' # generate data
##'   set.seed(5)
##'   n <- 500
##'   W <- matrix(rnorm(n*3), ncol=3)
##'   A <- rbinom(n,1, 1/(1+exp(-(.1*W[,1] - .1*W[,2] + .5*W[,3]))))
##'   Y <- A + 2*W[,1] + W[,3] + W[,2]^2 + rnorm(n)
##'   colnames(W) <- paste("W",1:3, sep="")
##' 
##' # Example 1. Simplest function invocation 
##' # If available, SuperLearner called to estimate Q, g
##' # otherwise, main terms regression using glm
##' # Delta defaults to 1 for all observations   
##'   result1 <- tmle(Y,A,W)
##'   summary(result1)
##' 
##' # Example 2: 
##' # User-supplied regression formulas to estimate Q and g
##' # binary outcome
##'   n <- 250
##'   W <- matrix(rnorm(n*3), ncol=3)
##'   colnames(W) <- paste("W",1:3, sep="")
##'   A <- rbinom(n,1, plogis(0.6*W[,1] +0.4*W[,2] + 0.5*W[,3]))
##'   Y <- rbinom(n,1, plogis(A + 0.2*W[,1] + 0.1*W[,2] + 0.2*W[,3]^2 ))
##'   result2 <- tmle(Y,A,W, family="binomial", Qform=Y~A+W1+W2+W3, gform=A~W1+W2+W3)
##'   summary(result2)
##' 
##' # Example 3: Population mean outcome
##' # User-supplied (misspecified) model for Q, 
##' # If available, Super learner called to estimate g, g.Delta, otherwise glm
##' # approx. 20% missing at random
##'   Delta <- rbinom(n, 1, 1/(1+exp(-(1.7-1*W[,1]))))
##'   result3 <- tmle(Y,A=NULL,W, Delta=Delta, Qform="Y~A+W1+W2+W3")
##'   summary(result3)
##' 
##' # Example 4: Controlled direct effect with missingness
##' # User-supplied models for Q, g, g.Z
##' # User-supplied values for pDelta1 
##'   n <- 1000
##'   W <- matrix(rnorm(n*3), ncol = 3)
##'   colnames(W) <- paste("W", 1:3, sep = "")
##'   A <- rbinom(n,1, plogis(0.6*W[,1] + 0.4*W[,2] + 0.5*W[,3]))
##'   Z <- rbinom(n,1, plogis(0.5 + A))
##'   Y <- A + A*Z+ 0.2*W[,1] + 0.1*W[,2] + 0.2*W[,3]^2 + rnorm(n)
##'   Delta <- rbinom(n,1, plogis(Z + A)) 
##'   pDelta1 <- cbind(rep(plogis(0), n), rep(plogis(1), n),
##'                    rep(plogis(1), n), rep(plogis(2), n))
##'   colnames(pDelta1) <- c("Z0A0", "Z0A1", "Z1A0", "Z1A1")
##'   Y[Delta == 0] <- NA
##'   result4 <- tmle(Y, A, W, Z, Delta = Delta, pDelta1= pDelta1, 
##'                   Qform = Y ~ 1, g.SL.library = "SL.glm") 
##'   result4
##' 
tmle <- function(Y,A,W,Z=NULL, Delta=rep(1,length(Y)),  
				Q=NULL, Q.Z1=NULL, Qform=NULL, Qbounds=NULL, 
				Q.SL.library=NULL, cvQinit=FALSE,
				g1W=NULL, gform=NULL, gbound=0.025, 
				pZ1=NULL, g.Zform=NULL,
				pDelta1=NULL, g.Deltaform=NULL, 
				g.SL.library=c("SL.glm", "SL.step", "SL.knn", "SL.DSA.2"),
				family="gaussian", fluctuation="logistic",
				alpha  = 0.995, id=1:length(Y),verbose=FALSE) {
	# Initializations
	psi.tmle <- varIC <- CI <- pvalue <- NA
	W <- as.matrix(W)
	colnames(W) <- .setColnames(colnames(W), NCOL(W), "W")
	if (identical(family, binomial)) {
        family == "binomial"
    }else if (identical(family, gaussian)) {
        family == "gaussian"
    }else if (identical(family, poisson)) {  
        family == "poisson"
        if(is.null(Qform)){ 		# only glm for Q when family=poisson
        	Qform <- paste("Y~A", paste(colnames(W), collapse="+"), sep="+")
        }
    }
    if(is.null(Q.SL.library)){
   		if(family=="binomial"){
   			Q.SL.library <- c("SL.glm", "SL.step", "SL.knn", "SL.DSA.2")
   		} else {
   			Q.SL.library <- c("SL.glm", "SL.step", "SL.DSA.2")
   		}
   	}

	if(is.null(A) | all(A==0)){
		A <- rep(1, length(Y))
	}
	
	if(!.verifyArgs(Y,Z,A,W,Delta, Qform, gform, g.Zform, g.Deltaform)){
		stop()
	}
   
 	maptoYstar <- fluctuation=="logistic"
 		    	
   	if(!is.null(Z) & !is.null(pZ1)){
   		if(NCOL(pZ1)==2){
   			colnames(pZ1) <- c("A0", "A1")
   		} else {
   			stop("pZ1 must be an nx2 matrix: [P(Z=1|A=0,W),P(Z=1|A=1,W)]\n")
   		}
   	}
   	
  	if(any(Delta!=1) & !is.null(pDelta1)){
   		if(NCOL(pDelta1)==2 & is.null(Z)){
   			pDelta1 <- cbind(pDelta1, pDelta1)
   			colnames(pDelta1) <- c("Z0A0", "Z0A1", "Z1A0", "Z1A1")  # won't use the second set
   		} else if(NCOL(pDelta1)==4){
   			colnames(pDelta1) <- c("Z0A0", "Z0A1", "Z1A0", "Z1A1") 
   		}else {
   			if(is.null(Z)){
   				stop("pDelta1 must be an nx2 matrix: [P(Delta=1|A=0,W), P(Delta=1|A=1,W)]\n")
   			} else {
   				stop("pDelta1 must be an nx4 matrix:\n  [P(Delta=1|Z=0,A=0,W), P(Delta=1|Z=0,A=1,W), P(Delta=1|Z=1,A=0,W), P(Delta=1|Z=1,A=1,W)]\n")
   			}
   		}
   	}
   	   	
   if(is.null(Z)){ 
   		Z <- rep(1, length(Y))	
   	}
   	CDE <- length(unique(Z))>1
	
	# Stage 1	
 	stage1 <- .initStage1(Y, A, Q, Q.Z1, Delta, Qbounds, alpha, maptoYstar, family)		
	Q <- suppressWarnings(estimateQ(Y=stage1$Ystar,Z,A,W, Delta, Q=stage1$Q, Qbounds=stage1$Qbounds, Qform, 
					maptoYstar=maptoYstar, SL.library=Q.SL.library, 
					cvQinit=cvQinit, family=family, id=id, verbose=verbose))
					
	# Stage 2
	if(length(gbound)==1){
		if(length(unique(A))==1 & length(unique(Z))==1){  # EY1 only, no controlled direct effect
			gbound <- c(gbound,1)
		} else {
			gbound <- c(gbound, 1-gbound)
		}
	}
 	g <- suppressWarnings(estimateG(d=data.frame(A,W), g1W, gform, g.SL.library, id=id, verbose, "treatment mechanism", outcome="A")) 
 	g$bound <- gbound
 	if(g$type=="try-error"){
 		stop("Error estimating treatment mechanism (hint: only numeric variables are allowed)") 
 	}
	
  	if(!CDE){
  		g.z <- NULL
  		g.z$type="No intermediate variable"
  		g.z$coef=NA
  		g.Delta <- suppressWarnings(estimateG(d=data.frame(Delta, Z=1, A, W), pDelta1, g.Deltaform, 
 	 		g.SL.library,id=id, verbose, "missingness mechanism", outcome="D")) 
 		g1W.total <- .bound(g$g1W*g.Delta$g1W[,"Z0A1"], gbound)
  		g0W.total <- .bound((1-g$g1W)*g.Delta$g1W[,"Z0A0"], gbound)  
  		if(all(g1W.total==0)){g1W.total <- rep(10^-9, length(g1W.total))}
  		if(all(g0W.total==0)){g0W.total <- rep(10^-9, length(g0W.total))}
  		H1W <- A/g1W.total
  		H0W <- (1-A)/g0W.total

  		suppressWarnings(
  			epsilon <- coef(glm(stage1$Ystar~-1 + offset(Q$Q[,"QAW"]) + H0W + H1W, family=Q$family, subset=Delta==1))
  		)
  		epsilon[is.na(epsilon)] <- 0  # needed for EY1 calculation
 		Qstar <- Q$Q + c((epsilon[1]*H0W + epsilon[2]*H1W), epsilon[1]/g0W.total, epsilon[2]/g1W.total)
		colnames(Qstar) <- c("QAW", "Q0W", "Q1W")
       Ystar <- stage1$Ystar
   		if (maptoYstar) {
			Qstar <- plogis(Qstar)*diff(stage1$ab)+stage1$ab[1] 
			Q$Q <- plogis(Q$Q)*diff(stage1$ab)+stage1$ab[1]
			Ystar <- Ystar*diff(stage1$ab)+stage1$ab[1]
		} else if (family == "poisson"){  	
    		Q$Q <- exp(Q$Q)				  
    		Qstar <- exp(Qstar)
    	}
    	colnames(Q$Q) <- c("QAW", "Q0W", "Q1W")
    	Q$Q <- Q$Q[,-1]
    	res <- calcParameters(Ystar, A, I.Z=rep(1, length(Ystar)), Delta, g1W.total, g0W.total, Qstar, 
   	   		  		mu1=mean(Qstar[,"Q1W"]), mu0=mean(Qstar[,"Q0W"]), id, family)
  		returnVal <- list(estimates=res, Qinit=Q, g=g, g.Z=g.z, g.Delta=g.Delta, Qstar=Qstar[,-1], epsilon=epsilon) 
  		class(returnVal) <- "tmle"
  	} else {
  		returnVal <- vector(mode="list", length=2)
  		g.z <- suppressWarnings(estimateG(d=data.frame(Z,A,W), pZ1, g.Zform, g.SL.library, id=id, 
  					  verbose, "intermediate variable", outcome="Z"))
  		g.Delta <- suppressWarnings(estimateG(d=data.frame(Delta,Z, A, W), pDelta1, g.Deltaform, 
  								 g.SL.library,id=id, verbose, "missingness mechanism", outcome="D")) 
    	ZAD <- cbind(D1Z0A0 = .bound((1-g$g1W)*(1-g.z$g1W[,"A0"])*g.Delta$g1W[,"Z0A0"], gbound),
  					  D1Z0A1 = .bound(g$g1W*(1-g.z$g1W[,"A1"])*g.Delta$g1W[,"Z0A1"], gbound),
  					  D1Z1A0 = .bound((1-g$g1W)*g.z$g1W[,"A0"]*g.Delta$g1W[,"Z1A0"], gbound),
  					  D1Z1A1 = .bound(g$g1W*g.z$g1W[,"A1"]*g.Delta$g1W[,"Z1A1"], gbound))
  	   adjustZero <- colSums(ZAD)==0
  	   ZAD[,adjustZero] <- 10^-9
  	   	for (z in 0:1){
  	     	H0W <- (1-A)*(Z==z)/ZAD[,z*2+1]
  			H1W <- A*(Z==z) /ZAD[,z*2+2]
  			suppressWarnings(
  				epsilon <- coef(glm(stage1$Ystar~-1 + offset(Q$Q[,"QAW"]) + H0W + H1W, family=Q$family,
  					  subset=(Delta==1 & Z==z)))
  			)  			

  			hCounter <- cbind(1/ZAD[,z*2+1], 1/ZAD[,z*2+2])
 			Qstar <- Q$Q[,c(1, z*2+2, z*2+3)] + c((epsilon[1]*H0W + epsilon[2]*H1W), 
 												        epsilon[1]*hCounter[,1], epsilon[2]*hCounter[,2])
			colnames(Qstar) <- c("QAW", "Q0W", "Q1W") 
			newYstar <- stage1$Ystar     
   			if (maptoYstar) {
				Qstar <- plogis(Qstar)*diff(stage1$ab)+stage1$ab[1] 
				Qinit.return <- plogis(Q$Q[,c(1, z*2+2, z*2+3)])*diff(stage1$ab)+stage1$ab[1]
				newYstar <- stage1$Ystar*diff(stage1$ab)+stage1$ab[1]
			} else if (family == "poisson"){ 
    			Qinit.return <- exp(Q$Q[,c(1, z*2+2, z*2+3)]) 
    			Qstar <- exp(Qstar)
    		}
    		colnames(Qinit.return) <- c("QAW", "Q0W", "Q1W")
    		res <- calcParameters(newYstar, A,I.Z=as.integer(Z==z), Delta, g1W=ZAD[,z*2+2], g0W=ZAD[,z*2+1],  
   	   					  Qstar, mu1=mean(Qstar[,"Q1W"]), mu0=mean(Qstar[,"Q0W"]), id, family)
   	   		Qreturn <- Q
   	   		Qreturn$Q <- Qinit.return[,-1]
   	   		returnVal[[z+1]] <- list(estimates=res, Qinit=Qreturn, g=g, g.Z=g.z, g.Delta=g.Delta, 
   	   									 Qstar=Qstar[,-1], epsilon=epsilon)
  		}
  		class(returnVal[[1]]) <- class(returnVal[[2]]) <- "tmle"
  		class(returnVal) <- "tmle.list"
  	}
  	return(returnVal)
}

