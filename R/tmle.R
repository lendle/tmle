# Targeted Maximum Likelihood Estimation
# for non-parametric estimation of the marginal effect of a binary point
# treatment, adjusting for treatment (g) and missingness (g.Delta) mechanisms
# Parameters include: additive treatment effect: E_W[E(Y|A=1,W) - E(Y|A=0,W)]
# and, for binary outcomes, relative risk(RR) and odds ratio(RR)
#   mu1 = E_W[E(Y|A=1,W)], mu0 = E_W[E(Y|A=0,W)]
#   RR = mu1/mu0
#   OR =  mu1/(1-mu1)/(mu0 / (1-mu0))
# Controlled direct effect estimation available for optional binary Z intermediate variable,
# (P(DeltaY|Z,A,W, Delta=1)P(Delta|Z,A,W)P(Z|A,W)P(A|W)P(W))
# EY1 parameter estimated when there is missinginess and no treatment assignment
# author: Susan Gruber, sgruber@berkeley.edu
# date:   September 25, 2010
# revised: December 20, 2010
# 
# Models or estimated values for Q, g, g.Z, g.Delta can be user-supplied or 
# estimated using super learner.  Cross-validated inital Q can be obtained
# by specifying cvQinit=TRUE

#-------------summary.tmle------------------
# create tmle summary object
# return values of psi, variance, epsilon, 
# and info on estimation of Q and g factors, if available
#---------------------------------------	 


#' Summarization of the results of a call to the tmle routine
#' 
#' These functions are all \link{methods} for class \code{tmle},
#' \code{tmle.list}, \code{summary.tmle}, \code{summary.tmle.list} objects
#' 
#' \code{print.tmle} prints the estimate, variance, p-value, and 95% confidence
#' interval only.  \code{print.summary.tmle}, called indirectly by entering the
#' command \kbd{summary(result)} (where \code{result} has class \code{tmle}),
#' outputs additional information.  Controlled direct effect estimates have
#' class \code{tmle.list}, a list of two objects of class \code{tmle}.  The
#' first item corresponds to \eqn{Z=0}, the second to \eqn{Z=1}
#' 
#' @aliases summary.tmle print.summary.tmle print.tmle summary.tmle.list
#' print.summary.tmle.list print.tmle.list
#' @param object an object of class \code{tmle} or \code{tmle.list}.
#' @param x an object of class \code{tmle} or \code{tmle.list} for summary
#' functions, class \code{summary.tmle} or \code{summary.tmle.list} for print
#' functions.
#' @param \dots currently ignored.
#' @return \item{estimates}{list of parameter estimates, pvalues, and 95%
#' confidence intervals} \item{Qmodel}{working model used to obtain initial
#' estimate of \code{Q} portion of the likelihood, if \code{glm} used}
#' \item{Qterms}{terms in the model for \code{Q}} \item{Qcoef}{coefficient of
#' each term in model for \code{Q}} \item{gmodel}{model used to estimate
#' treatment mechanism \code{g}} \item{gterms}{terms in the treatment mechanism
#' model} \item{gcoef}{coefficient of each term in model for treatment
#' mechanism} \item{gtype}{description of estimation procedure for treatment
#' mechanism, e.g. "glm"} \item{g.Zmodel}{model used to estimate intermediate
#' variable assignment mechanism \code{g.Z}} \item{g.Zterms}{terms in the
#' intermediate mechanism model} \item{g.Zcoef}{coefficient of each term in
#' model for intermediate mechanism} \item{g.Ztype}{description of estimation
#' procedure for intermediate variable} \item{g.Deltamodel}{model used to
#' estimate missingness mechanism \code{g.Delta}} \item{g.Deltaterms}{terms in
#' the missingness mechanism model} \item{g.Deltacoef}{coefficient of each term
#' in model for missingness mechanism} \item{g.Deltatype}{description of
#' estimation procedure for missingness}
#' @author Susan Gruber
#' @seealso \code{\link{tmle}}
#' @examples
#' 
#' # generate data
#'   n <- 500
#'   W <- matrix(rnorm(n*3), ncol=3)
#'   A <- rbinom(n,1, 1/(1+exp(-(.1*W[,1] - .1*W[,2] + .5*W[,3]))))
#'   Y <- A + 2*W[,1] + W[,3] + W[,2]^2 + rnorm(n)
#'   colnames(W) <- paste("W",1:3, sep="")
#' 
#'   result <- tmle(Y,A,W, Qform="Y~A+W1", g1W=rep(.5, n))
#'   summary(result)
#' 
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

#-------------.verifyArgs------------------
# initial checks on data passed in
#-------------------------------------------
.verifyArgs <- function(Y,Z,A,W,Delta, Qform, gform, g.Zform,g.Deltaform){
	formulas <- list(Qform, gform, g.Zform, g.Deltaform)
	validFormula <- sapply(formulas, function(x){identical(class(try(as.formula(x))), "formula")})
	validNames <- c("Z", "A", "Delta", colnames(W))
	validTerms <- rep(TRUE, length(formulas))
	validTerms[validFormula] <- sapply(formulas[which(validFormula)], 
		function(x){is.null(x) || all(attr(terms.formula(as.formula(x), data=data.frame(W)), "term.labels") %in% validNames)})
	ok <- c(length(Y) == length(A) & length(A) == NROW(W) & length(A)==NROW(Delta),
			sum(is.na(A), is.na(W)) == 0,
	 	 	all(A[!is.na(A)] %in% 0:1),
	 	 	is.null(Z) || all(Z[!is.na(Z)] %in% 0:1),
	 	 	is.null(Z) || length(unique(Z)) == 1 | (length(unique(Z)) > 1 & length(unique(A))>1),
			validFormula,
			validTerms
			)
	warning_messages <- c("\t'Y', 'A', 'W', 'Delta', must contain the same number of observations\n",
				"\tNo missing values allowed in 'A' or 'W'\n",
				"\t'A' must be binary (0,1)\n",
				"\t'Z' must be binary (0,1)\n",
				"\tIntermediate variable (Z) not allowed when there is no experimentation in A",
				"\tInvalid regression formula for 'Qform'",
				"\tInvalid regression formula for 'gform'",
				"\tInvalid regression formula for 'g.Zform'",
				"\tInvalid regression formula for 'g.Deltaform'",
				"\tInvalid term name in regression formula for 'Qform'",
				"\tInvalid term name in regression formula for 'gform'",
				"\tInvalid term name in regression formula for 'g.Zform'",
				"\tInvalid term name in regression formula for 'g.Deltaform'"
				)
	if(!all(ok)){
		warning("\n", warning_messages[!ok], immediate. = TRUE)
	}
	return(all(ok))
}

#---------- function .setColnames ---------------
# assign names to every unnamed column of x
# arguments
# 	x.colnames - current column names
#	x.ncols - current number of columns
# 	firstChar - prefix for internally assigned name
# return the names
#-----------------------------------------
.setColnames <- function(x.colnames, x.ncols, firstChar){
	if(is.null(x.colnames)) {
		if(x.ncols > 1){
			x.colnames <- paste(firstChar,1:x.ncols, sep="")
		} else {
			x.colnames <- firstChar
		}
	} else {
		invalid.name <- nchar(x.colnames) == 0
		if(any(invalid.name)){
			x.colnames[invalid.name] <- paste(".internal",firstChar, which(invalid.name), sep="")
		}
	}
	return(x.colnames)
}

#---------- function .bound ---------------
# set outliers to min/max allowable values
# assumes x contains only numerical data
#-----------------------------------------
.bound <- function(x, bounds){
	x[x>max(bounds)] <- max(bounds)
	x[x<min(bounds)] <- min(bounds)
	return(x)
}

#---------- function .initStage1 ---------------
# Bound Y, map to Ystar if applicable, and
# set boundson on Q and enforce on user-specified values
# returns
#   Ystar - outcome values (between [0,1] if maptoYstar=TRUE)
#   Q - matrix of user-specified values
#   Qbounds - bounds on predicted values for Q 
#			(-Inf,+Inf) is default for linear regression
#   ab - bounding levels used to transform Y to Ystar
#-----------------------------------------------
.initStage1 <- function(Y,A, Q, Q.Z1=NULL, Delta, Qbounds, alpha, maptoYstar, family){
	if(family=="binomial") {Qbounds <- c(0,1)}
 	if(is.null(Qbounds)) {
 		if(maptoYstar){ 
 			Qbounds <- range(Y[Delta==1])
 		} else {
 			Qbounds <- c(-Inf, Inf)
 		}
   	}
	if(!is.null(Q)){
   		QAW <- (1-A)*Q[,1] + A*Q[,2]
   		Q <- cbind(QAW, Q0W=Q[,1], Q1W=Q[,2])
   	} 
   	if(!is.null(Q.Z1)){
   		Q <- cbind(Q, Q0W.Z1=Q.Z1[,1], Q1W.Z1=Q.Z1[,2])
   	}
  	ab <- c(0,1)
   	Ystar <- Y
   	if(maptoYstar){ 
   		Ystar <- .bound(Y, Qbounds)
   		if(!is.null(Q)){
   			Q <- .bound(Q, Qbounds)
   		}
   		if(0 >= alpha | 1 <= alpha){
			alpha <- .995
			warning(paste("\n\talpha must be between 0 and 1, alpha reset to",alpha,"\n"),
							immediate. = TRUE)
		}
		ab <- range(Ystar, na.rm=TRUE)
		Ystar[is.na(Ystar)] <- 0  
		Ystar <- (Ystar-ab[1])/diff(ab)	
		if(!is.null(Q)){Q <- (Q-ab[1])/diff(ab)}
		Qbounds <- c(alpha, 1-alpha)
	}	
	return(list(Ystar=Ystar, Q=Q, Qbounds=Qbounds, ab=ab))
} 


#----- function .estQcvSL ----
# purpose: Obtain cross-validated estimates for initial Q using discrete SL.  
# 	This function will return cross-validated predicted initial values 
#   corresponding to the best algorithm in the library, or the convex combination
#   calculated by SL itself, as chosen by cvRSS.
#   The fitted value for each observation i, is predicted from a fit based on a training set
#	excluding the fold containing observation i. Observations with same id are grouped in the
#   same fold.
# arguments: 
# 	Y - outcome
# 	X - design matrix with (Z,A,W)
# 	SL.library - prediction algorithms
# 	V - number of outer cv-folds for disscrete SL
# 	V_SL - number of folds for internal SL cross-validation
# 	family - binomial or gaussian
#   Delta - missingness indicator
#   Qbounds - bounds on predicted values for Q
#   id - subject identifier
# returns:
#  Q - nx5 matrix of predicted values on linear scale (e.g. logit if Qfamily=binomial)
#-----------------------------------------------
.estQcvSL <- function(Y,X,SL.library=NULL, V=5, V_SL=5, family="gaussian", Delta, Qbounds, id, verbose){
	SL.version <- 2
	Q <- cvRSS <- best_alg <- NULL
	n <- length(Y)
	u.id <- unique(id)
	n.id <- length(u.id)
	fold <- by(sample(1:n.id),rep(1:V, length.out=n.id),function(x){which(id %in% u.id[x])})
	n_predictors <-length(SL.library)
	CDE <- length(unique(X[,1])) > 1

	if(NCOL(X) > 1 & require(SuperLearner)){
	  if(packageDescription("SuperLearner")$Version < SL.version){
		newX <- rbind(X,X,X,X,X)
		newX[(n+1):(3*n),1]   <- 0 
		newX[(3*n+1):(5*n),1] <- 1
		newX[(n+1):(2*n),2]   <- 0
		newX[(2*n+1):(3*n),2] <- 1
		newX[(3*n+1):(4*n),2] <- 0
		newX[(4*n+1):(5*n),2] <- 1
		# We'll create a matrix of predictions - one column for each predictor in the library
		# plus one more for SL itself, with 5*n predicted values per column, corresponding to newX.
   		predictions <- matrix(nrow=5*n, ncol=length(SL.library)+1)
   		m_SL <- NULL
    	for (v in 1:V) {
    		fold_rows <- c(fold[[v]], fold[[v]]+n, fold[[v]]+2*n, fold[[v]]+3*n, fold[[v]]+4*n)
    		if(class(m_SL) != "try-error"){
    			train.observed <- (1:n)[-fold[[v]]][Delta[-fold[[v]]]==1]
    			suppressWarnings(
    				m_SL <- try(SuperLearner(Y=Y[train.observed], X=X[train.observed,], newX=newX[fold_rows,],
    			 	V=V_SL, save.fit.library=FALSE, family=family,SL.library=SL.library,id=id[train.observed]))
    			 )
    		}
    		if(class(m_SL) != "try-error"){
    			predictions[fold_rows,1] <- m_SL$SL.predict
    			for (s in 1:n_predictors){
    				predictions[fold_rows,s+1] <- m_SL$library.predict[,s]
    			}
    			predictions <- .bound(predictions, Qbounds)
    		}
    	} 
    	cvRSS <- colSums(Delta*(Y-predictions[1:n,])^2)
    	names(cvRSS) <- c("SL", SL.library)
    	best <- which.min(cvRSS)
    	best_alg <- c("SL", SL.library)[best]
    	Q <- matrix(data=predictions[,best], nrow=n, ncol=5, byrow=FALSE)
    	colnames(Q) <- c("QAW", "Q0W", "Q1W", "Q0W.Z1", "Q1W.Z1")
    } 	}
	if(verbose){cat("\tDiscrete SL: best algorithm = ", best_alg,"\n")}
	if (is.null(Q) | class(m_SL) == "try-error"){
		Q <- 0
		class(Q) <- "try-error"
	}
	Qinit <- list(Q=Q, family=family, SL.library=SL.library, cvRSS=cvRSS, best_alg=best_alg)
    return(Qinit)
 }

#-----------estimateQ----------------
# purpose: estimate Q=E(Y |Z, A,W) data-adaptively,
# unless super learner not available, or user specifies 
# initial values or a regression formula
# arguments: 
# 	Y - outcome 
# 	Z - intermediate variable between A and Y (default= 0 when no int. var.) 
#	A - treatment indicator (1=treatment, 0=control)
# 	W - baseline covariates
#	Delta - missingness indicator
#	Q - optional externally estimated values for Q
#	Qbounds - bounds for predicted values 
#  	Qform - optional regression formula to use for glm if 
#	        non-data adaptive estimation specified
# 	maptoYstar - if TRUE, using logistic fluctuation for bounded, continuous outcomes
# 		estimation inital Q on linear scale, bounded by (0,1),and return on logit scale
#		(will work if family=poisson)
#	SL.library - library of prediction algorithms for Super Learner
#   cvQinit - flag, if TRUE, cross-validate SL.
# 	family - regression family
#	id - subject identifier
# returns matrix of linear predictors for Q(A,W), Q(0,W), Q(1,W),
#   (for controlled direct effect, 2 additional columns: Q(Z=1,A=0,W), Q(Z=1,A=1,W)) 
#		family for stage 2 targeting
#		coef, NA, unless Q is estimated using a parametric model
# 		type, estimation method for Q
#----------------------------------------


#' Initial Estimation of Q portion of the Likelihood
#' 
#' An internal function called by the \code{tmle} function to obtain an initial
#' estimate of the \eqn{Q} portion of the likelihood based on user-supplied
#' matrix values for predicted values of (counterfactual outcomes)
#' \code{Q(0,W),Q(1,W)}, or a user-supplied regression formula, or based on a
#' data-adaptively selected \code{SuperLearner} fit.  If \code{SuperLearner} is
#' not available and values for \code{Q} are not user-supplied, estimation is
#' based on a main terms regression using \code{glm}.
#' 
#' 
#' @param Y continuous or binary outcome variable
#' @param Z optional binary indicator for intermediate covariate for conrolled
#' direct effect estimation
#' @param A binary treatment indicator, \code{1} - treatment, \code{0} -
#' control
#' @param W vector, matrix, or dataframe containing baseline covariates
#' @param Delta indicator of missing outcome. \code{1} - observed, \code{0} -
#' missing
#' @param Q 3-column matrix \code{(Q(A,W), Q(0,W), Q(1,W))}
#' @param Qbounds Bounds on predicted values for \code{Q}, set to \code{alpha}
#' for logistic fluctuation, or \code{range(Y)} if not user-supplied
#' @param Qform regression formula of the form \code{Y~A+W}
#' @param maptoYstar if \code{TRUE} indicates continuous \code{Y} values should
#' be shifted and scaled to fall between (0,1)
#' @param SL.library specification of prediction algorithms, default is
#' (\sQuote{SL.glm}, \sQuote{SL.step}, \sQuote{SL.DSA.2}) when
#' family=\sQuote{gaussian}, (\sQuote{SL.glm}, \sQuote{SL.step},
#' \sQuote{SL.knn}, \sQuote{SL.DSA.2}) when family = \sQuote{binomial}.  In
#' practice, including more prediction algorithms in the library improves
#' results.
#' @param cvQinit logical, whether or not to estimate cross-validated values
#' for initial \code{Q}, default=\code{FALSE}
#' @param family family specification for regressions, generally
#' \sQuote{gaussian} for continuous oucomes, \sQuote{binomial} for binary
#' outcomes
#' @param id subject identifier
#' @param verbose status message printed if set to \code{TRUE}
#' @return \item{Q}{\eqn{nx3} matrix, columns contain the initial estimate of
#' \eqn{[Q(A,W)=E(Y|A=a,W), Q(0,W)=E(Y|A=0,W), Q(1,W)=E(Y|A=1,W)]}. For
#' controlled direct estimation, \eqn{nx5} matrix, \eqn{E(Y|Z,A,W)}, evaluated
#' at \eqn{(z,a), (0,0), (0,1), (1,0), (1,1)} on scale of linear predictors}
#' \item{Qfamily}{\sQuote{binomial} for targeting with logistic fluctuation,
#' \sQuote{gaussian} for linear fluctuation} \item{coef}{coefficients for each
#' term in working model used for initial estimation of \code{Q} if \code{glm}
#' used.} \item{type}{type of estimation procedure}
#' @author Susan Gruber
#' @seealso \code{\link{tmle}}, \code{\link{estimateG}},
#' \code{\link{calcParameters}}
estimateQ <- function (Y,Z,A,W, Delta, Q, Qbounds, Qform, maptoYstar, 
		SL.library, cvQinit, family, id, verbose) {
	SL.version <- 2
	Qfamily <- family
	m <- NULL
	coef <- NA
	CDE <- length(unique(Z)) > 1
	type <- "user-supplied values"
	if(is.null(Q)){
		if(verbose) { cat("\tEstimating initial regression of Y on A and W\n")}
		Q <- matrix(NA, nrow=length(Y), ncol = 5)
  	  	colnames(Q)<- c("QAW", "Q0W", "Q1W", "Q0W.Z1", "Q1W.Z1")
  	  	if(!(is.null(Qform))){
  	  		if(identical(as.character(as.formula(Qform)), c("~","Y", "."))){
  	  			if(CDE){
  	  				Qform <- paste("Y~Z+A+", paste(colnames(W), collapse="+"))
  	  			} else {
  	   	  			Qform <- paste("Y~A+", paste(colnames(W), collapse="+"))
				}
  	  		}
  	  		m <- suppressWarnings(glm(Qform, data=data.frame(Y,Z,A,W, Delta), family=family, subset=Delta==1))
  	  		Q[,"QAW"] <- predict(m, newdata=data.frame(Y,Z,A,W), type="response")
  	  		Q[,"Q0W"] <- predict(m, newdata=data.frame(Y,Z=0,A=0,W), type="response")
  	  		Q[,"Q1W"] <- predict(m, newdata=data.frame(Y,Z=0,A=1,W), type="response")
  	  		Q[,"Q0W.Z1"] <- predict(m, newdata=data.frame(Y,Z=1,A=0,W), type="response")
  	  		Q[,"Q1W.Z1"] <- predict(m, newdata=data.frame(Y,Z=1,A=1,W), type="response")
  	  		coef <- coef(m)
  	  		type="glm, user-supplied model"	
  	  	} else {
  	  		if(cvQinit){
  	  			m <- try(.estQcvSL(Y,X=cbind(Z,A,W),SL.library, family=family, 
  	  					Delta=Delta, Qbounds=Qbounds,id=id, verbose=verbose))
  	  			if(!(identical(class(m), "try-error"))){
  	  				type <- "cross-validated SL"
  	  				Qinit <- m
  	  				Q <- Qinit$Q
  	  			}
  	  		} else {
  	  			if(require(SuperLearner)){
  	  			  if(packageDescription("SuperLearner")$Version < SL.version){
  					if(verbose) {cat("\t using SuperLearner\n")}
  					n <- length(Y)
  					X <- data.frame(Z,A,W)
  					X00 <- data.frame(Z=0,A=0, W)
  					X01 <- data.frame(Z=0,A=1, W)
  					newX <- rbind(X, X00, X01)
  					if(CDE) {
  						X10 <- data.frame(Z=1,A=0, W)
  						X11 <- data.frame(Z=1,A=1, W)
  						newX <- rbind(newX, X10, X11)
  					}
  					suppressWarnings(
  						m <- try(SuperLearner(Y[Delta==1],X[Delta==1,], newX=newX, SL.library=SL.library,
  							V=5, family=family, save.fit.library=FALSE, id=id[Delta==1]))
  					)
  					if(identical(class(m),"SuperLearner")){
  						Q[,"QAW"] <- m$SL.predict[1:n]
  						Q[,"Q0W"] <- m$SL.predict[(n+1):(2*n)]
  						Q[,"Q1W"] <- m$SL.predict[(2*n+1):(3*n)]
  						if(CDE){
  							Q[,"Q0W.Z1"] <- m$SL.predict[(3*n+1):(4*n)]
  				  			Q[,"Q1W.Z1"] <- m$SL.predict[(4*n+1):(5*n)]
  						}
  						type <- "SuperLearner"
  	  	} } } }}
  	} 
  	if(is.na(Q[1,1]) | identical(class(m), "try-error")){
  	  		if(verbose) {cat("\t Running main terms regression for 'Q' using glm\n")}
  	  		Qform <- paste("Y~Z+A+", paste(colnames(W), collapse="+"))
  	  		m <- glm(Qform, data=data.frame(Y,Z,A,W, Delta), family=family, subset=Delta==1)
  	  		Q[,"QAW"] <- predict(m, newdata=data.frame(Y,Z,A,W), type="response")
  	  		Q[,"Q1W"] <- predict(m, newdata=data.frame(Y,Z=0,A=1,W), type="response")
  	  		Q[,"Q0W"] <- predict(m, newdata=data.frame(Y,Z=0,A=0,W), type="response")
	  		Q[,"Q0W.Z1"] <- predict(m, newdata=data.frame(Y,Z=1,A=0,W), type="response")
  	  		Q[,"Q1W.Z1"] <- predict(m, newdata=data.frame(Y,Z=1,A=1,W), type="response")

  	  		coef <- coef(m)
  	  		type="glm, main terms model"
  	 }
	Q <- .bound(Q, Qbounds)
	if(maptoYstar | identical(Qfamily,"binomial") | identical(Qfamily, binomial)){
			Q <- qlogis(Q)
			Qfamily <- "binomial"	
	} else if (identical(Qfamily, "poisson") | identical(Qfamily, poisson)) {
			Q <- log(Q)
			Qfamily <- "poisson"
	}
	if(!CDE){
		Q <- Q[,1:3]
	}
	if(cvQinit){
		Qinit$Q <- Q
	} else {
		Qinit <- list(Q=Q, family=Qfamily, coef=coef, type=type)
		if(type=="SuperLearner"){
			 Qinit$SL.library=SL.library
			 Qinit$coef=m$coef
		}
	}
	return(Qinit)
}

#-----------estimateG----------------
# Estimate factors of g
# 		P(A=1|W), P(Z=1|A,W), P(Delta=1|Z,A,W)
# d - dataframe (A,W), (Z,A,W), or (Delta,Z,A,W)
# g1W - optional vector/matrix  of externally estimated values
# gform - optionalformula to use for glm
# SL.library - algorithms to use for super learner estimation
#   id - subject identifier
# verbose - flag, whether or not to print messages
# message - printed when verbose=TRUE 
# outcome - "A" for treatment, "Z" for intermediate variable,
#           "D" for Delta (missingness)
# newdata - optional values to predict on (needed by tmleMSM function)
# d = [A,W] for treatment
# d = [Z,A,W] for intermediate
# d = [Delta, Z,A,W for missingness]
#----------------------------------------


#' Estimate Treatment or Missingness Mechanism
#' 
#' An internal function called by the \code{tmle} function to obtain an
#' estimate of conditional treatment assignment probabiliites \eqn{P(A=1|W)},
#' and conditional probabilites for missingness, \eqn{P(Delta=1|A,W)}.  The
#' estimate can be based on user-supplied values, a user-supplied regression
#' formula, or a data-adaptive super learner fit.  If the \code{SuperLearner}
#' package is not available, and there are no user-specifications, estimation
#' is carried out using main terms regression with \code{glm}.  These main
#' terms-based estimates may yield poor results.
#' 
#' 
#' @param d dataframe with binary dependent variable in the first column,
#' predictors in remaining columns
#' @param g1W vector of values for \eqn{P(A=1|W)}, \eqn{P(Z=1|A,W)}, or
#' \eqn{P(Delta=1|Z,A,W)}
#' @param gform regression formula of the form \code{A~W1}, (dependent variable
#' is one of \eqn{A,Z,D}) if specified this overrides the call to
#' \code{SuperLearner}
#' @param SL.library vector of prediction algorithms used by
#' \code{SuperLearner}, default value is (\sQuote{SL.glm}, \sQuote{SL.step},
#' \sQuote{SL.knn}, \sQuote{SL.DSA.2})
#' @param id subject identifier
#' @param verbose status messages printed if set to TRUE
#' @param message text specifies whether treatment or missingness mechanism is
#' being estimated
#' @param outcome \code{A, D, Z} to indicate which quantity is being estimated.
#' @param newdata optional dataset to be used for prediction after fitting on
#' \code{d}.
#' @return \item{g1W}{a vector containing values for \eqn{P(A=1|W)}, matrix for
#' \eqn{P(Z=1|A,W)}, evaluated at A=0, A=1, or matrix \eqn{P(Delta=1|Z,A,W))}
#' evaluated at (0,0), (0,1), (1,0), (1,1)} \item{coef}{coefficients for each
#' term in the working model used for estimation if \code{glm} was used}
#' \item{type}{estimation procedure}
#' @author Susan Gruber
#' @seealso \code{\link{tmle}}, \code{\link{estimateQ}},
#' \code{\link{calcParameters}}
estimateG <- function (d,g1W, gform,SL.library, id, verbose, message, outcome="A", newdata=d)  {
  SL.version <- 2
  SL.ok <- FALSE
  m <- NULL
  coef <- NA
  type <- NULL
  if (is.null(g1W)){
  	if(verbose){cat("\tEstimating", message, "\n")}
	if (length(unique(d[,1]))==1) {
		g1W <- rep(1,nrow(d))
		type <- paste("No", strsplit(message, " ")[[1]][1])
		if(outcome=="Z"){
			g1W <- cbind(A0=g1W, A1=g1W)
		} else if (outcome=="D"){
			g1W <- cbind(Z0A0=g1W, Z0A1=g1W, Z1A0=g1W, Z1A1=g1W)
		}
	} else {
  	  if (is.null(gform)){
  		if(require(SuperLearner)){
  		 SL.ok <- packageDescription("SuperLearner")$Version < SL.version
  		 if(SL.ok){
  			suppressWarnings(
  				m <- try(SuperLearner(Y=d[,1], X=d[,-1], newX=newdata[,-1], family="binomial", SL.library=SL.library,
  							V=5, id=id))
  			)
  			if(identical(class(m),"SuperLearner")) {
  				g1W <- as.vector(m$SL.predict)
  			}
  		 }
  	    } 
  	    if (!SL.ok){
  			if(verbose){cat("\tRunning main terms regression for 'g' using glm\n")}
			form <- paste(paste(colnames(d)[1],"~1"), paste(colnames(d)[-1], collapse = "+"), sep="+")  
			m <- glm(form, data=d, family="binomial")
			g1W <- predict(m, newdata=newdata, type="response")
			coef <- coef(m)
  		}
  	  } else {
  	  	form <- try(as.formula(gform))
  	  	if(class(form)== "formula") {
  	  		m <- try(glm(form,  data=d, family="binomial"))
  	  		if (class(m)[1]=="try-error"){
  				if(verbose){cat("\tInvalid formula supplied. Running glm using main terms\n")}
				form <- paste(colnames(d)[1],"~1", paste(colnames(d)[-1], collapse = "+"), sep="+")  
  	  			m <- glm(form, data=d, family="binomial")
  	  		} else {
  	  			type <- "user-supplied regression formula"
  	  		}
  	  	} else {
  	  	  	if(verbose){cat("\tRunning main terms regression for 'g' using glm\n")}
			form <- paste(colnames(d)[1],"~1", paste(colnames(d)[-1], collapse = "+"), sep="+")  
  	  		m <- glm(form, data=d, family="binomial")
  	  	}
  	  	 g1W <- predict(m, newdata=newdata, type="response")
  	  	 coef <- coef(m)
  	  }
  	  # Get counterfactual predicted values
  	  if(outcome=="Z"){
  	  	if(identical(class(m),"SuperLearner")){
  	  		g1W <- cbind(predict(m, newdata=data.frame(A=0, newdata[,-(1:2), drop=FALSE]), type="response", 
  	  								X=d[,-1, drop=FALSE], Y=d[,1])$fit, 
  	  				      predict(m, newdata=data.frame(A=1, newdata[,-(1:2), drop=FALSE]), type="response", 
  	  				   				X=d[,-1, drop=FALSE], Y=newdata[,1])$fit)
  	  	} else {
  	  		g1W <- cbind(predict(m, newdata=data.frame(A=0, newdata[,-(1:2), drop=FALSE]), type="response"), 
  	  				   predict(m, newdata=data.frame(A=1, newdata[,-(1:2), drop=FALSE]), type="response"))
  	  	}	
  	  	colnames(g1W) <- c("A0", "A1")
		
  	  } else if (outcome=="D"){
  	  	if(identical(class(m),"SuperLearner")){
  	  		g1W <- cbind(predict(m, newdata=data.frame(Z=0, A=0, newdata[,-(1:3), drop=FALSE]), type="response", 
  	  							X=d[,-1,drop=FALSE], Y=d[,1])$fit, 
  	  			 		  predict(m, newdata=data.frame(Z=0, A=1, newdata[,-(1:3), drop=FALSE]), type="response", 
  	  			 				X=d[,-1, drop=FALSE], Y=d[,1])$fit,
  	  		  	 		  predict(m, newdata=data.frame(Z=1, A=0, newdata[,-(1:3), drop=FALSE]), type="response", 
  	  		  	 				X=d[,-1, drop=FALSE], Y=d[,1])$fit,	     	     		  predict(m, newdata=data.frame(Z=1, A=1, newdata[,-(1:3), drop=FALSE]), type="response", 
  	  		  	 				X=d[,-1, drop=FALSE], Y=d[,1])$fit)
  	   } else{
  	   	 	g1W <- cbind(predict(m, newdata=data.frame(Z=0, A=0, newdata[,-(1:3), drop=FALSE]), type="response"), 
  	  			 		  predict(m, newdata=data.frame(Z=0, A=1, newdata[,-(1:3), drop=FALSE]), type="response"),
  	  		  	 		  predict(m, newdata=data.frame(Z=1, A=0, newdata[,-(1:3), drop=FALSE]), type="response"),	     	     		  predict(m, newdata=data.frame(Z=1, A=1, newdata[,-(1:3), drop=FALSE]), type="response"))
  	   }
  	   colnames(g1W) <- c("Z0A0", "Z0A1", "Z1A0", "Z1A1")
  	  }
  	}	
} else {
  		type <- "user-supplied values"
  		if(outcome=="Z") {
  			colnames(g1W) <- c("A0", "A1")
  		} else if (outcome=="D"){
  			colnames(g1W) <- c("Z0A0", "Z0A1", "Z1A0", "Z1A1")[1:ncol(g1W)]
  		}
  	}
  	if(is.null(type)){ type <- class(m)[1]}
  	returnVal <- list(g1W=g1W, coef=coef, type=type)
  	if(type=="SuperLearner"){
			 returnVal$SL.library=SL.library
			 returnVal$coef=m$coef
		}
 return(returnVal)
}


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


#' Calculate Parameter Estimates (calcParameters)
#' 
#' An internal function called by the \code{tmle} function to calculate the
#' population mean effect when there is missingness in the data, but no
#' treatment assignment.  When observations are in treatment and control
#' groups, estimates the additive treatment effect, and if the outcome is
#' binary, also the relative risk and odds ratio parameters.  P-values and 95%
#' confidence intervals are also calculated (on the log scale for RR and OR).
#' 
#' 
#' @param Y continuous or binary outcome variable
#' @param A binary treatment indicator, \code{1} - treatment, \code{0} -
#' control
#' @param I.Z Indicator Z=z, needed for CDE estimation
#' @param Delta indicator of missing outcome. \code{1} - observed, \code{0} -
#' missing
#' @param g1W censoring mechanism estimates, \eqn{P(A=1|W)*P(Delta=1|A,W)}
#' @param g0W censoring mechanism estimates, \eqn{P(A=0|W)*P(Delta=1|A,W)}
#' @param Q a 3-column matrix \code{(Q(A,W), Q(1,W), Q(0,W))}
#' @param mu1 targeted estimate of \eqn{E(Y|A=1,W)}
#' @param mu0 targeted estimate of \eqn{E(Y|A=0,W)}
#' @param id subject identifier
#' @param family family specification for regressions, generally
#' \sQuote{gaussian} for continuous oucomes, \sQuote{binomial} for binary
#' outcomes
#' @return \item{EY1}{Population mean outcome estimate, variance, p-value, 95%
#' confidence interval (missingness only, no treatment assingment), or
#' \code{NULL}} \item{ATE}{additive treatment effect estimate, variance,
#' p-value, 95% confidence interval, or \code{NULL}} \item{RR}{relative risk
#' estimate, p-value, 95% confidence interval, log(RR), variance(log(RR)), or
#' \code{NULL}} \item{OR}{odds ratio estimate, p-value, 95% confidence
#' interval, log(OR), variance(log(OR)), or \code{NULL}}
#' @author Susan Gruber
#' @seealso \code{\link{tmle}}, \code{\link{estimateQ}},
#' \code{\link{estimateG}}
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


#' Targeted Maximum Likelihood Estimation
#' 
#' Carries out targeted maximum likelihood estimation of marginal additive
#' treatment effect of a binary point treatment on an outcome.  For binary
#' outcomes risk ratio and odds ratio estimates are also reported. The
#' \code{tmle} function is minimally called with arguments \code{(Y,A,W)},
#' where \code{Y} is a continuous or binary outcome variable, \code{A} is a
#' binary treatment variable, (\code{A=1} for treatment, \code{A=0} for
#' control), and \code{W} is a matrix or dataframe of baseline covariates. The
#' population mean outcome is calculated when there is no variation in
#' \code{A}. If values of binary mediating variable \code{Z} are supplied,
#' estimates are returned at each level of \code{Z}. Missingness in the outcome
#' is accounted for in the estimation procedure if missingness indicator
#' \code{Delta} is 0 for some observations.  Repeated measures can be
#' identified using the \code{id} argument.
#' 
#' \code{gbounds} defaults to (0.025, 0.975).  To specify symmetrical
#' truncation levels, only one value has to be provided. The upper bound is set
#' to 1 for EY1 parameter estimated when there is missingness only, no
#' treatment assignment.
#' 
#' Controlled direct effects are estimated when binary covariate \code{Z} is
#' non-null. The tmle function returns an object of class \code{tmle.list}, a
#' list of two items of class \code{tmle}.  The first corresponds to estimates
#' obtained when \code{Z} is fixed at \eqn{0}, the second correspondes to
#' estimates obtained when \code{Z} is fixed at \eqn{1}.
#' 
#' \code{Q.SL.library} defaults to (\sQuote{SL.glm}, \sQuote{SL.step},
#' \sQuote{SL.DSA.2}) when family=\sQuote{gaussian}, (\sQuote{SL.glm},
#' \sQuote{SL.step}, \sQuote{SL.knn}, \sQuote{SL.DSA.2}) when family =
#' \sQuote{binomial}
#' 
#' \code{g.SL.library} Defaults to (\sQuote{SL.glm}, \sQuote{SL.step},
#' \sQuote{SL.knn}, \sQuote{SL.DSA.2})
#' 
#' See \code{SuperLearner} help files for further information.
#' 
#' @param Y continuous or binary outcome variable
#' @param A binary treatment indicator, \code{1} - treatment, \code{0} -
#' control
#' @param W vector, matrix, or dataframe containing baseline covariates
#' @param Z optional binary indicator for intermediate covariate for conrolled
#' direct effect estimation
#' @param Delta indicator of missing outcome or treatment assignment.  \code{1}
#' - observed, \code{0} - missing
#' @param Q \eqn{E(Y|A,W)}, optional \eqn{nx2} matrix of initial values for
#' \eqn{Q} portion of the likelihood, \eqn{(E(Y|A=0,W), E(Y|A=1,W))}
#' @param Q.Z1 \eqn{E(Y|Z=1,A,W)}, optional \eqn{nx2} matrix of initial values
#' for \eqn{Q} portion of the likelihood, \eqn{(E(Y|Z=1,A=0,W),
#' E(Y|Z=1,A=1,W))}. (When specified, values for \eqn{E(Y|Z=0,A=0,W),
#' E(Y|Z=0,A=1,W)} are passed in using the \code{Q} argument
#' @param Qform optional regression formula for estimation of \eqn{E(Y|A, W)},
#' suitable for call to \code{glm}
#' @param Qbounds vector of upper and lower bounds on \code{Y} and predicted
#' values for initial \code{Q}
#' @param Q.SL.library optional vector of prediction algorithms to use for
#' \code{SuperLearner} estimation of initial \code{Q}
#' @param cvQinit logical, if \code{TRUE}, estimates cross-validated predicted
#' values using discrete super learning, default=\code{FALSE}
#' @param g1W optional vector of conditional treatment assingment
#' probabilities, \eqn{P(A=1|W)}
#' @param gform optional regression formula of the form \code{A~W}, if
#' specified this overrides the call to \code{SuperLearner}
#' @param gbound value between (0,1) for truncation of predicted probabilities.
#' See \code{Details} section for more information
#' @param pZ1 optional\eqn{nx2} matrix of conditional probabilities
#' \eqn{P(Z=1|A=0,W), P(Z=1|A=1,W)}
#' @param g.Zform optionalregression formula of the form \code{Z~A+W}, if
#' specified this overrides the call to \code{SuperLearner}
#' @param pDelta1 optional matrix of conditional probabilities for missingness
#' mechanism, \eqn{nx2} when \code{Z} is \code{NULL} \eqn{P(Delta=1|A=0,W),
#' P(Delta=1|A=1,W)}.  \eqn{nx4} otherwise, \eqn{P(Delta=1|Z=0,A=0,W),
#' P(Delta=1|Z=0,A=1,W), P(Delta=1|Z=1,A=0,W), P(Delta=1|Z=1,A=1,W)}
#' @param g.Deltaform optional regression formula of the form \code{Delta~A+W},
#' if specified this overrides the call to \code{SuperLearner}
#' @param g.SL.library optional vector of prediction algorithms to use for
#' \code{SuperLearner} estimation of \code{g1W} or \code{pDelta1}
#' @param family family specification for working regression models, generally
#' \sQuote{gaussian} for continuous outcomes (default), \sQuote{binomial} for
#' binary outcomes
#' @param fluctuation \sQuote{logistic} (default), or \sQuote{linear}
#' @param alpha used to keep predicted initial values bounded away from (0,1)
#' for logistic fluctuation
#' @param id optional subject identifier
#' @param verbose status messages printed if set to \code{TRUE}
#' (default=\code{FALSE})
#' @return \item{estimates}{list with elements EY1 (population mean), ATE
#' (additive treatment effect), RR (relative risk), OR (odds ratio). Each
#' element in the estimates of these is itself a list containing \itemize{
#' \item psi - parameter estimate \item pvalue - two-sided p-value \item CI -
#' 95% confidence interval \item var.psi - Influence-curve based variance of
#' estimate (ATE parameter only) \item log.psi - Parameter estimate on log
#' scale (RR and OR parameters) \item var.log.psi - Influence-curve based
#' variance of estimate on log scale (RR and OR parameters) }}
#' \item{Qinit}{initial estimate of \code{Q}. \code{Qinit$coef} are the
#' coefficients for a \code{glm} model for \code{Q}, if applicable.
#' \code{Qinit$Q} is an \eqn{nx2} matrix, where \code{n} is the number of
#' observations.  Columns contain predicted values for \code{Q(0,W),Q(1,W)}
#' using the initial fit.  \code{Qinit$type} is method for estimating \code{Q}}
#' \item{Qstar}{targeted estimate of \code{Q}, an \eqn{nx2} matrix with
#' predicted values for \code{Q(0,W), Q(1,W)} using the updated fit}
#' \item{g}{treatment mechanism estimate. A list with three items: \code{g$g1W}
#' contains estimates of \eqn{P(A=1|W)} for each observation, \code{g$coef} the
#' coefficients for the model for \eqn{g} when \code{glm} used, \code{g$type}
#' estimation procedure} \item{g.Z}{intermediate covariate assignment estimate
#' (when applicable). A list with three items: \code{g.Z$g1W} an \eqn{nx2}
#' matrix containing values of \eqn{P(Z=1|A=1,W), P(Z=1|A=0,W)} for each
#' observation, \code{g.Z$coef} the coefficients for the model for \eqn{g} when
#' \code{glm} used, \code{g.Z$type} estimation procedure}
#' \item{g.Delta}{missingness mechanism estimate. A list with three items:
#' \code{g.Delta$g1W} an \eqn{nx4} matrix containing values of
#' \eqn{P(Delta=1|Z,A,W)} for each observation, with (Z=0,A=0), (Z=0,A=1),
#' (Z=1,A=0),(Z=1,A=1). (When \code{Z} is \code{NULL}, columns 3 and 4 are
#' duplicates of 1 and 2.) \code{g.Delta$coef} the coefficients for the model
#' for \eqn{g} when \code{glm} used, \code{g.Delta$type} estimation procedure}
#' @author Susan Gruber \email{sgruber@@berkeley.edu}, in collaboration with
#' Mark van der Laan.
#' @seealso \code{\link{summary.tmle}}, \code{\link{estimateQ}},
#' \code{\link{estimateG}}, \code{\link{calcParameters}}
#' @references 1. Gruber, S. and van der Laan, M.J. (2009), Targeted Maximum
#' Likelihood Estimation: A Gentle Introduction. \emph{U.C. Berkeley Division
#' of Biostatistics Working Paper Series}.  Working Paper 252.
#' \url{http://www.bepress.com/ucbbiostat/paper252}
#' 
#' 2. Gruber, S. and van der Laan, M.J.  (2010), A Targeted Maximum Likelihood
#' Estimator of a Causal Effect on a Bounded Continuous Outcome. \emph{The
#' International Journal of Biostatistics}, 6(1), 2010.
#' 
#' 3. van der Laan, M.J. and Rubin, D. (2006), Targeted Maximum Likelihood
#' Learning. \emph{The International Journal of Biostatistics}, 2(1).
#' \url{http://www.bepress.com/ijb/vol2/iss1/11/}
#' 
#' 4. van der Laan, M.J., Rose, S., and Gruber,S., editors, (2009) Readings in
#' Targeted Maximum Likelihood Estimation . \emph{U.C. Berkeley Division of
#' Biostatistics Working Paper Series}.  Working Paper 254.
#' \url{http://www.bepress.com/ucbbiostat/paper254}
#' @examples
#' 
#' library(tmle)
#' 
#' # generate data
#'   set.seed(5)
#'   n <- 500
#'   W <- matrix(rnorm(n*3), ncol=3)
#'   A <- rbinom(n,1, 1/(1+exp(-(.1*W[,1] - .1*W[,2] + .5*W[,3]))))
#'   Y <- A + 2*W[,1] + W[,3] + W[,2]^2 + rnorm(n)
#'   colnames(W) <- paste("W",1:3, sep="")
#' 
#' # Example 1. Simplest function invocation 
#' # If available, SuperLearner called to estimate Q, g
#' # otherwise, main terms regression using glm
#' # Delta defaults to 1 for all observations   
#'   result1 <- tmle(Y,A,W)
#'   summary(result1)
#' 
#' # Example 2: 
#' # User-supplied regression formulas to estimate Q and g
#' # binary outcome
#'   n <- 250
#'   W <- matrix(rnorm(n*3), ncol=3)
#'   colnames(W) <- paste("W",1:3, sep="")
#'   A <- rbinom(n,1, plogis(0.6*W[,1] +0.4*W[,2] + 0.5*W[,3]))
#'   Y <- rbinom(n,1, plogis(A + 0.2*W[,1] + 0.1*W[,2] + 0.2*W[,3]^2 ))
#'   result2 <- tmle(Y,A,W, family="binomial", Qform=Y~A+W1+W2+W3, gform=A~W1+W2+W3)
#'   summary(result2)
#' 
#' # Example 3: Population mean outcome
#' # User-supplied (misspecified) model for Q, 
#' # If available, Super learner called to estimate g, g.Delta, otherwise glm
#' # approx. 20% missing at random
#'   Delta <- rbinom(n, 1, 1/(1+exp(-(1.7-1*W[,1]))))
#'   result3 <- tmle(Y,A=NULL,W, Delta=Delta, Qform="Y~A+W1+W2+W3")
#'   summary(result3)
#' 
#' # Example 4: Controlled direct effect with missingness
#' # User-supplied models for Q, g, g.Z
#' # User-supplied values for pDelta1 
#'   n <- 1000
#'   W <- matrix(rnorm(n*3), ncol = 3)
#'   colnames(W) <- paste("W", 1:3, sep = "")
#'   A <- rbinom(n,1, plogis(0.6*W[,1] + 0.4*W[,2] + 0.5*W[,3]))
#'   Z <- rbinom(n,1, plogis(0.5 + A))
#'   Y <- A + A*Z+ 0.2*W[,1] + 0.1*W[,2] + 0.2*W[,3]^2 + rnorm(n)
#'   Delta <- rbinom(n,1, plogis(Z + A)) 
#'   pDelta1 <- cbind(rep(plogis(0), n), rep(plogis(1), n),
#'                    rep(plogis(1), n), rep(plogis(2), n))
#'   colnames(pDelta1) <- c("Z0A0", "Z0A1", "Z1A0", "Z1A1")
#'   Y[Delta == 0] <- NA
#'   result4 <- tmle(Y, A, W, Z, Delta = Delta, pDelta1= pDelta1, 
#'                   Qform = Y ~ 1, g.SL.library = "SL.glm") 
#'   result4
#' 
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

