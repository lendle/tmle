#-------------summary.tmle------------------
# create tmle summary object
# return values of psi, variance, epsilon, 
# and info on estimation of Q and g factors, if available
#---------------------------------------	 
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

