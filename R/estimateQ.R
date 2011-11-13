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