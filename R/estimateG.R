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

##' Estimate Treatment or Missingness Mechanism
##' 
##' An internal function called by the \code{tmle} function to obtain an
##' estimate of conditional treatment assignment probabiliites \eqn{P(A=1|W)},
##' and conditional probabilites for missingness, \eqn{P(Delta=1|A,W)}.  The
##' estimate can be based on user-supplied values, a user-supplied regression
##' formula, or a data-adaptive super learner fit.  If the \code{SuperLearner}
##' package is not available, and there are no user-specifications, estimation
##' is carried out using main terms regression with \code{glm}.  These main
##' terms-based estimates may yield poor results.
##' 
##' 
##' @param d dataframe with binary dependent variable in the first column,
##' predictors in remaining columns
##' @param g1W vector of values for \eqn{P(A=1|W)}, \eqn{P(Z=1|A,W)}, or
##' \eqn{P(Delta=1|Z,A,W)}
##' @param gform regression formula of the form \code{A~W1}, (dependent
##' variable is one of \eqn{A,Z,D}) if specified this overrides the call to
##' \code{SuperLearner}
##' @param SL.library vector of prediction algorithms used by
##' \code{SuperLearner}, default value is (\sQuote{SL.glm}, \sQuote{SL.step},
##' \sQuote{SL.knn}, \sQuote{SL.DSA.2})
##' @param id subject identifier
##' @param verbose status messages printed if set to TRUE
##' @param message text specifies whether treatment or missingness mechanism is
##' being estimated
##' @param outcome \code{A, D, Z} to indicate which quantity is being
##' estimated.
##' @param newdata optional dataset to be used for prediction after fitting on
##' \code{d}.
##' @return \item{g1W}{a vector containing values for \eqn{P(A=1|W)}, matrix
##' for \eqn{P(Z=1|A,W)}, evaluated at A=0, A=1, or matrix
##' \eqn{P(Delta=1|Z,A,W))} evaluated at (0,0), (0,1), (1,0), (1,1)}
##' \item{coef}{coefficients for each term in the working model used for
##' estimation if \code{glm} was used} \item{type}{estimation procedure}
##' @author Susan Gruber
##' @seealso \code{\link{tmle}}, \code{\link{estimateQ}},
##' \code{\link{calcParameters}}
##' @export
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
