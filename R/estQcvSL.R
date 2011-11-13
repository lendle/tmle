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
