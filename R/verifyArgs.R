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
