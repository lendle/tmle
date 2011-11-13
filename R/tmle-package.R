

##' Targeted Maximum Likelihood Estimation for Binary Point Treatment Effects
##' 
##' Targeted maximum likelihood estimation of marginal treatment effect of a
##' binary point treatment on a continuous or binary outcome, adjusting for
##' baseline covariates. Missingness in the outcome is accounted for in the
##' estimation procedure. The population mean outcome is calculated when there
##' is missingness and no treatment.  Controlled direct effect estimation is
##' available. Optional data-adaptive estimation of \emph{Q} and \emph{g}
##' portions of the likelihood using the \code{SuperLearner} package is
##' strongly encouraged.
##' 
##' \tabular{ll}{ Package: \tab tmle\cr Type: \tab Package\cr Version: \tab
##' 1.1-0\cr Date: \tab 2010-12-20\cr License: \tab GPL-2 \cr LazyLoad: \tab
##' no\cr }
##' 
##' @name tmle-package
##' @docType package
##' @author Susan Gruber, in collaboration with Mark van der Laan.
##' 
##' Maintainer: Susan Gruber, \email{sgruber@@berkeley.edu}
##' @seealso \code{\link{tmle}}
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
NULL



