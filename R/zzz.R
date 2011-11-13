.onAttach = function(libname, pkgname) {
        cat("Welcome to the tmle package, version 1.1\n")
        cat("This package is provided for academic use only.\n")
        cat("\nData-adaptive prediction algorithms are required to obtain reliable results.\nSuperLearner and DSA packages are strongly recommended.\n")
        cat("\nSuperLearner is available at http://www.stat.berkeley.edu/users/ecpolley/SL/\n")
 cat("DSA: http://www.stat.berkeley.edu/users/laan/Software/\n")
 cat("These packages depend on
 \tmodelUtils:  http://www.stat.berkeley.edu/users/laan/Software/
 \tnnls:        http://cran.r-project.org/web/packages/nnls/index.html
 \tquadprog:    http://cran.r-project.org/web/packages/quadprog/index.html\n")
 cat("\nNew in version 1.1:\n")
 cat("\t Population mean  parameter calculated when there is missing data, no variability in treatment assignment\n") 
 cat("\t Poisson regression available for count data with user-supplied regression formula\n")
 cat("\t Binary mediating variable allowed\n")
 cat("\t Bug in missingness mechanism fixed\n")
 cat("\t Some arguments and return values have been modified. See help files for details\n")
}
