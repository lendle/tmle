#---------- function .bound ---------------
# set outliers to min/max allowable values
# assumes x contains only numerical data
#-----------------------------------------
.bound <- function(x, bounds){
	x[x>max(bounds)] <- max(bounds)
	x[x<min(bounds)] <- min(bounds)
	return(x)
}
