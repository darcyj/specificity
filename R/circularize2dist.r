
#' circularize2dist
#'
#' Circularizes a vector into a dist object. For example, a vector of days
#' of the year, where the distance between 365 and 2 should be less than the 
#' distance between 350 and 365. Another example may be direction, where 0.1 and
#' 2pi radians are close together.
#'
#' @author John L. Darcy
#'
#' @param x a numeric vector. All values should be >0. 
#' @param maxx the maximum theoretical value (also the zero value!) of variable x.
#'   In the exampe of months of the year, maxx would be 12, even if you only had
#'   data for months 1-8. For degrees, maxx=360. For radians, maxx=2*pi. Must be
#'   greater than or equal to values of x.
#'
#' @return a vector of differences, ordered identically to a "dist" object.
#'
#' @examples
#' library(specificity)
#' # make some fake data to represent months of the year
#' months <- c(1, 4, 11)
#' # run circularize2dist() on the months. Must specify that
#' # maxx = 12, since december is both 12 and 0 for these data.
#' circularize2dist(months, 12)
#' # output is a distance matrix. 
#' # rows and cols of months_circdm are months - it's ordered.
#' # notice the distance between 11 and 1 is 2, not 10!
#'
#' @export
circularize2dist <- function(x, maxx){
	if(!is.numeric(x)){stop("x must be numeric")}
	if(!is.atomic(x) || !is.vector(x)){stop("x must be a 1-dimensional atomic vector")}
	if(!is.numeric(maxx)){stop("maxx must be numeric")}
	if(length(maxx) != 1){stop("maxx must be scalar")}
	if(any(x > maxx)){stop("some values of x are greater than maxx")}
	half <- maxx / 2
	cdis <- function(x1, x2){
		dis <- abs(x1 - x2)
		if(dis > half){
			return(abs(maxx - dis))
		}else{
			return(dis)
		}
	}
	return(as.dist(sapply(x, function(x1){sapply(x, function(x2){cdis(x1,x2)})})))
}