#' mat_prop
#'
#' mat_prop(m) takes matrix m and transforms it such that rowsums are 1 (unless
#' the row had a total of zero, in which case it will be left as all zeroes).
#' 
#'
#' @author John L. Darcy
#' @param m matrix or data frame of numeric values. Columns represent species, rows
#'   are observations corresponding to hosts. Negative values forbidden.
#'
#' @return matrix of same dimensions as m
#'
#' @examples
#'  set.seed(12345)
#'	m <- matrix(sample(c(0,1), 12, replace=T), nrow=3, ncol=4)
#'  mat_prop(m)
#'  rowSums(mat_prop(m))
#'
#' @export
	mat_prop <- function(m){
		m <- as.matrix(m)
		prop <- function(x){
			if(sum(x) != 0){
				return(x/sum(x))
			}else{
				return(x)
			}
		}
		return( t(apply(X=m, FUN=prop, MAR=1)) )
	}
#'
#' mat_passthru
#'
#' mat_passthru(m) simply returns matrix m, unaltered.
#' 
#'
#' @author John L. Darcy
#' @param m matrix or data frame of numeric values. Columns represent species, rows
#'   are observations corresponding to hosts. Negative values forbidden.
#'
#' @return matrix of same dimensions as m
#'
#' @examples
#'  set.seed(12345)
#'	m <- matrix(sample(c(0,1), 12, replace=T), nrow=3, ncol=4)
#'  mat_passthru(m)
#'  mat_passthru(m) == m
#' 
#' @export
	mat_passthru <- function(m){ 
		return( as.matrix(m) )
	}













