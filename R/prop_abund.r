#' prop_abund
#'
#' Calculates proportional abundance of each species (columns) across samples
#' (rows) in community data matrix m. Row sums of output matrix will all be 1. 
#'
#' @author John L. Darcy
#'
#' @param m matrix or data frame of numeric values. Columns represent 
#'   species, rows are samples. 
#' @param to_int logical. Should output matrix be transformed into integers from
#'   0 to max_int? Integers take up half as much space as doubles, and as weights
#'   are equivalent for calculating specificity. The tradeoff is a little bit of
#'   precision (default: FALSE).
#' @param max_int integer. Maximum integer value used for to_int. If pairwise
#'   geometric means will be calculated with these data, it is nice to keep this
#'   value as the square root of the maximum integer size, which is the default.
#' @param speciesRows logical. Do rows represent species (instead of samples)? 
#'   (default:FALSE)
#' @return matrix of proportional abundances.
#'
#' @examples
#' library(specificity)
#' attach(endophyte)
#' m_dbl <- prop_abund(otutable)
#' m_int <- prop_abund(otutable, to_int=TRUE)
#' head(rowSums(m_dbl))
#' head(rowSums(m_int))
#' # note that they are off by a little bit. This small loss in precision is OK.
#' object.size(m_dbl)
#' object.size(m_int)
#' random_positions <- random_rep_positions(m_dbl, 100)
#' plot(m_int[random_positions] ~ m_dbl[random_positions])
#' 
#' @export
prop_abund <- function(m, to_int=FALSE, max_int=floor(sqrt(.Machine$integer.max)), speciesRows=FALSE){
	if(speciesRows){m <- t(m)}
	totals <- matrix(rep(rowSums(m), ncol(m)), nrow=nrow(m))
	m_out <- m / totals
	if(to_int){
		m_out <- matrix(as.integer(m_out * max_int), nrow=nrow(m_out))
		dimnames(m_out) <- dimnames(m)
	}
	if(speciesRows){m_out <- t(m_out)}
	return(m_out)
}