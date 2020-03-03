#' daughter_seeds
#' 
#' Makes n daughter seeds from seed s. This is useful for processes one wishes to be
#' deterministic, but may not be executed in the same order every time.
#' 
#' @author John L. Darcy
#' 
#' @param n integer. Number of daughter seeds to make.
#' @param s integer. A seed (DEFAULT: 12345).
#' 
#' @return vector of length n containing integer seeds.
#' 
#' @export
daughter_seeds <- function(n,s=12345){
	set.seed(s)
	return(replicate(n, as.integer(paste(sample(0:9, nchar(s), replace=TRUE), collapse=""))))
}
