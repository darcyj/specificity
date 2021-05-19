#' random_rep_positions
#' 
#' Finds positions in a vector (or matrix) that are randomly located within n_bins
#' evenly sized bins. This is useful for 1:1 comparisons of large vectors where
#' plotting or comparing all points is prohibitive. Only used in an example for the
#' prop_abund() function.
#' 
#' @author John L. Darcy
#' 
#' @param x vector
#' @param nbins number of bins to use
#' 
#' @return integer vector of positions that were selected
#' 
#' @examples
#' # library(specificity)
#' # attach(endophyte)
#' # m_dbl <- prop_abund(otutable)
#' # m_int <- prop_abund(otutable, to_int=TRUE)
#' # head(rowSums(m_dbl))
#' # head(rowSums(m_int))
#' # # note that they are off by a little bit. This small loss in precision is OK.
#' # object.size(m_dbl)
#' # object.size(m_int)
#' # random_positions <- random_rep_positions(m_dbl, 100)
#' # plot(m_int[random_positions] ~ m_dbl[random_positions])
#' 
#' @export
random_rep_positions <- function(x, nbins=50){
	breaks <- seq(from=min(x), to=max(x), length.out=50)
	output <- rep(0, nbins - 1)
	for(i in 1:(nbins -1)){
		bin_TF <- ( x >= breaks[i] & x <= breaks[i+1] )
		bin_pos <- (1:length(x))[bin_TF]
		output[i] <- sample(bin_pos, 1)
	}
	return(output)
}