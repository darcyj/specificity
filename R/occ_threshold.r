#' occ_threshold
#' 
#' removes species (columns) from a matrix that don't meet a minimum occupancy, 
#' defined as the number of samples in which that species was observed. 
#' 
#' @author John L. Darcy
#' 
#' @param m matrix or data frame of numeric values. Columns represent 
#'   species, rows are samples. 
#' @param threshold integer. Minimum number of samples a species can occupy
#'   without being removed.
#' @param max_absent float. Maximum abundance value at which a species will be
#'   considered absent (default: 0). 
#' 
#' @return matrix with low-occupancy species removed.
#' 
#' @examples
#' library(specificity)
#' attach(endophyte)
#' dim(otutable)
#' otutable_over25 <- occ_threshold(otutable, 25)
#' dim(otutable_over25)
#' 
#' @export 
occ_threshold <- function(m, threshold, max_absent=0){
	occs <- colSums(m > max_absent)
	goodspecies <- occs >= threshold
	return(m[,goodspecies])
}