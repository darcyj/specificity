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
#'   precision (DEFAULT: FALSE).
#' @param max_int integer. Maximum integer value used for to_int. If pairwise
#'   geometric means will be calculated with these data, it is nice to keep this
#'   value as the square root of the maximum integer size, which is the default.
#' @param speciesRows logical. Do rows represent species (instead of samples)? 
#'   (DEFAULT:FALSE)
#' @return matrix of proportional abundances.
#'
#' @examples
#' attach(endophyte)
#' m_dbl <- prop_abund(zotutable)
#' m_int <- prop_abund(zotutable, to_int=TRUE)
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
#'   considered absent (DEFAULT: 0). 
#' 
#' @return matrix with low-occupancy species removed.
#' 
#' @examples
#' attach(endophyte)
#' dim(zotutable)
#' zotutable_over25 <- occ_threshold(zotutable, 25)
#' dim(zotutable_over25)
#' 
#' @export 
occ_threshold <- function(m, threshold, max_absent=0){
	occs <- colSums(m > max_absent)
	goodspecies <- occs >= threshold
	return(m[,goodspecies])
}




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


#' geo_distcalc
#' 
#' Calculates pairwise geographic distance between locations on earth. Just a
#' convenient wrapper for fields::rdist.earth. 
#' 
#' @author John L. Darcy
#' 
#' @param lat Numeric vector. Latitudes in decimal degree format.
#' @param lng Numeric vector. Longitudes in decimal degree format.
#' @param sampIDs Character vector. Sample identifiers. Only required
#'   if output dist should have names associated.
#' 
#' @return matrix containing all pairwise geographic distances in km. 
#' 
#' @examples
#' data(endophyte)
#' geo_dists <- distcalc(metadata$Lat, metadata$Lon, metadata$SampleID)
#' all(rownames(geo_dists) == metadata$SampleID)
#' 
#' @export
distcalc <- function(lat, lng, sampIDs=NULL){
	require(fields)
	longlats <- data.frame(lng, lat)
	if(!is.null(sampIDs)){
		rownames(longlats) <- sampIDs
	}
	distmat=rdist.earth(as.matrix(longlats), miles=F, R=NULL) 
	return((distmat))
}


