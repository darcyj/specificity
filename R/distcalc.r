#' distcalc
#' 
#' Calculates pairwise geographic distance between locations on earth. Just a
#' convenient wrapper for fields::rdist.earth().
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
#' # library(specificity)
#' # attach(endophyte)
#' # geo_dists <- distcalc(metadata$Lat, metadata$Lon, metadata$SampleID)
#' # all(rownames(geo_dists) == metadata$SampleID)
#' 
#' @export
distcalc <- function(lat, lng, sampIDs=NULL){
	longlats <- data.frame(lng, lat)
	if(!is.null(sampIDs)){
		rownames(longlats) <- sampIDs
	}
	distmat=fields::rdist.earth(as.matrix(longlats), miles=F, R=NULL) 
	return((distmat))
}
