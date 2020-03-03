#' gamma_mode
#'
#' Estimates mode by fitting a gamma distribution.
#'
#' @author John L. Darcy
#'
#' @param x numeric vector of values
#' @param mincor minnimum acceptable correlation between binned frequencies of
#'   x (as in hist(x)) and the predicted probability densities for those bins
#'   from the fit gamma distribution. 
#' @param fallback function to use on x if gamma fails to fit. Default is mean.
#'   A warning will be given if fallback is used. 
#' @return mode of x, single value.
#'
#' @examples
#' none yet written.
#' 
#' @export
gamma_mode <- function(x, mincor=0.85, fallback=mean){
	require("fitdistrplus")

	output <- tryCatch({
		# fit gamma
	   	gfit <- suppressWarnings(fitdistrplus::fitdist(data=x, distr="gamma", lower=c(1,0)))
		shp <- gfit$estimate[names(gfit$estimate) == "shape"]
		rte <- gfit$estimate[names(gfit$estimate) == "rate"]

		# user-interpretable goodness of fit
		# basically, does fit dist's predictions correlate with hist(x)?
		hx <- hist(x, plot=F)
		hg <- dgamma(hx$mids, rate=rte, shape=shp)

		cxg <- cor(hx$counts, hg)

		# check if fit is good enough
		if(cxg < mincor){
			warning("Gamma fit not good enough, using fallback.")
		}

		list(warns=NULL, value=as.numeric((shp - 1) / rte))
		
	}, warning = function(w) {
		return(list(warn=w, value=fallback(x)))
	}, error = function(e) {
		return(list(warn="Gamma fit error, using fallback.", value=fallback(x)))
	}, finally = { })

	if(! is.null(output$warn)){
		warning(output$warn)
	}
	return(output$value)
}

