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
#'   set.seed(12345)
#'   # example that should work:
#'   a <- (rgamma(100, 12, 1))
#'   hist(a)
#'   gamma_mode(a)
#'   abline(v=gamma_mode(a), col="red")
#'   # example that should NOT work:
#'   a <- runif(100, 0, 1)
#'   gamma_mode(a)
#'   
#' 
#' @export
gamma_mode <- function(x, mincor=0.85, fallback=mean){

	output <- tryCatch({
		# fit gamma
		gfit <- gamma_fit(x)
		shp <- gfit[names(gfit) == "shape"]
		rte <- gfit[names(gfit) == "rate"]

		# user-interpretable goodness of fit
		# basically, does fit dist's predictions correlate with hist(x)?
		hx <- hist(x, plot=F)
		hg <- dgamma(hx$mids, rate=rte, shape=shp)

		cxg <- cor(hx$counts, hg)

		# check if fit is good enough
		if(cxg < mincor){
			wn <- paste("Gamma fit not good enough, using fallback:", deparse(substitute(fallback)))
			warning(wn)
		}

		list(warns=NULL, value=as.numeric((shp - 1) / rte))
		
	}, warning = function(w) {
		return(list(warn=w, value=fallback(x)))
	}, error = function(e) {
		wn <- paste("Gamma fit error, using fallback:", deparse(substitute(fallback)))
		return(list(warn=wn, value=fallback(x)))
	}, finally = { })

	if(! is.null(output$warn)){
		warning(output$warn)
	}
	return(output$value)
}

