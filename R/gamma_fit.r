#' gamma_fit
#'
#' Fits a gamma distribution to data. This function is really just a shorter, 
#' faster, far less felxible version of fitdist() from the "fitdistrplus"
#' package.
#'
#' @author John L. Darcy
#'
#' @param x numeric vector of values
#' @param lower_prop proportion of starting values to use for the lower bound
#'   of gamma distribution parameters (DEFAULT=0.10). Starting values are
#'   approximated as:
#'   \describe{
#'     \item{"shape":}{
#'       (m^2)/v
#'     }
#'     \item{"rate":}{
#'       m/v
#'     }
#'   }
#'   where m is the mean of x and v is the population variance of x. These values
#'   are multiplied by lower_prop in order to set shape and rate minnima for
#'   constrained optimization. The actual minima are c(0,0), and setting 
#'   lower_prop=0 will acheive this. But since the approximation is pretty
#'   good, using a lower_prop higher than 0 will give a faster optimization. 
#' @param upper_prop same as lower_prop, but for upper bound (Default = 10).
#' 
#' @return Vector of estimated values; c(shape, rate).
#'
#' @examples
#' none yet written.
#' 
#' @export
gamma_fit <- function(x, lower_prop=0.10, upper_prop=10){
	# get starting values
	n <- length(x)
	m <- mean(x)
	v <- (n-1)/n * var(x)
	start <- c(shape=(m^2)/v, rate=m/v)
	# calculate lower vals and upper vals
	lv <- start * abs(lower_prop)
	uv <- start * abs(upper_prop)
	# do optimization
	ans <- optim(
		par=start, 
		fn=function(par, obs){
			a <- dgamma(x=obs, shape=par[1], rate=par[2])
			# if optim tries to search for ridiculous values, it can make
			# zeroes. this just replaces them with the next best thing.
			#if(0 %in% a){
			#	a[a==0] <- min(a[a > 0])
			#}
			# next line for diagnostic purposes
			# print(paste(min(a), par[1], par[2]))
			#return( -1 * sum(log(a)) )
			return( -1 * sum(log(a + 1E-10)) )

		}, 
		method="L-BFGS-B", 
		lower=lv, 
		upper=uv, 
		obs=x
	)
	# just return estimate
	return(ans$par)
}
