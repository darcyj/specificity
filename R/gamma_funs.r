
# function to fit a gamma distribution to data
# returns a list with par, log likelihood, and type - forward or reverse. for reverse,
# also returns mx which is max(x), needed for reversing quantiles.
fit_gamma_fwd_rev <- function(x, lower_prop=0.10, upper_prop=10){
	# lower_prop <- 0.10
	# upper_prop <- 10
	n<-length(x); xc<-x-mean(x); 
	skewx <- (sqrt(n) * sum(xc^3)/(sum(xc^2)^(3/2))) * ((1 - 1/n))^(3/2) 

	if(is.null(skewx) || is.na(skewx) || is.nan(skewx)){
		# likely case where skew can't be calculated because var=0.
		# return fit that tells downstream functions to use meaan(x)
		# instead of mode (because they are the same) and to use
		# quantile method for p (because p=0 or 1).
		output <- (list(type="nofit"))
	}else if(skewx >= 0){ # if positive skew, fit regular gamma.

		output <- tryCatch({
			# try:
			# get starting values
			n <- length(x)
			m <- mean(x)
			v <- (n-1)/n * var(x)
			start <- c(shape=(m^2)/v, rate=m/v)
			# calculate lower vals and upper vals
			lv <- start * abs(lower_prop)
			uv <- start * abs(upper_prop)
			# do optimization
			fit <- optim(
				par=start, 
				fn=function(par){
					if(par[1] <= 0){par[1] <- 1E-10} # soft boundry for nm
					if(par[2] <= 0){par[2] <- 1E-10} # soft boundry for nm
					a <- dgamma(x, shape=par[1], rate=par[2])
					return( -1 * sum(log(a + 1E-10)) )

				}, 
				method="Nelder-Mead"
				#method="L-BFGS-B",  # nm is twice as fast.
				#lower=lv, 
				#upper=uv
			)
			# just return estimate
			list(
				par=fit$par,
				ll=-1 * fit$value,
				type="forward"
			)

		}, error=function(cond){
			return(list(type="nofit"))
		}, warning=function(cond){
			return(list(type="nofit"))
		})

	}else if(skewx < 0){  # if negative skew, fit reversed gamma.

		output <- tryCatch({
			# try:
			# reverse x
			mx <- max(x)
			mx <- mx + 0.001 * mx # add a little bit to make sure we never have 0 values
			rx <- mx - x
			n <- length(rx)
			v <- ((n-1)/n) * var(rx)
			
			# get starting values as gamma mles
			start <- c(shape=(mean(rx)^2)/v, rate=mean(rx)/v)
			fit <- optim(
				par=start,
				fn <- function(par){
					if(par[1] <= 0){par[1] <- 1E-10} # soft boundry for nm
					if(par[2] <= 0){par[2] <- 1E-10} # soft boundry for nm
					a <- dgamma(rx, shape=par[1], rate=par[2])
					return(-1 * sum(log(a+1E-10)))
				},
				method="Nelder-Mead"
				#method="L-BFGS-B",  # nm is twice as fast.
				#lower=lv, 
				#upper=uv
			)
			list(
				par=fit$par,
				ll=-1 * fit$value,
				mx=mx,
				type="reverse"
			)
		}, error=function(cond){
			return(list(type="nofit"))
		}, warning=function(cond){
			return(list(type="nofit"))
		})

	}else{
		# likely case where skew can't be calculated because var=0.
		# return fit that tells downstream functions to use meaan(x)
		# instead of mode (because they are the same) and to use
		# quantile method for p (because p=0 or 1).
		output <- (list(type="nofit"))
	}

	# final sanity check
	xhcor <- gamma_hist_cor(x, output)
	if(xhcor < 0.70){
		output <- (list(type="nofit"))		
	}
	output$xhcor <- xhcor

	return(output)
}



# user-interpretable check of goodness of fit for distribution
# x is perms, fit is result of fit_gamma_fwd_rev
gamma_hist_cor <- function(x, fit){
	# because of bug in hist() with fake numbers, here is a fake hist.
	fakehist <- function(x){
		histtab <- table(cut(x, breaks=ceiling(1+3.322*log10(length(x)))))
		nm <- gsub("[()]", x=gsub(pattern="\\]", x=names(histtab), ""), "")
		mids <- sapply(
			X=lapply(X=strsplit(nm, ","), FUN=as.numeric), 
			FUN=function(x){x[1] + 0.5 * (x[2] - x[1])}
		)
		return(list(mids=mids, density=as.numeric(histtab)))
	}
	if(fit$type == "forward"){
		hx <- fakehist(x)
		dx <- dgamma(hx$mids, fit$par[1], fit$par[2])
		return(cor(hx$density, dx))
	}else if(fit$type == "reverse"){
		rx <- fit$mx - x
		hrx <- fakehist(rx)
		drx <- dgamma(hrx$mids, fit$par[1], fit$par[2])
		return(cor(hrx$density, drx))
	}else{
		return(0)
	}
}


# calculates mode for a given gamma distribution fit
# can handle fwd or reverse fits.
# fit is just the object returned by fit_gamma_fwd_rev.
# fallback is what to return if fit failed
mode_gamma_fwd_rev <- function(fit, fallback){
	mo <- (fit$par[1] - 1) / fit$par[2]
	if(fit$type == "reverse"){
		# reverse mo
		mo <- fit$mx - mo 
	}else if(fit$type=="forward"){
		return(mo)
	}else{
		return(fallback)
	}
}


# same as pgamma but for custom fit objects returned by fit_gamma_fwd_rev.
# q - a quantile (in non-reversed orientation)
# fit - an object from fit_gamma_fwd_rev\
# fallback - what to return if fit failed
p_gamma_fwd_rev <- function(q, fit, fallback){
	if(fit$type == "reverse"){
		q_r <- fit$mx - q  # reverse the quantile, then plug in
		return(1 - pgamma(q_r, fit$par[1], fit$par[2]))
	}else if(fit$type=="forward"){
		return( pgamma(q, fit$par[1], fit$par[2]) )
	}else{
		return( fallback )
	}
}

# calculates a p-value either for permutations or from a gfit object (from fit_gamma_fwd_rev )
# emp - an empirical value
# perm - values from permutation
# gfit - gamma fit from fit_gamma_fwd_rev (use this OR perm, not both)
# tails - 1=left, 2=both, 3=right
# threshold - minimum number of perm values to return a P
# fallback - p to use instead of gamma p if gamma fit didn't work
p_from_perms_or_gfit <- function(emp, perm=NULL, gfit=NULL, fallback=1, 
	tails=1, threshold=30){

	if(!is.null(perm)){
		n <- sum(!is.na(perm))
		if(n < threshold){
			LTP <- NA
			RTP <- NA
		}else{
			# remove NAs
			perm <- perm[!is.na(perm)]
			# calculate P for right and left tails
			LTP <- (sum(perm <= emp) + 1)/n
			RTP <- 1 - LTP
		}
	}else if(!is.null(gfit)){
		LTP <- p_gamma_fwd_rev(emp, gfit, fallback)
		RTP <- 1 - LTP
	}else{
		stop("Invalid p-value method argument.")
	}

	# fix P > 1 in cases where all values of perm == emp
	if(LTP > 1){LTP <- 1}
	if(RTP > 1){RTP <- 1}
	# calculate P for tailed situations
	if(tails == 1){
		return(LTP)
	}else if(tails == 3){
		return(RTP)
	}else if(tails == 2){
		return(min(c(LTP, RTP)))
	}else{
		stop("Invalid tails argument.")
	}
}
