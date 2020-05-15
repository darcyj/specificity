#' pval_from_perms
#'
#' Calculates P-value for permutation tests. 
#'
#' @author John L. Darcy
#'
#' @param emp Numeric scalar. An empirical test statistic value. 
#' @param perm Numeric vector. Test statistic values similar to emp, but calculated
#'   from permuted data. 
#' @param tails integer. 
#'   \describe{
#'     \item{1:}{Left tail only.}
#'     \item{2:}{2-tailed test.}
#'     \item{3:}{Right tail only.}
#'     \item{0:}{No test, P=1.}
#'   }
#' @param method string. Method by which P should be calculated from perms:
#'   \describe{
#'     \item{"raw":}{
#'       P is calculated as the sum of sim values more extreme than the empirical
#'       value plus one, divided by the number of sim values. 
#'     }
#'     \item{"dens_fit":}{
#'       P is calculated via kernel density estimation. No better than "raw".
#'     }
#'     \item{"gamma_fit":}{
#'       P is calculated by fitting a gamma distribution to sim values and calculating 
#'       area under the curve from (-inf,emp] or [emp,inf) depending on tailedness.
#'     }
#'   }
#' @param threshold integer. Minimum number n of non-NA values in perm that are
#'   acceptable. If n < threshold, P=NA (DEFAULT: 50).
#' @param rounding integer. Number of decimal places to round emp and perm This is only
#'   useful when emp and perm are expected to contain the exact same value, but the number
#'   of decimal places in that value is different between emp and perm Use a number less
#'   then zero to disable rounding (DEFAULT: -1).
#'
#' @return a P-value.
#' 
#' @export
pval_from_perms <- function(emp, perm, tails, method="raw", threshold=30, rounding=-1){
	# do rounding...?
	if(rounding >= 0){
		emp <- round(emp, rounding)
		perm <- round(perm, rounding)
	}
	# check if number of perms is below threshold
	n <- sum(!is.na(perm))
	if(n < threshold){
		LTP <- NA
		RTP <- NA
	}else if(method == "dens_fit"){
		# remove NAs
		perm <- perm[!is.na(perm)]
		frm <- min(c(mean(perm) - (mean(perm) - min(perm)) * 4, emp))
		too <- max(c(mean(perm) + (max(perm) - mean(perm)) * 4, emp))
		dns <- density(perm, from=frm, to=too, bw="sj")
		# Riemann sum
		dx <- dns$x[2] - dns$x[1]
		tot <- sum(dns$y * dx)
		LTP <- sum(dns$y[dns$x < emp] * dx) / tot
		RTP <- 1-LTP
	}else if(method == "raw"){
		# remove NAs
		perm <- perm[!is.na(perm)]
		# calculate P for right and left tails
		LTP <- (sum(perm <= emp) + 1)/n
		RTP <- (sum(perm >= emp) + 1)/n
	}else if(method == "gamma_fit"){
		gfit <- gamma_fit(perm)
		shp <- gfit[1]
		rte <- gfit[2]
		LTP <- pgamma(q=emp, shape=shp, rate=rte)
		RTP <- 1 - LTP
	}else{
		stop("Invalid method argument.")
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
