#' randomgrid
#'
#' Generates a random spatial sampling using a bivariate random uniform distribution.
#'
#' @author John L. Darcy
#'
#' @param n_samp number of sampling locations to output (DEFAULT=1000).
#' @param xmin minimum x-axis coordinate (DEFAULT=-100).
#' @param xmax maximum x-axis coordinate (DEFAULT=100).
#' @param ymin minimum y-axis coordinate (DEFAULT=-100).
#' @param ymax maximum y-axis coordinate (DEFAULT=100).
#' @param seed integer, seed for randomization.
#'
#' @return data.frame object with x and y columns, with n_samp rows.
#'
#' @examples
#'   g <- randomgrid()
#'   plot(g)
#'   g2 <- randomgrid(nsamp=50, xmin=0, ymin=0)
#'   plot(g2)
#'
#' @export
randomgrid <- function(n_samp=1000, xmin=-100, xmax=100, ymin=-100, ymax=100, seed=123456){
	set.seed(seed)
	x <- runif(n=n_samp, min=xmin, max=xmax)
	y <- runif(n=n_samp, min=ymin, max=ymax)
	return(data.frame(x,y))
}
