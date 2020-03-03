#' plot_grid_abunds
#'
#' plots species abundances across spatial sampling locations
#'
#' @author John L. Darcy
#' 
#' @param grid data frame with columns x and y, representing cartesian coordinates
#'   of sample locations. Can be artificial (generate with randomgrid()) or real.
#' @param abunds abundances of a species, corresponding to rows in grid.
#' @param pch pch character code to use for bottom of each abundance line (DEFAULT="")
#' @param ... arguments to be passed to plot.
#'
#' @return returns nothing, just makes a plot. 
#'
#' @examples
#'   g1 <- randomgrid()
#'   plot(g1)
#'   a1 <- geo_spec_sim(sdev=c(30, 30, 30, 30), n_obs=1000, grid=g1, up=c(0, 0.20, 0.40, 0.60))
#'   par(mfrow=c(2,2))
#'   plot_grid_abunds(g1, a1$matrix[,1])
#'   plot_grid_abunds(g1, a1$matrix[,2])
#'   plot_grid_abunds(g1, a1$matrix[,3])
#'   plot_grid_abunds(g1, a1$matrix[,4])
#'
#' @export
plot_grid_abunds <- function(grid, abunds, pch="", ...){
	plot(grid, pch=pch, ...)
	qtr <- (max(grid$x) - min(grid$x)) / 4
	heights <- (abunds / max(abunds)) * qtr
	segments(x0=grid$x, y0=grid$y, x1=grid$x, y1=grid$y+heights, col="blue", lwd=2)
}