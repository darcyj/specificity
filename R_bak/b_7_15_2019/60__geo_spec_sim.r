#' geo_spec_sim
#'
#' Simulates inputs for phy_or_env_spec, by creating a species distribution over
#' artificial (or real) geographic space. That distribution has a bivariate mean
#' at the "ideal" location inspace for the simulated species, and the standard
#' deviation of that (normal) distribution controls the extent to which the 
#' species specific to geographic space. A high SD means less specificity, and a
#' low SD means more specificity. 
#'
#' @author John L. Darcy
#'
#' @param sdev numeric vector. Standard deviation of the probability distribution
#'   P(species), in the same units as grid. P(species) is a function of the distance
#'   between a sample site and its closest ideal location (specified with ideal_x/2/3
#'   and ideal_y/2/3). Low values mean that the species is found in abundance within
#'   only short distances of ideal locations, high values mean the species is found
#'   across a wider area. Multiple values can be input in order to simulate a range
#'   of specificities simultaneously. Can be length 1 or n.
#' @param n_obs number of observations to make, i.e. number of times species is
#'   observed. Will be the sum of the species' output column. Can be length 1 or n.
#' @param grid data frame with columns x and y, representing cartesian coordinates
#'   of sample locations. Can be artificial (generate with randomgrid()) or real.
#' @param ideal_x x-coordinate of the ideal spatial location for species (DEFAULT=0).
#' @param ideal_y y-coordinate of the ideal spatial location for species (DEFAULT=0).
#' @param ideal_x2 x-coordinate for secondary ideal location. Only used if n_ideal<1 (DEFAULT=0).
#' @param ideal_y2 y-coordinate for secondary ideal location. Only used if n_ideal<1 (DEFAULT=0).
#' @param ideal_x3 x-coordinate for secondary ideal location. Only used if n_ideal<2 (DEFAULT=0).
#' @param ideal_x3 x-coordinate for secondary ideal location. Only used if n_ideal<2 (DEFAULT=0).
#' @param n_ideal number of ideal locations to use. Must be 1, 2, or 3 (DEFAULT=1).
#' @param up numeric vector. up=uniform proportion. This is the proportion of
#'   the probability distribution P(species) that is composed of a uniform
#'   distribution, if desired. If set to a value above zero (and blow 1), 
#'   P(species) will be a weighted sum of the normal distribution described above,
#'   and a uniform distribution. The weight for the uniform distribution will be
#'   up, and the weight for the normal distribution will be 1-up (default: 0).
#' @param seed integer. Seed for randomization. Daughter seeds will be generated for
#'   parallel computations, each with the same number of digits as seed 
#'   (DEFAULT: 1234567).
#' @param n_cores integer. Number of CPU cores for parallel computation (DEFAULT: 2).
#'
#' @return List object containing "matrix" and "params" objects:
#'   \item{matrix}{
#'     matrix where each column is a vector of simulated observations for each row in grid;
#'    each column of matrix represents a simulated species.
#'   }
#'   \item{params}{
#'     data.frame of parameters (columns) used to simulate each species (rows).
#'   }
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
#'   a2 <- geo_spec_sim(sdev=c(10, 20, 30, 40), n_obs=1000, grid=g1, ideal_x=-50, ideal_x2=50, n_ideal=2)
#'   par(mfrow=c(2,2))
#'   plot_grid_abunds(g1, a2$matrix[,1], main="sd=10")
#'   plot_grid_abunds(g1, a2$matrix[,2], main="sd=20")
#'   plot_grid_abunds(g1, a2$matrix[,3], main="sd=30")
#'   plot_grid_abunds(g1, a2$matrix[,4], main="sd=40")
#'
#' @export
	geo_spec_sim <- function(sdev, n_obs, grid,
		ideal_x=0, ideal_y=0, ideal_x2=0, ideal_y2=0, ideal_x3=0, ideal_y3=0,
		n_ideal=1, up=0, seed=123456, n_cores=2){

		require(parallel)

		# for testing:
		if(FALSE){
			sdev=1:10; n_obs=1000; grid=randomgrid();
			ideal_x=0; ideal_y=0; ideal_x2=0; ideal_y2=0; ideal_x3=0; ideal_y3=0;
			n_ideal=1; up=0.10; seed=123456
		}

		# deal with variable inputs by constructing table for each simulated species
		# each row is the 4 parameters for a species.
		var_lens <- c(length(sdev), length(n_obs), length(ideal_x), length(ideal_y),
			length(ideal_x2), length(ideal_y2), length(ideal_x3), length(ideal_y3), length(n_ideal), length(up))
		if( ! all(var_lens) %in% c(1, max(var_lens))){
			stop("sdev, n_obs, ideal_x/y, n_ideal, and up must be either the length of the longest input variable or length 1.")
		}

		# check that all n_ideal values are in 1,2,3
		if(! all(n_ideal %in% 1:3)){
			stop("all n_ideal values must be either 1, 2, or 3.")
		}

		# generate n_sim daughter seeds
		set.seed(seed)
		seeds <- replicate(n=max(var_lens), as.integer(paste(sample(0:9, nchar(seed), replace=TRUE), collapse="")))
		seeds <- formatC(seeds, width=nchar(seed), format="d", flag="0")

		# mini 1-row df for each species that will be simulated.
		var_df <- data.frame(index=1:max(var_lens), sdev, n_obs,
			ideal_x, ideal_y, ideal_x2, ideal_y2, ideal_x3, ideal_y3, n_ideal, up, seed=seeds)
		var_list <- lapply(X=as.list(1:nrow(var_df)), FUN=function(x){ var_df[x,] })

		# function that takes parameters from var_list and generates counts
		geo_spec_sim_1species <- function(par, g=grid){
			# first, get differences from grid to ideal
			d2i <- sqrt( ((par$ideal_x - g$x)^2) + ((par$ideal_y - g$y)^2) )
			if(par$n_ideal == 2){
				d2i2 <- sqrt( ((par$ideal_x2 - g$x)^2) + ((par$ideal_y2 - g$y)^2) )
				d2i <- apply(X=rbind(d2i, d2i2), MAR=2, FUN=min)
			}else if(par$n_ideal == 3){
				d2i2 <- sqrt( ((par$ideal_x2 - g$x)^2) + ((par$ideal_y2 - g$y)^2) )
				d2i3 <- sqrt( ((par$ideal_x3 - g$x)^2) + ((par$ideal_y3 - generate$y)^2) )
				d2i <- apply(X=rbind(d2i, d2i2, d2i3), MAR=2, FUN=min)
			}

			# turn distances into probabilities
			unif_part <- rep(1/length(d2i), length(d2i)) * par$up
			norm_part <- dnorm(x=d2i, mean=0, sd=par$sdev)
			norm_part <- (norm_part / sum(norm_part)) * (1-par$up)
			probs <- unif_part + norm_part

			# sample using multinomial
			return(as.vector(rmultinom(n=1, size=par$n_obs, prob=probs)))
		}

		# apply geo_spec_sim_1species to all the pars

		out_matrix <- simplify2array(mclapply(X=var_list, FUN=geo_spec_sim_1species, mc.cores=n_cores))

		return(list(
			matrix=out_matrix, 
			params=var_df
		))
	}


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
