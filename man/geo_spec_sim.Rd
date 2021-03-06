% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geo_spec_sim.r
\name{geo_spec_sim}
\alias{geo_spec_sim}
\title{geo_spec_sim}
\usage{
geo_spec_sim(
  sdev,
  n_obs,
  grid,
  ideal_x = 0,
  ideal_y = 0,
  ideal_x2 = 0,
  ideal_y2 = 0,
  ideal_x3 = 0,
  ideal_y3 = 0,
  n_ideal = 1,
  up = 0,
  seed = 123456,
  n_cores = 2
)
}
\arguments{
\item{sdev}{numeric vector. Standard deviation of the probability distribution
P(species), in the same units as grid. P(species) is a function of the distance
between a sample site and its closest ideal location (specified with ideal_x/2/3
and ideal_y/2/3). Low values mean that the species is found in abundance within
only short distances of ideal locations, high values mean the species is found
across a wider area. Multiple values can be input in order to simulate a range
of specificities simultaneously. Can be length 1 or n.}

\item{n_obs}{integer vector. Number of observations to make, i.e. number of times species is
observed. Will be the sum of the species' output column. Can be length 1 or n.}

\item{grid}{data frame with columns x and y, representing cartesian coordinates
of sample locations. Can be artificial (generate with randomgrid()) or real.}

\item{ideal_x}{numeric vector. x-coordinate of the ideal spatial location for species (default:0).}

\item{ideal_y}{numeric vector. y-coordinate of the ideal spatial location for species (default:0).}

\item{ideal_x2}{numeric vector. x-coordinate for secondary ideal location. Only used if n_ideal<1 (default:0).}

\item{ideal_y2}{numeric vector. y-coordinate for secondary ideal location. Only used if n_ideal<1 (default:0).}

\item{ideal_x3}{numeric vector. x-coordinate for secondary ideal location. Only used if n_ideal<2 (default:0).}

\item{ideal_y3}{numeric vector. y-coordinate for secondary ideal location. Only used if n_ideal<2 (default:0).}

\item{n_ideal}{integer vector. Number of ideal locations to use. Must be 1, 2, or 3 (default:1).}

\item{up}{numeric vector. up=uniform proportion. This is the proportion of
the probability distribution P(species) that is composed of a uniform
distribution, if desired. If set to a value above zero (and blow 1), 
P(species) will be a weighted sum of the normal distribution described above,
and a uniform distribution. The weight for the uniform distribution will be
up, and the weight for the normal distribution will be 1-up (default: 0).}

\item{seed}{integer. Seed for randomization. Daughter seeds will be generated for
parallel computations, each with the same number of digits as seed 
(default: 1234567).}

\item{n_cores}{integer. Number of CPU cores for parallel computation (default: 2).}
}
\value{
List object containing "matrix" and "params" objects:
  \describe{
    \item{matrix}{
      matrix where each column is a vector of simulated observations for each row in grid;
      each column of matrix represents a simulated species.
    }
    \item{params}{
      data.frame of parameters (columns) used to simulate each species (rows).
    }
  }
}
\description{
Simulates inputs for phy_or_env_spec, by creating a species distribution over
artificial (or real) geographic space. That distribution has a bivariate mean
at the "ideal" location inspace for the simulated species, and the standard
deviation of that (normal) distribution controls the extent to which the 
species specific to geographic space. A high SD means less specificity, and a
low SD means more specificity.
}
\examples{
library(specificity)
g1 <- randomgrid()
plot(g1)
a1 <- geo_spec_sim(sdev=c(30, 30, 30, 30), n_obs=1000, grid=g1, up=c(0, 0.20, 0.40, 0.60))
par(mfrow=c(2,2))
plot_grid_abunds(g1, a1$matrix[,1])
plot_grid_abunds(g1, a1$matrix[,2])
plot_grid_abunds(g1, a1$matrix[,3])
plot_grid_abunds(g1, a1$matrix[,4])
a2 <- geo_spec_sim(sdev=c(10, 20, 30, 40), n_obs=1000, grid=g1, ideal_x=-50, ideal_x2=50, n_ideal=2)
par(mfrow=c(2,2))
plot_grid_abunds(g1, a2$matrix[,1], main="sd=10")
plot_grid_abunds(g1, a2$matrix[,2], main="sd=20")
plot_grid_abunds(g1, a2$matrix[,3], main="sd=30")
plot_grid_abunds(g1, a2$matrix[,4], main="sd=40")

}
\author{
John L. Darcy
}
