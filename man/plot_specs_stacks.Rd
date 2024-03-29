% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_specificities.r
\name{plot_specs_stacks}
\alias{plot_specs_stacks}
\alias{plot_specificities}
\title{plot_specs_stacks}
\usage{
plot_specs_stacks(
  specs_list,
  n_bins = 20,
  col_sig = "black",
  col_nsig = "gray",
  col_bord = NA,
  alpha = 0.05,
  label_cex = 1
)
}
\arguments{
\item{specs_list}{list of data.frames. Each data.frame must be an output from
phy_or_env_spec; must have columns "Spec" and "Pval".}

\item{n_bins}{integer. Number of bins for stacked violins (default: 20).}

\item{col_sig}{string. Color name or hex code for species where Pval <= alpha
(default: "black").}

\item{col_nsig}{string. Color name or hex code for species where Pval > alpha
(default: "gray").}

\item{col_bord}{string. Color name or hex code for border color. Use NA for no
border (default: NA).}

\item{alpha}{float. alpha value for determining statistical significance; see
col_sig and col_nsig above (default: 0.05).}

\item{label_cex}{float. Used to change size of x-axis labels (default: 1).}
}
\value{
returns nothing (a plot is made).
}
\description{
Visualizes results from phy_or_env_spec as stacked histograms. Aliased to 
plot_specificities() for backward compatibility.
}
\examples{
# library(specificity)
# attach(endophyte)
# # only analyze species with occupancy >= 20
# m <- occ_threshold(prop_abund(otutable), 20)
# # create list to hold phy_or_env_spec outputs
# specs_list <- list()
# specs_list$NDVI <- phy_or_env_spec(m, env=metadata$NDVI, 
#   n_cores=10, n_sim=50, p_method="gamma_fit")
# specs_list$Evapotranspiration <- phy_or_env_spec(m,
#   env=metadata$Evapotranspiration, n_cores=10, 
#   n_sim=100, p_method="gamma_fit")
# specs_list$Rainfall <- phy_or_env_spec(m, env=metadata$Rainfall,
#   n_cores=10, n_sim=50, p_method="gamma_fit")
# plot_specs_stacks(specs_list)

}
\author{
John L. Darcy
}
