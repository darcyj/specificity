% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_spec_and_pval.r
\name{calculate_spec_and_pval}
\alias{calculate_spec_and_pval}
\title{calculate_spec_and_pval}
\usage{
calculate_spec_and_pval(
  emp_raos,
  sim_raos,
  abunds_mat,
  env,
  p_adj = "fdr",
  tails = 1,
  n_cores = 2,
  verbose = TRUE,
  p_method = "raw",
  center = "mean",
  denom_type = "index",
  diagnostic = FALSE
)
}
\arguments{
\item{emp_raos}{vector. Empirical rao values, one per species ("feature").}

\item{sim_raos}{list of numeric vectors. Sim rao values, generated under null
hypothesis. Each item in list corresponds to an entry in emp_raos. As such,
length(emp_raos) must equal length(sim_raos). Each item within sim_raos is
a vector or rao values (length=n_sim in the case of phy_or_env_spec()).}

\item{abunds_mat}{site x species matrix. See ?phy_or_env_spec.}

\item{env}{MUST BE A dist OBJECT!!!! VERY IMPORTANT!!!! See ?phy_or_env_spec.}

\item{p_adj}{string. Type of multiple hypothesis testing correction performed
on P-values. Can take any valid method argument to p.adjust, including "none",
"bonferroni", "holm", "fdr", and others (default: "fdr").}

\item{tails}{integer. 1 = 1-tailed, test for specificity only. 2 = 2-tailed.
3 = 1-tailed, test for cosmopolitanism only. 0 = no test, P=1.0 (default: 1).}

\item{n_cores}{integer. Number of CPU cores to use for parallel operations. If
set to 1, lapply will be used instead of mclapply (default: 2).}

\item{verbose}{logical. Should status messages be displayed? (default: TRUE).}

\item{p_method}{string. "raw" for quantile method, or "gamma_fit" for calculating P
by fitting a gamma distribution (default: "raw").}

\item{center}{string. Type of central tendency to use for simulated RQE values.
Options are "mean", "median", and "mode". If mode is chosen, a reversible gamma
distribution is fit and mode is calculated using that distribution (default: mean).}

\item{denom_type}{string. Type of denominator (d) to use (default: "index"). Note
that denominator type does NOT affect P-values.
\describe{
  \item{"ses":}{
    d for species s is calculated as the standard deviation of RQE values
    calculated from permuted species weights. This makes the output specificity
    a standardized effect size (SES). Unfortunately, this makes SES 
    counterintuitively sensitive to occupancy, where species with high occupancy
    have more extreme SES than rare species, due to their more deterministic sim
    specificities. Included for comparative purposes, not suggested.
  }
  \item{"raw":}{
    d is 1 for all species, so output specificity has units of distance, i.e. the
    raw difference between empirical and simulated RQE. This means that results
    from different variables are not comparable, since it is not scale-invariant to
    env or hosts_phylo. It not scale-invariant to the species weights in aunds_mat,
    either. Not sensitive to number of samples. Not suggested because units are
    strange, and isn't comparable between variables. 
  }
  \item{"index":}{
    d is the mean of simulated (permuted) RQE values for species that have stronger
    specificity than expected by chance, resulting in specificity values with range
    [-1, 0), with 0 as the null hypothesis. In this case, -1 indicates perfect
    specificity, where a species is associated with zero environmental variability.
    In the euclidean sense, this could be a species that is always found at the
    exact same elevation or the exact same pH. For species that have weaker specificity
    than expected by chance, d is x minus the center (see above) of simulated RQE 
    values, where x is the maximum possible dissimilarity observable given species
    weights. This d has other useful properties: scale invariance to env/hosts_phylo,
    insensitivity to the number of samples, insensitivity to occupancy, and strong 
    sensitivity to specificity (default).
  }
}}

\item{diagnostic}{logical. If true, changes output to include different parts of SES. 
This includes Pval, SES, raw, denom, emp, and all sim values with column labels as
simN where N is the number of sims (default: FALSE)}
}
\value{
data.frame where each row is an input species. First column is P-value
  ($Pval), second column is specificity ($Spec).
}
\description{
This function is called by phy_or_env_spec(). It is made available as a
standalone function in the (rare) case a user wishes to calculate specificity
using their own null model. calculate_spec_and_pval() takes empirical rao
values and sim rao values (from a null model) and calculates specificity and
P-values. To do that, use your own null model to make species data, and use
rao1sp() and/or raoperms() to get raw rao values. This function expects a
vector of empirical values, and a list of vectors of sim values (see below).
Most of the inputs for this function are the same as phy_or_env_spec(). Think 
of this function as the final component of "build your own phy_or_env_spec()".
Note that for this custom approach, the environmental variable must be a dist.
}
\examples{
# None yet. Forthcoming examples:
# 1. calculating regular old elevational specificity the hard way
# 2. same thing, but using vazquez null model from bipartite package

}
\author{
John L. Darcy
}
