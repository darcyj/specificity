#' multi_phy_or_env_spec
#'
#' Calculates standardized effect size (SES) of host phylogenetic specificity (per 
#' Poulin et al. 2011) or standardized effect size of environmental specificity
#' for each species (column) in a matrix. Arguments are mostly the same as 
#' phy_or_env_spec(), except abunds_mat must be a matrix or data frame where each 
#' column is a vector of species abundances corresponding to hosts. This function
#' can be run in parallel using the n_cores argument. 
#'
#' @author John L. Darcy
#' @references 
#' Poulin et al. (2011) Host specificity in phylogenetic and geographic
#'   space. Trends Parasitol 8:355-361. doi: 10.1016/j.pt.2011.05.003
#'
#' @param abunds_mat matrix or data frame of numeric values. Columns represent 
#'   species, rows are observations corresponding to hosts. Thus, each column is 
#'   the 'abunds' input to phy_or_env_spec().
#' @param hosts character vector. Host identities corresponding to abunds. Only
#' required if calculating SES for phylogenetic specificity.
#' @param hosts_phylo phylo object. Tree containing all unique hosts as tips. Only
#' required if calculating SES for phylogenetic specificity.
#' @param env numeric vector. Environmental variable corresponding to abunds. For
#'   example, temperature. Only required if calculating SES for environmental
#'   specificity.
#' @param unweighted logical. If TRUE, will use PD instead of WPD for phylogenetic
#'   specificity (DEFAULT: FALSE).
#' @param n_sim integer. number of simulations to do under the null hypothesis that
#'   host or environmental association is random (DEFAULT: 1000)
#' @param seed integer. Seed to use so that this is repeatable (DEFAULT: 12345)
#' @param tails integer. 1 = 1-tailed, test for specificity only. 2 = 2-tailed.
#'   3 = 1-tailed, test for cosmopolitanism only. 4 = no test, P=1.0 (DEFAULT: 1).
#' @param occ_stab logical. If TRUE, output statistic is stabilized for occupancy.
#'   if FALSE, not stabilized as in Poulin et al. 2011 (DEFAULT=TRUE).
#' @param occ_stab_frac float. Fractional occupancy between 0 and 1 for occupancy-
#'   stable calculation of denominator. Values between 0.10 and 0.40 work best
#'   (DEFAULT: 0.20).
#' @param force_denom float. If nonzero and occ_stab==TRUE, occupancy stabilization
#'   is done by forcing SES denominator to be a given value. Use 1 for non-
#'   standardized effect size (units of env or pd) or use 0 for denominator to be
#'   simulated using occ_stab_frac. Useful to make the code run faster by pre-
#'   computing the denominator if it's constant (DEFAULT: 0).
#'
#' @return data.frame where each row is an input species. First column is P-value,
#'   second column is SES.
#'
#' @examples
#' # none yet written.
#'
#' @export
	multi_phy_or_env_spec <- function(abunds_mat, hosts=NULL, hosts_phylo=NULL, 
		env=NULL, unweighted=FALSE, n_sim=1000, seed=12345, tails=1, occ_stab=TRUE, 
		occ_stab_frac=0.20, force_denom=0, n_cores=2){
		require("parallel")
		# set seed and generate a daughter seed for each column
		set.seed(seed)
		seeds <- replicate(n=ncol(abunds_mat), as.integer(paste(sample(0:9, 5, replace=TRUE), collapse="")))

		# if occ_stab, calculate denominator once so this can go faster
		if(occ_stab == TRUE){
			denom2force <- phy_or_env_spec(abunds=abunds_mat[,1], hosts=hosts,
			hosts_phylo=hosts_phylo, env=env, unweighted=unweighted, n_sim=n_sim, seed=seed, 
			tails=4, occ_stab=TRUE, occ_stab_frac=occ_stab_frac, 
			force_denom=force_denom, denom_only=TRUE)
		}else{
			denom2force <- 0
		}

		# wrapper function for parallelization - takes index as input
		col_inds <- 1:ncol(abunds_mat)
		phy_spec_wrapper <- function(i){ 
			phy_or_env_spec(abunds=abunds_mat[,i], hosts, hosts_phylo, env,
			unweighted, n_sim, seed=seeds[i], tails, occ_stab, occ_stab_frac,
			force_denom=denom2force, denom_only=FALSE)
		}

		results <- do.call(rbind, mclapply(X=col_inds, FUN=phy_spec_wrapper, mc.cores=n_cores))

		rownames(results) <- colnames(abunds_mat)
		results <- as.data.frame(results)
		return(results)

	}
