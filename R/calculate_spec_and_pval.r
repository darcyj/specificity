#' calculate_spec_and_pval
#'
#' This function is called by phy_or_env_spec(). It is made available as a
#' standalone function in the (rare) case a user wishes to calculate specificity
#' using their own null model. calculate_spec_and_pval() takes empirical rao
#' values and sim rao values (from a null model) and calculates specificity and
#' P-values. To do that, use your own null model to make species data, and use
#' rao1sp() and/or raoperms() to get raw rao values. This function expects a
#' vector of empirical values, and a list of vectors of sim values (see below).
#' Most of the inputs for this function are the same as phy_or_env_spec(). Think 
#' of this function as the final component of "build your own phy_or_env_spec()".
#' Note that for this custom approach, the environmental variable must be a dist.
#'
#' @author John L. Darcy
#'
#' @param emp_raos vector. Empirical rao values, one per species ("feature").
#' @param sim_raos list of numeric vectors. Sim rao values, generated under null
#'   hypothesis. Each item in list corresponds to an entry in emp_raos. As such,
#'   length(emp_raos) must equal length(sim_raos). Each item within sim_raos is
#'   a vector or rao values (length=n_sim in the case of phy_or_env_spec()).
#' @param abunds_mat site x species matrix. See ?phy_or_env_spec.
#' @param env MUST BE A dist OBJECT!!!! VERY IMPORTANT!!!! See ?phy_or_env_spec.
#' @param p_adj string. Type of multiple hypothesis testing correction performed
#'   on P-values. Can take any valid method argument to p.adjust, including "none",
#'   "bonferroni", "holm", "fdr", and others (default: "fdr").
#' @param tails integer. 1 = 1-tailed, test for specificity only. 2 = 2-tailed.
#'   3 = 1-tailed, test for cosmopolitanism only. 0 = no test, P=1.0 (default: 1).
#' @param n_cores integer. Number of CPU cores to use for parallel operations. If
#'   set to 1, lapply will be used instead of mclapply (default: 2).
#' @param verbose logical. Should status messages be displayed? (default: TRUE).
#' @param p_method string. "raw" for quantile method, or "gamma_fit" for calculating P
#'   by fitting a gamma distribution (default: "raw").
#' @param center string. Type of central tendency to use for simulated RQE values.
#'   Options are "mean", "median", and "mode". If mode is chosen, a reversible gamma
#'   distribution is fit and mode is calculated using that distribution (default: mean).
#' @param denom_type string. Type of denominator (d) to use (default: "index"). Note
#'   that denominator type does NOT affect P-values.
#'   \describe{
#'     \item{"ses":}{
#'       d for species s is calculated as the standard deviation of RQE values
#'       calculated from permuted species weights. This makes the output specificity
#'       a standardized effect size (SES). Unfortunately, this makes SES 
#'       counterintuitively sensitive to occupancy, where species with high occupancy
#'       have more extreme SES than rare species, due to their more deterministic sim
#'       specificities. Included for comparative purposes, not suggested.
#'     }
#'     \item{"raw":}{
#'       d is 1 for all species, so output specificity has units of distance, i.e. the
#'       raw difference between empirical and simulated RQE. This means that results
#'       from different variables are not comparable, since it is not scale-invariant to
#'       env or hosts_phylo. It not scale-invariant to the species weights in aunds_mat,
#'       either. Not sensitive to number of samples. Not suggested because units are
#'       strange, and isn't comparable between variables. 
#'     }
#'     \item{"index":}{
#'       d is the mean of simulated (permuted) RQE values for species that have stronger
#'       specificity than expected by chance, resulting in specificity values with range
#'       [-1, 0), with 0 as the null hypothesis. In this case, -1 indicates perfect
#'       specificity, where a species is associated with zero environmental variability.
#'       In the euclidean sense, this could be a species that is always found at the
#'       exact same elevation or the exact same pH. For species that have weaker specificity
#'       than expected by chance, d is x minus the center (see above) of simulated RQE 
#'       values, where x is the maximum possible dissimilarity observable given species
#'       weights. This d has other useful properties: scale invariance to env/hosts_phylo,
#'       insensitivity to the number of samples, insensitivity to occupancy, and strong 
#'       sensitivity to specificity (default).
#'     }
#'   }
#' @param diagnostic logical. If true, changes output to include different parts of SES. 
#'   This includes Pval, SES, raw, denom, emp, and all sim values with column labels as
#'   simN where N is the number of sims (default: FALSE)
#'
#'
#' @return data.frame where each row is an input species. First column is P-value
#'   ($Pval), second column is specificity ($Spec).
#'
#' @examples
#' # None yet. Forthcoming examples:
#' # 1. calculating regular old elevational specificity the hard way
#' # 2. same thing, but using vazquez null model from bipartite package
#'
#' @export
calculate_spec_and_pval <- function(emp_raos, sim_raos, abunds_mat, env,
	p_adj="fdr", tails=1, n_cores=2, verbose=TRUE, p_method="raw", center="mean", 
	denom_type="index", diagnostic=FALSE){

	# little message function, to clean up code a bit
	msg <- function(x){if(verbose){message(x)}}

	# function to pretty format time stuff or whatever
	fms <- function(x){
		secs <- x[3]
		if(secs < 60){
			return(paste(round(secs, 1), "seconds"))
		}else if(secs < 60*60){
			return(paste(round(secs/60, 1), "minutes"))
		}else{
			return(paste(round(secs/(60*60), 1), "hours"))
		}
	}



	# apply function for parallelism or not:
	if(n_cores > 1){
		lapply_fun <- function(X, FUN, ...){parallel::mclapply(X, FUN, mc.cores=n_cores, ...)}
		mapply_fun <- function(FUN, ...){parallel::mcmapply(FUN, mc.cores=n_cores, ...)}
	}else{
		lapply_fun <- function(X, FUN, ...){lapply(X, FUN, ...)}
		mapply_fun <- function(FUN, ...){mapply(FUN, ...)}
	}

	# make sure p_adj is a real method
	if(! p_adj %in% p.adjust.methods){
		stop("p_adj not in p.adjust.methods.")
	}
	# make sure center is a valid argument
	if(! center %in% c("mean", "median", "mode")){
		stop(paste(center, "is not a valid argument for center. Use either \"mean\", \"median\", or \"mode\"."))
	}
	# make sure p_method is a valid argument
	if(! p_method %in% c("raw", "gamma_fit")){
		stop(paste(p_method, "is not a valid argument for p_method. Use \"raw\" or \"gamma_fit\"."))
	}

	# fit gamma distribution (incl reverse) to sim_raos
	gfits <- NULL
	if(p_method == "gamma_fit" || center == "mode"){
		msg("Fitting gamma dists to sim RAO...")
		tt <- proc.time()
		gfits <- lapply_fun(X=as.list(sim_raos), FUN=fit_gamma_fwd_rev)
		msg(paste0("...done (took ", fms(proc.time()-tt), ")"))
	}
	
	# calculate p-values
	msg("Calculating and adjusting P-values...")
	tt <- proc.time()
	Pval_perm <- mapply_fun(FUN=p_from_perms_or_gfit, emp=emp_raos, perm=sim_raos, tails=tails)
	if(p_method == "raw"){
		Pval <- Pval_perm
	}else if(p_method == "gamma_fit"){
		Pval <- mapply_fun(FUN=p_from_perms_or_gfit, emp=emp_raos, gfit=gfits, 
			fallback=Pval_perm, tails=tails)
	}
	# adjust p-values
	Pval <- p.adjust(Pval, method=p_adj)
	msg(paste0("...done (took ", fms(proc.time()-tt), ")"))

	# calculate means as fallback for mode
	msg("...Calculating sim RAO centers")
	spec_sim_means <- sapply(X=sim_raos, FUN=mean)
	if(center == "mean"){
		sim_rao_centers <- spec_sim_means
	}else if(center == "median"){
		sim_rao_centers <- sapply(X=sim_raos, FUN=median)
	}else if(center == "mode"){
		sim_rao_centers <- mapply_fun(FUN=mode_gamma_fwd_rev, fit=gfits, fallback=spec_sim_means)
	}

	# get denominator
	if(denom_type == "ses"){
		denom <- sapply(X=sim_raos, FUN=sd)
	}else if(denom_type == "raw"){
		denom <- rep(1, length(sim_raos))
	}else if(denom_type == "index"){
		# initialize denominator to sim_rao_centers, which is correct denominator for spec <= 0 otus
		denom <- sim_rao_centers
		msg("...Approximating max RAO values")
		# which species are overdispersed and need maxs calculated? This gives their indices
		species_inds_4_max <- which(emp_raos >  sim_rao_centers)
		# approximate max rao for each species indexed by species_inds_4_max
		max_raos <- apply(X=abunds_mat[,species_inds_4_max, drop=F], MARGIN=2, 
			FUN=rao_sort_max, D=env)
		# put newly calculated maxs where they belong (for spec > 0 otus)
		denom[species_inds_4_max] <- max_raos - sim_rao_centers[species_inds_4_max]
	}else if(denom_type == "sim_mode"){
		denom <- sim_rao_centers
	}else{
		stop("Invalid denom_type.")
	}
	msg(paste0("...done (took ", fms(proc.time()-tt), ")"))

	msg("Building output...")
	tt <- proc.time()
	out_specs <- ((emp_raos - sim_rao_centers)/ denom)
	out_specs <- round(out_specs, 4)
	
	# format output object
	output <- data.frame(Pval, Spec=out_specs)
	# if col names were input, use them for output too.
	if( !is.null(colnames(abunds_mat))){
		rownames(output) <- colnames(abunds_mat)
	}

	# if diagnostic output is desired, add columns to output
	if(diagnostic){
		specs_sim_mat_named <- do.call("rbind", sim_raos)
		colnames(specs_sim_mat_named) <- paste0("sim", 1:ncol(specs_sim_mat_named))

		output <- data.frame(output, 
			raw=emp_raos - sim_rao_centers,
			emp=emp_raos,
			sim_mean=sim_rao_centers,
			denom=denom,
			specs_sim_mat_named
		)
	}
	msg(paste0("...done (took ", fms(proc.time()-tt), ")"))

	return(output)
}