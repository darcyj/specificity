#' calculate_spec_and_pval
#'
#' This function is called by phy_or_env_spec(). It is made available as a
#' standalone function in the (rare) case a user wishes to calculate Spec
#' using their own null model. calculate_spec_and_pval() takes empirical rao
#' values and sim rao values (from a null model) and calculates Spec and
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
#'   set to 1, lapply will be used instead of mclapply. A warning will be shown if
#'   n_cores > 1 on Windows, which does not support forked parallelism (default: 2).
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
#'       d is the center of simulated (permuted) RQE values for species that have stronger
#'       specificity than expected by chance, resulting in specificity values with range
#'       [-1, 0), with 0 as the null hypothesis. In this case, -1 indicates perfect
#'       specificity, where a species is associated with zero environmental variability.
#'       In the euclidean sense, this could be a species that is always found at the
#'       exact same elevation or the exact same pH. For species that have weaker specificity
#'       than expected by chance, d is x minus the center (see above) of simulated RQE 
#'       values, where x is the maximum possible dissimilarity observable given species
#'       weights. x is estimated using a genetic algorithm. This d has
#'       other useful properties: scale invariance to env/hosts_phylo, insensitivity to
#'       the number of samples, insensitivity to occupancy, and strong sensitivity to 
#'       specificity (default).
#'     }
#'     \item{"sim_center":}{
#'       d is always the center of simulated (permuted) RQE values. For species that have
#'       stronger specificity than expected by chance, this will return the same Spec 
#'       values as "index". For species with weaker specificity than expected by chance,
#'       instead of values that range between 0 and 1, they will range between 0 and Inf.
#'       This is much faster than "index" because the genetic algorithm is not used. So 
#'       if species with weaker specificity than expected by chance are not interesting
#'       to you, this may be a good option.
#'     }
#'   }
#' @param diagnostic logical. If true, changes output to include different parts of SES. 
#'   This includes Pval, SES, raw, denom, emp, and all sim values with column labels as
#'   simN where N is the number of sims (default: FALSE)
#' @param ga_params list. Parameters for genetic algorithm that maximizes RQE. Only used
#'   with denom_type="index". Default is the output of get_ga_defaults(). 
#'   If different parameters are desired, start with output of get_ga_defaults and modify
#'   accordingly.
#'
#'
#' @return data.frame where each row is an input species. First column is P-value
#'   ($Pval), second column is specificity ($Spec).
#'
#' @examples
#' # # calculating regular old elevational specificity the hard way
#' # attach(endophyte)
#' # library(parallel)
#' # otutable <- occ_threshold(prop_abund(otutable), 10)
#' # env <- dist(metadata$Elevation)
#' # emp_raos <- apply(X=otutable, MARGIN=2, FUN=rao1sp, 
#' #   D=env, perm=F, seed=12345)
#' # sim_raos <- mclapply(X=as.data.frame(otutable), FUN=function(p){
#' #   replicate(200, rao1sp(p, D=env, perm=TRUE, seed=0))}, mc.cores=20)
#' # calculate_spec_and_pval(emp_raos, sim_raos, otutable, env, 
#' #   n_cores=20)
#' 
#' 
#' @export
calculate_spec_and_pval <- function(emp_raos, sim_raos, abunds_mat, env,
	p_adj="fdr", tails=1, n_cores=2, verbose=TRUE, p_method="raw", center="mean", 
	denom_type="index", diagnostic=FALSE, ga_params=get_ga_defaults()){

	# warn if ncores > 1 and platform isn't  "unix" (windows can't do forked parallelism)
	if(n_cores > 1 && .Platform$OS.type != "unix"){
		warning("Windows is incompatible with n_cores > 1.")
	}

	# check if env is a dist
	if( ! inherits(env, "dist")){
		stop("env must be of class dist")
	}

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

	# calculate centers
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
	}else if(denom_type %in% c("index", "index_full")){
		msg("...Approximating max RAO values")
		# initialize denominator to sim_rao_centers, which is correct denominator for spec <= 0 otus
		denom <- sim_rao_centers
		# which species are overdispersed and need maxs calculated? This gives their indices
		species_inds_4_max <- which(emp_raos > sim_rao_centers)

		# if there are twice as many species_inds_4_max than there are available cores,
		# do order optimization first. otherwise, just do it in parallel, unless there's
		# just one, then don't.
		if(length(species_inds_4_max) <= 0){
			max_raos <- NULL
		}else if(length(species_inds_4_max) == 1){
			# use GA to get maxrao for that one species
			max_raos <- rao_genetic_max(
				p=abunds_mat[species_inds_4_max[1]],
				D=env,
				swap_freq=ga_params$swap_freq,
				term_cycles=ga_params$term_cycles, 
				popsize_swap=ga_params$popsize_swap,
				popsize_perm=ga_params$popsize_perm,
				keep=ga_params$keep,
				cross=ga_params$cross,
				prc=ga_params$prc, 
				maxiters=ga_params$maxiters
			)$best_rao
		}else{
			max_raos <- unlist(lapply_fun(
				X=species_inds_4_max,
				FUN=function(j, par){
					rao_genetic_max(
						p=abunds_mat[,j],
						D=env,
						swap_freq=par$swap_freq,
						term_cycles=par$term_cycles, 
						popsize_swap=par$popsize_swap,
						popsize_perm=par$popsize_perm,
						keep=par$keep, 
						cross=par$cross,
						prc=par$prc, 
						maxiters=par$maxiters
					)$best_rao
				}, 
				par=ga_params
			))
		}
		# put newly calculated maxs where they belong (for spec > 0 otus)
		denom[species_inds_4_max] <- max_raos - sim_rao_centers[species_inds_4_max]
	}else if(denom_type %in% c("sim_center", "sim_mode")){
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



