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
#'     \item{"index_full":}{
#'       d is the mean of simulated (permuted) RQE values for species that have stronger
#'       specificity than expected by chance, resulting in specificity values with range
#'       [-1, 0), with 0 as the null hypothesis. In this case, -1 indicates perfect
#'       specificity, where a species is associated with zero environmental variability.
#'       In the euclidean sense, this could be a species that is always found at the
#'       exact same elevation or the exact same pH. For species that have weaker specificity
#'       than expected by chance, d is x minus the center (see above) of simulated RQE 
#'       values, where x is the maximum possible dissimilarity observable given species
#'       weights. In "index_full", x is estimated using a genetic algorithm. This d has
#'       other useful properties: scale invariance to env/hosts_phylo, insensitivity to
#'       the number of samples, insensitivity to occupancy, and strong sensitivity to 
#'       specificity (default).
#'     }
#'     \item{"index_rough":}{
#'       Same as "index_full", but for species where specificity is weaker than expected by
#'       chance, a rough approximation is used to vastly reduce computational load. 
#'       This approximation does NOT give values where 1 is maximized generality. instead,
#'       the values will tightly correlate to the correct values, but the slope of that
#'       relationship will not be known a priori. In other words, if you don't care
#'       about species that are more general than expected by chance, this is fine. 
#'     }
#'     \item{"index_fast":}{
#'       Same as "index_full", but results as per "index_rough" are used together with
#'       a subset of results from "index_full" to create a model relating the two. Thus,
#'       speed is a compromise between the two approaches, with medium-high accuracy.
#'     }
#'   }
#' @param diagnostic logical. If true, changes output to include different parts of SES. 
#'   This includes Pval, SES, raw, denom, emp, and all sim values with column labels as
#'   simN where N is the number of sims (default: FALSE)
#' @param ga_params list. Parameters for genetic algorithm that maximizes RQE. Only used
#'   with denom_type="index_full/rough/fast". Default is the output of get_ga_defaults(). If different
#'   parameters are desired, start with output of get_ga_defaults and modify accordingly.
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
	denom_type="index_full", diagnostic=FALSE, ga_params=get_ga_defaults()){

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
	}else if(denom_type == "index_rough"){
		# initialize denominator to sim_rao_centers, which is correct denominator for spec <= 0 otus
		denom <- sim_rao_centers
		msg("...Approximating max RAO values (rough)")
		# which species are overdispersed and need maxs calculated? This gives their indices
		species_inds_4_max <- which(emp_raos > sim_rao_centers)

		# approximate max rao for each species indexed by species_inds_4_max
		# max_raos <- apply(X=abunds_mat[,species_inds_4_max, drop=F], MARGIN=2, 
		# 	FUN=rao_sort_max, D=env)
		max_raos <- unlist(lapply_fun(X=species_inds_4_max, FUN=function(j){
			rao_sort_max(abunds_mat[,j], D=env)}))

		# put newly calculated maxs where they belong (for spec > 0 otus)
		denom[species_inds_4_max] <- max_raos - sim_rao_centers[species_inds_4_max]
	}else if(denom_type %in% c("index", "index_full")){
		msg("...Approximating max RAO values (full)")
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
			max_raos <- rao_genetic_max(p=abunds_mat[species_inds_4_max[1]],
				D=env, term_cycles=ga_params$term_cycles, popsize=ga_params$popsize,
				keep=ga_params$keep, prc=ga_params$prc, maxiters=ga_params$maxiters)$best_rao
		}else if(length(species_inds_4_max) < (2 * n_cores)){
			max_raos <- unlist(lapply_fun(
				X=species_inds_4_max,
				FUN=function(j, par){
					rao_genetic_max(p=abunds_mat[,j], D=env, 
						term_cycles=par$term_cycles, popsize=par$popsize,
						keep=par$keep, prc=par$prc, maxiters=par$maxiters)$best_rao
				}, 
				par=ga_params
			))
		}else{
			seed_order <- rao_genetic_max(p=1:nrow(abunds_mat), D=env, 
				term_cycles=ga_params$term_cycles, popsize=ga_params$popsize,
				keep=ga_params$keep, prc=ga_params$prc, maxiters=ga_params$maxiters)$best_p
		
			max_raos <- unlist(lapply_fun(
				X=species_inds_4_max,
				FUN=function(j, par, ord){
					rao_genetic_max(p=sort(abunds_mat[,j])[ord], D=env, 
						term_cycles=par$term_cycles, popsize=par$popsize,
						keep=par$keep, prc=par$prc)$best_rao
				}, 
				par=ga_params, ord=seed_order
			))
		}

		# put newly calculated maxs where they belong (for spec > 0 otus)
		denom[species_inds_4_max] <- max_raos - sim_rao_centers[species_inds_4_max]
	}else if(denom_type == "index_fast"){
		msg("...Approximating max RAO values (fast)")
		# initialize denominator to sim_rao_centers, which is correct denominator for spec <= 0 otus
		denom <- sim_rao_centers
		# which species are overdispersed and need maxs calculated? This gives their indices
		species_inds_4_max <- which(emp_raos > sim_rao_centers)

		# determine how many values (=NR =number of representatives) to use for modeling
		NR <- 5
		if(NR < n_cores){NR <- n_cores}
		if(length(species_inds_4_max) <= NR){
			msg("......Approximation not neccesary (full instead)")

			# just use GA on each species instead
			# this also catches the case where NR > length(species_inds_4_max)
			max_raos <- unlist(lapply_fun(
				X=species_inds_4_max,
				FUN=function(j, par){
					rao_genetic_max(p=sort(abunds_mat[,j]), D=env, 
						term_cycles=par$term_cycles, popsize=par$popsize,
						keep=par$keep, prc=par$prc)$best_rao
				}, 
				par=ga_params
			))
			# put new values in
			denom[species_inds_4_max] <- max_raos - sim_rao_centers[species_inds_4_max]
		}else{
			msg("......1/3: Naive approximation")
			# calculate "naive" max value for all species in abunds_mat
			max_raos_n_log <- log(unlist(lapply_fun(X=1:ncol(abunds_mat), FUN=function(j){
				rao_sort_max(abunds_mat[,j], D=env)})))
			# choose NR representatives as even as possible
			ideals <- seq(from=min(max_raos_n_log), to=max(max_raos_n_log), length.out=NR)
			inds2get <- unique(sapply(X=ideals, FUN=function(i){which.min(abs(i-max_raos_n_log))}))
			if(length(inds2get) < NR){
				unusedinds <- (1:length(max_raos_n_log))[!(1:length(max_raos_n_log)) %in% inds2get]
				inds2get <- c(inds2get, sample(unusedinds, NR - length(inds2get)))
			}
			msg("......2/3: Gen. alg. subset")
			# use genetic algo on species whose indices are in inds2get
			ga_raos <- unlist(lapply_fun(
				X=inds2get,
				FUN=function(j, par){
					rao_genetic_max(p=sort(abunds_mat[,j]), D=env, 
						term_cycles=par$term_cycles, popsize=par$popsize,
						keep=par$keep, prc=par$prc)$best_rao
				}, 
				par=ga_params
			))
			msg("......3/3: Linear model and predict")
			# make linear model ga_raos ~ max_raos_all[inds2get]
			ga_mod <- lm(log(ga_raos) ~ max_raos_n_log[inds2get])
			msg(paste0(".........Rsquared=", summary.lm(ga_mod)$r.squared))
			preds <- exp( coef(ga_mod)[1] + (coef(ga_mod)[2] * max_raos_n_log) )

			# put new values in
			denom[species_inds_4_max] <- preds[species_inds_4_max] - sim_rao_centers[species_inds_4_max]
		}
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



