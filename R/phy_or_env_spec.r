#' phy_or_env_spec
#'
#' Calculates species' specificities to either a 1-dimensional variable (vector), 
#' 2-dimensional variable (matrix), or to a phylogeny. Transforms all variable
#' input types into a matrix D, and calculates specificity by comparing empirical
#' Rao's Quadratic Entropy to simulated RQE (same but with permuted abundances). 
#' By default (denom_type = "index"), an index is calculated from emp and sim
#' values such that Spec=0 indicates random assortment (null hypothesis), and more
#' negative values indicate stronger specificity.
#'
#' @author John L. Darcy
#' @references
#' \itemize{
#'   \item Poulin et al. (2011) Host specificity in phylogenetic and geographic
#'     space. Trends Parasitol 8:355-361. doi: 10.1016/j.pt.2011.05.003
#'   \item Rao CR (2010) Quadratic entropy and analysis of diversity. Sankhya 
#'     72:70-80. doi: 10.1007/s13171-010-0016-3
#'   \item Rao CR (1982) Diversity and dissimilarity measurements: A unified 
#'     approach. Theor Popul Biol 21:24-43.
#' }
#'
#' @param abunds_mat matrix or data frame of numeric values. Columns represent 
#'   species, rows are samples. For columns where the value is nonzero for two or
#'   fewer data points, specificity cannot be calculated, and NAs will be 
#'   returned. Negative values in abunds_mat are not allowed (REQUIRED).
#' @param env numeric vector, dist, or square matrix. Environmental variable 
#'   corresponding to abunds. For example, temperature, or geographic distance.
#'   Not required for computing phylogenetic specificity. If square matrix provided,
#'   note that only the lower triangle will be used (default: NULL). 
#' @param hosts character vector. Host identities corresponding to abunds. Only
#' required if calculating phylogenetic specificity (default: NULL).
#' @param hosts_phylo phylo object. Tree containing all unique hosts as tips. Only
#' required if calculating phylogenetic specificity (default: NULL).
#' @param n_sim integer. Number of simulations of abunds_mat to do under the null 
#'   hypothesis that host or environmental association is random. P-values will not
#'   be calculated if n_sim < 100 (default: 500).
#' @param p_adj string. Type of multiple hypothesis testing correction performed
#'   on P-values. Can take any valid method argument to p.adjust, including "none",
#'   "bonferroni", "holm", "fdr", and others (default: "fdr").
#' @param seed integer. Seed to use so that this is repeatable. Same seed will be
#'   used for each species in abunds_mat, so all species will experience the same 
#'   permutations. This can be disabled by setting seed=0, which will make permutation
#'   is both non deterministic (not repeatable) AND each species will experience
#'   different permutations (default: 1234557).
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
#' @param diagnostic logical. If true, changes output to include different parts of Spec. 
#'   This includes Pval, Spec, raw, denom, emp, and all sim values with column labels as
#'   simN where N is the number of sims (default: FALSE)
#' @param chunksize integer. If greater than zero, computation of sim RAO values will be
#'   done using chunked evaluation, which lowers memory use considerably for larger data
#'   sets. Can be disabled by setting to 0. Default value is 1000 species per chunk
#'   (default: 1000).
#' @param ga_params list. Parameters for genetic algorithm that maximizes RQE. Only used
#'   with denom_type="index". Default is the output of get_ga_defaults(). If different
#'   parameters are desired, start with output of get_ga_defaults and modify accordingly.
#'
#' @return data.frame where each row is an input species. First column is P-value
#'   ($Pval), second column is specificity ($Spec).
#'
#' @examples
#' # library(specificity)
#' # attach(endophyte)
#' # # only analyze species with occupancy >= 20
#' # m <- occ_threshold(prop_abund(otutable), 20)
#' # # create list to hold phy_or_env_spec outputs
#' # specs_list <- list()
#' # 
#' # # phylogenetic specificity using endophyte data set
#' # specs_list$host <- phy_or_env_spec(
#' #     abunds_mat=m,
#' #     hosts=metadata$PlantGenus, 
#' #     hosts_phylo=supertree,
#' #     n_sim=100, p_method="gamma_fit",
#' #     n_cores=4
#' # )
#' # 
#' # # environmental specificity using elevation from endophyte data set:
#' # specs_list$elev <- phy_or_env_spec(
#' #     abunds_mat=m,
#' #     env=metadata$Elevation,
#' #     n_sim=100, p_method="gamma_fit",
#' #     n_cores=4
#' # )
#' # 
#' # # geographic specificity using spatial data from endophyte data set:
#' # specs_list$geo <- phy_or_env_spec(
#' #     abunds_mat=m,
#' #     env=distcalc(metadata$Lat, metadata$Lon),
#' #     n_sim=100, p_method="gamma_fit",
#' #     n_cores=4
#' # )
#' # 
#' # plot_specs_violin(specs_list, cols=c("forestgreen", "red", "black"))
#'
#' @export
phy_or_env_spec <- function(abunds_mat, env=NULL, hosts=NULL, 
	hosts_phylo=NULL, n_sim=1000, p_adj="fdr", seed=1234567, 
	tails=1, n_cores=2, verbose=TRUE, p_method="raw", center="mean", 
	denom_type="index_full", diagnostic=F, chunksize=1000,
	ga_params=get_ga_defaults()){

	# debugging stuff:
	# library(specificity); attach(endophyte)
	# abunds_mat <- otutable[,colSums(otutable > 0) >= 10]
	# abunds_mat <- prop_abund(abunds_mat)
	# env = metadata$Elevation
	# n_sim=100; n_cores=10
	# hosts=NULL; hosts_phylo=NULL
	# p_adj="fdr"; seed=1234567; tails=1; verbose=TRUE
	# denom_type="index"
	# p_method <- "gamma_fit"
	# diagnostic <- TRUE
	# center <- "mean"
	# chunksize <- 50

    # warn if ncores > 1 and platform isn't  "unix" (windows can't do forked parallelism)
    if(n_cores > 1 && .Platform$OS.type != "unix"){
        warning("Windows is incompatible with n_cores > 1.")
    }
 
	# little message function, to clean up code a bit
	msg <- function(x){if(verbose){message(x)}}

	# function to pretty format seconds from proc.time()
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

	msg("Checking inputs...")
	tt <- proc.time()

	# replace NA/NaNs in abunds_mat with zeroes
	nNAs <- sum(is.na(abunds_mat))
	if(nNAs > 0){
		abunds_mat[is.na(abunds_mat)] <- 0
		warning(paste("Found", nNAs, "NA or NaN values in abunds_mat. Zero will be used in their place."))
	}

	# these checks are redundant with checks in calculate_spec_and_pval(), but they are done here
	# so that user knows there's an error before waiting for calculations to be done!
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

	# possible input types are: error, mat, dist, vec, phy
	data_type <- check_pes_inputs(abunds_mat, env, hosts, hosts_phylo, verbose)
	if(data_type == "error"){stop()}

	# turn it all into a dist, no matter what.
	if(data_type == "phy"){
		msg("...Converting tree to dist...")
		env <- tree2mat(tree=hosts_phylo, x=hosts, n_cores=n_cores)
	}else if(data_type == "vec"){
		msg("...Converting env vector to dist...")
		env <- dist(env)
	}else if(data_type == "mat"){
		env <- as.dist(env)
	}else if(data_type == "dist"){
		# do nothing, all good
	}else{
		stop("unknown error")
	}

	# check for integers in abunds_mat. integers can cause issues because
	# .Machine$integer.max is pretty small, and is easily outpaced by rao computation.
	# no, apply(abunds_mat, FUN=is.integer, MARGIN=2) does NOT work on dfs, don't know why.
	# test <- data.frame(a=1:2, b=c(0.1, 0.2))
	# apply(test, FUN=is.integer, MARGIN=2)                         # WRONG
	# lapply(test, FUN=is.integer )                                 # CORRECT (will work on df but not matrix)
	# sapply(X=1:ncol(test), FUN=function(j){is.integer(test[,j])}) # CORRECT (works on both)
	abunds_mat_intcols <- sapply(X=1:ncol(abunds_mat), FUN=function(j){is.integer(abunds_mat[j])})
	if(any(abunds_mat_intcols)){
		msg("...Changing abunds_mat ints to doubles...")
		abunds_mat <- apply(X=abunds_mat, MARGIN=2, FUN=as.double)
	}
	# not necessary to do the same for env since double * int = double.
	
	# below is old code used to check if this function could produce
	# integers that were too large. new code (above) converts integer
	# values to doubles, so it should be fine. 
		# check if max rao is an int an if it's too large
		# col2check <- which.max(colSums(abunds_mat))
		# max2check <- rao_sort_max(abunds_mat[,col2check], D=env)
		# if(max2check > .Machine$integer.max){
		# 	stop("Maximum possible RAO value is greater than max integer value.
		# 		Use prop_abund() on abunds_mat, and/or divide env by a scalar 
		# 		to remedy this.")
		# }
	msg(paste0("...done (took ", fms(proc.time()-tt), ")"))


	# make wrapper for lapply/mclapply so with cores=1 it just uses lapply
	# also for mapply
	if(n_cores > 1){
		lapply_fun <- function(X, FUN, ...){parallel::mclapply(X, FUN, mc.cores=n_cores, ...)}
		mapply_fun <- function(FUN, ...){parallel::mcmapply(FUN, mc.cores=n_cores, ...)}
	}else{
		lapply_fun <- function(X, FUN, ...){lapply(X, FUN, ...)}
		mapply_fun <- function(FUN, ...){mapply(FUN, ...)}
	}

	# calculate empirical rao
	msg(paste("Calculating emp RAO..."))
	tt <- proc.time()
	# emp_raos <- apply(X=abunds_mat, MARGIN=2, FUN=rao1sp, D=env, perm=FALSE)
	emp_raos <- unlist(lapply_fun(X=1:ncol(abunds_mat), FUN=function(x){
		rao1sp(abunds_mat[,x], D=env, perm=FALSE)
	}))
	msg(paste0("...done (took ", fms(proc.time()-tt), ")"))

	# calculate sim rao
	msg(paste("Calculating sim RAO..."))
	tt <- proc.time()
	# handle chunksize - don't chunk if chunksize is too large or is 0
	if(chunksize > ncol(abunds_mat) || chunksize<=0){
		sim_raos <- lapply_fun(X=1:ncol(abunds_mat), FUN=function(p){
			raoperms(p=abunds_mat[,p], D=env, n_sim=n_sim, seed=seed)})
	}else{
		# list of chunks by column index, i.e. first item in list is all indices of
		# features (columns) that are to be used in first chunk
		msg(paste("...Chunking..."))
		nchunks <- ceiling(ncol(abunds_mat) / chunksize)
		chunks <- split(1:ncol(abunds_mat), rep(1:nchunks, each=chunksize,
			length.out=ncol(abunds_mat)))

		msg(paste("...Calculating..."))
		sim_raos <- list()
		if(verbose){pb <- txtProgressBar(min=0, max=length(chunks), style=3)}
		for(ch in 1:length(chunks)){
			sim_raos[[ch]] <- lapply_fun(X=chunks[[ch]], FUN=function(p){
				raoperms(p=abunds_mat[,p], D=env, n_sim=n_sim, seed=seed)})
			if(verbose){setTxtProgressBar(pb, ch)}
		}
		if(verbose){close(pb)}

		msg(paste("...Unchunking..."))
		# de-chunk
		sim_raos <- unlist(sim_raos, recursive=FALSE)
	}

	msg(paste0("...done (took ", fms(proc.time()-tt), ")"))

	# calculate output specificity
	output <- calculate_spec_and_pval(emp_raos=emp_raos, sim_raos=sim_raos, 
		abunds_mat=abunds_mat, env=env, p_adj=p_adj, tails=tails, n_cores=n_cores,
		verbose=verbose, p_method=p_method, center=center, denom_type=denom_type,
		diagnostic=diagnostic, ga_params=ga_params)
	return(output)
}

