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
#' Poulin et al. (2011) Host specificity in phylogenetic and geographic space. 
#'   Trends Parasitol 8:355-361. doi: 10.1016/j.pt.2011.05.003
#' Rao CR (2010) Quadratic entropy and analysis of diversity. Sankhya 72:70-80.
#'   doi: 10.1007/s13171-010-0016-3
#' Rao CR (1982) Diversity and dissimilarity measurements: A unified approach.
#'   Theor Popul Biol 21:24-43.
#'
#' @param abunds_mat matrix or data frame of numeric values. Columns represent 
#'   species, rows are samples. For columns where the value is nonzero for two or
#'   fewer data points, environmental SES cannot be calculated, and NAs will be 
#'   returned. Negative values in abunds_mat are not allowed (REQUIRED).
#' @param env numeric vector, dist, or square matrix. Environmental variable 
#'   corresponding to abunds. For example, temperature, or geographic distance.
#'   Not required for computing phylogenetic specificity (DEFAULT: NULL).
#' @param hosts character vector. Host identities corresponding to abunds. Only
#' required if calculating SES for phylogenetic specificity (DEFAULT: NULL).
#' @param hosts_phylo phylo object. Tree containing all unique hosts as tips. Only
#' required if calculating SES for phylogenetic specificity (DEFAULT: NULL).
#' @param n_sim integer. Number of simulations of abunds_mat to do under the null 
#'   hypothesis that host or environmental association is random. P-values will not
#'   be calculated if n_sim < 100 (DEFAULT: 500).
#' @param sim_fun function. A function f where f(abunds_mat) returns a matrix
#'   object with the same number of rows and columns as abunds_mat. Default is
#'   f=function(m){ m[sample(1:nrow(m)),]}, which just permutes the order of 
#'   rows in abunds_mat. Users may wish to use a null model that is able to
#'   preserve row and column totals such as the function permatswap() from the
#'   vegan package or the function vaznull() from the bipartite package. Either
#'   of these can be easily adapted to return only a single matrix (see examples). 
#'   However, neither can accomodate non-integer matrices. 
#' @param p_adj string. Type of multiple hypothesis testing correction performed
#'   on P-values. Can take any valid method argument to p.adjust, including "none",
#'   "bonferroni", "holm", "fdr", and others (DEFAULT: "fdr").
#' @param seed integer. Seed to use so that this is repeatable (DEFAULT: 1234557).
#' @param tails integer. 1 = 1-tailed, test for specificity only. 2 = 2-tailed.
#'   3 = 1-tailed, test for cosmopolitanism only. 0 = no test, P=1.0 (DEFAULT: 1).
#' @param n_cores integer. Number of CPU cores to use for parallel operations. If
#'   set to 1, lapply will be used instead of mclapply (DEFAULT: 2).
#' @param verbose logical. Should status messages be displayed? (DEFAULT: TRUE).
#' @param p_method string. method argument to pval_from_perms (DEFAULT: "raw").
#' @param denom_type string. Type of denominator (d) to use (DEFAULT: "index"). Note
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
#'       than expected by chance, d is x minus the mean of simulated RQE values, where x is
#'       the maximum possible dissimilarity observable given species weights. This d has 
#'       other useful properties: scale invariance to env/hosts_phylo, insensitivity to 
#'       the number of samples, insensitivity to occupancy, and strong sensitivity to 
#'       specificity (DEFAULT).
#'     }
#'   }
#' @param diagnostic logical. If true, changes output to include different parts of SES. 
#'   This includes Pval, SES, raw, denom, emp, and all sim values with column labels as
#'   simN where N is the number of sims (DEFAULT: FALSE)
#'
#'
#' @return data.frame where each row is an input species. First column is P-value
#'   ($Pval), second column is specificity ($Spec).
#'
#' @examples
#' # phylogenetic specificity using endophyte data set
#' attach(endophyte)
#' # only analyze species with occupancy >= 20
#' m <- occ_threshold(prop_abund(zotutable), 20)
#' ses_host <- phy_or_env_spec(
#'     abunds_mat=m,
#'     hosts=metadata$PlantGenus, 
#'     hosts_phylo=supertree,
#'     n_cores=12
#' )
#'
#' # using vazquez null model from bipartite package as an alternate permutation:
#' # note that the "creating permuted matrices" step will be slow.
#' library(bipartite)
#' ses_host_vaz <- phy_or_env_spec(
#'     abunds_mat=m,
#'     hosts=metadata$PlantGenus, 
#'     hosts_phylo=supertree,
#'     n_cores=12,
#'     sim_fun=function(m){bipartite::vaznull(1, m)[[1]]},
#' )
#'
#' # compare naive permutation vs. vazquez:
#' plot(ses_host$Spec, ses_host_vaz$Spec, ylab="bipartite::vaznull", xlab="naive")
#' abline(h=0);abline(v=0)
#' hist(ses_host_vaz$Spec)
#' hist(ses_host$Spec)
#'
#' # environmental specificity using elevation from endophyte data set:
#' ses_elev <- phy_or_env_spec(
#'     abunds_mat=m,
#'     env=metadata$Elevation,
#'     n_cores=12
#' )
#' 
#' # geographic specificity using spatial data from endophyte data set:
#' ses_geo <- phy_or_env_spec(
#'     abunds_mat=m,
#'     env=distcalc(metadata$Lat, metadata$Lon),
#'     n_cores=12
#' )
#'
#' @export
phy_or_env_spec <- function(abunds_mat, env=NULL, hosts=NULL, 
	hosts_phylo=NULL, n_sim=1000,  sim_fun=function(m){ m[sample(1:nrow(m)), ] }, 
	p_adj="fdr", seed=1234567, tails=1, n_cores=2, verbose=TRUE,
	p_method="raw", denom_type="index", diagnostic=F){

	# debugging stuff:
	# library(specificity); attach(endophyte)
	# abunds_mat <- otutable[,colSums(otutable > 0) >= 10]
	# env = metadata$Elevation
	# n_sim=100; n_cores=10
	# hosts=NULL; hosts_phylo=NULL
	# sim_fun=function(m){ m[sample(1:nrow(m)), ] }
	# p_adj="fdr"; seed=1234567; tails=1; verbose=TRUE
	# denom_type="index"
	# p_method <- "gamma_fit"
	# diagnostic <- TRUE


	require(ape)
	require(geiger)

	# little message function, to clean up code a bit
	msg <- function(x){if(verbose){message(x)}}

	msg("Checking inputs.")
	# make sure p_adj is a real method
	if(! p_adj %in% p.adjust.methods){
		stop("p_adj not in p.adjust.methods.")
	}
	# possible input types are: error, mat, dist, vec, phy
	data_type <- check_pes_inputs(abunds_mat, env, hosts, hosts_phylo, verbose)
	if(data_type == "error"){stop()}

	# turn it all into a dist, no matter what.
	if(data_type == "phy"){
		msg("Converting tree to dist.")
		env <- tree2mat(tree=hosts_phylo, x=hosts, n_cores=n_cores, delim=";")
	}else if(data_type == "vec"){
		msg("Converting env vector to dist.")
		env <- dist(env)
	}else if(data_type == "mat"){
		env <- as.dist(env)
	}else if(data_type == "dist"){
		# do nothing, all good
	}else{
		stop("unknown error")
	}

	# make wrapper for lapply/mclapply so with cores=1 it just uses lapply
	if(n_cores > 1){
		require("parallel")
		lapply_fun <- function(X, FUN, ...){mclapply(X, FUN, mc.cores=n_cores, ...)}
	}else{
		lapply_fun <- function(X, FUN, ...){lapply(X, FUN, ...)}
	}

	# generate n_sim daughter seeds for generation of each permuted matrix
	seeds <- daughter_seeds(n=n_sim, s=seed)

	# function to generate matrix given seed:
	perm_mat <- function(s){
		set.seed(s)
		return((sim_fun(abunds_mat)))
	}
	# create n_sim permuted matrices, but convert them into a list of column vectors.
	msg(paste("Creating", n_sim, "permuted matrices."))
	perm_cols <- do.call("cbind", lapply_fun(X=seeds, FUN=perm_mat))

	msg(paste("Calculating sim RAO."))
	# split perm cols into a list of length n_cores
	split_matrix <- function(m, g){
		sp <- sort(rep(1:g, length.out=ncol(m)))
		return(lapply(X=unique(sp), FUN=function(x){m[,sp==x, drop=FALSE]}))
	}
	sim_submats <- split_matrix(as.matrix(perm_cols), n_cores)
	# calculate specs for perm_cols in parallel
	specs_sim <- unlist(lapply_fun(X=sim_submats, FUN=spec_core, D=env))
	specs_sim_mat <- matrix(specs_sim, ncol=ncol(abunds_mat), byrow = TRUE )

	# use spec_fun on empirical data
	msg("Calculating emp RAO.")
	emp_submats <- split_matrix(as.matrix(abunds_mat), n_cores)
	specs_emp <- unlist(lapply_fun(X=emp_submats, FUN=spec_core, D=env))

	# calculate p-values
	msg("Calculating P-values.")
	Pval <- rep(-1, length(specs_emp))
	for(i in 1:length(Pval)){
		Pval[i] <- pval_from_perms(emp=specs_emp[i], perm=specs_sim_mat[,i], 
			tails=tails, method=p_method, rounding=4)
	}
	# adjust p-values
	Pval <- p.adjust(Pval, method=p_adj)

	# calculate output specificity
	msg("Calculating specificities.")
	# note these aren't necessarily means! Trying to get mode instead, to avoid
	# issue of skewed distributions. 
	spec_sim_means <- unlist(lapply_fun(X=as.data.frame(specs_sim_mat), FUN=gamma_mode))
	# spec_sim_means <- apply(X=specs_sim_mat, MAR=2, FUN=mean)

	# get denominator
	if(denom_type == "ses"){
		denom <- apply(X=specs_sim_mat, MAR=2, FUN=sd)
	}else if(denom_type == "raw"){
		denom <- rep(1, ncol(specs_sim_mat))
	}else if(denom_type == "index"){
		# initialize denominator to spec_sim_means, which is correct denominator for spec <= 0 otus
		denom <- spec_sim_means
		msg("Approximating max RAO values.")
		# which otus are overdispersed and need maxs calculated? This gives their indices
		otu_inds_4_max <- which(specs_emp >  spec_sim_means)
		# use above indices to make a list of corresponding column vectors from abunds mat
		emp_col_list <- lapply(X=otu_inds_4_max, FUN=function(x){abunds_mat[,x]})
		# calculate max rao values for those cols in parallel
		maxraos <- unlist(lapply_fun(X=emp_col_list, FUN=rao_sort_max, D=env))
		# put newly calculated maxs where they belong (for spec > 0 otus)
		denom[specs_emp > spec_sim_means] <- maxraos
		# for(j in otu_inds_4_max){
		# 	denom[j] <- rao_sort_max(w=abunds_mat[,j], D=env)
		# }
	}else if(denom_type == "sim_mode"){
		denom <- spec_sim_means
	}else{
		stop("Invalid denom_type.")
	}

	out_specs <- ((specs_emp - spec_sim_means)/ denom)
	out_specs <- round(out_specs, 4)
	
	# format output object
	output <- data.frame(Pval, Spec=out_specs)
	# if col names were input, use them for output too.
	if( !is.null(colnames(abunds_mat))){
		rownames(output) <- colnames(abunds_mat)
	}

	# if diagnostic output is desired, add columns to output
	if(diagnostic){
		specs_sim_mat_named <- t(specs_sim_mat)
		colnames(specs_sim_mat_named) <- paste0("sim", 1:ncol(specs_sim_mat_named))

		output <- data.frame(output, 
			raw=specs_emp - spec_sim_means,
			emp=specs_emp,
			sim_mean=spec_sim_means,
			denom=denom,
			specs_sim_mat_named
		)
	}

	msg("Done.")
	return(output)
}
