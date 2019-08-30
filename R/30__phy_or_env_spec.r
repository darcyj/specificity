#' phy_or_env_spec
#'
#' Calculates standardized effect size (SES) of species' specificity to either
#' a 1-dimensional variable (vector), 2-dimensional variable (matrix), or to
#' a phylogeny. Default operation is to transform all variable input types into a 
#' matrix M, and calculate specificity as SES of the empirical sum of M
#' weighted by weights matrix W. W is a matrix of pair-wise products of species
#' species abundances (Rao's quadratic entropy). Optionally, host phylogenetic
#' specificity can be calculated per Poulin et al. (2011) using weighted 
#' phylogenetic entropy, and specificity to 1-dimensional data can be calculated
#' using a weighted variance approach. These approaches are included mainly for
#' comparative methodological purposes. 
#'
#' @author John L. Darcy
#' @references 
#' Poulin et al. (2011) Host specificity in phylogenetic and geographic
#'   space. Trends Parasitol 8:355-361. doi: 10.1016/j.pt.2011.05.003
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
#' @param lowmem logical. Should this function be run in a way that saves memory?
#'   If TRUE, will run FAR slower but use almost no ram. If FALSE, lots of ram is
#'   used but it will run quickly. For example, on an 800x1600 matrix, with nperm=
#'   1000, n_cores=12 and lowmem=F, will use roughly 90 GB ram but finish in one
#'   minute. With the same settings but lowmem=T, will use less than 1 GB ram but 
#'   take over an hour. A progress bar is shown if TRUE (DEFAULT: FALSE).
#' @param vec_wvar logical. For vector input, should old weighted-variance approach
#'   be used instead of transforming the vector into a euclidean distance matrix?
#'   (DEFAULT: FALSE).
#' @param phy_ent logical. For phylogenetic input, should SES of phylogenetic
#'   entropy be used instead of trandforming the phylogeny into patristic distances
#'   using tree2mat? (DEFAULT: FALSE).
#' @param denom_type string. Type of denominator (d) to use for SES. SES is the 
#'   difference between empirical and mean simulated specificity, divided by d
#'   (DEFAULT: global_unif).
#'   \describe{
#'     \item{species_sim:}{
#'       d for species s is calculated as the standard deviation of specificities
#'       calculated from permuted abundances of s. This makes SES of specificity
#'       counterintuitively sensitive to occupancy, where species with high occupancy
#'       have more extreme SES of specificity than rare species, due to their more 
#'       deterministic sim specificities. Not suggested.
#'     }
#'     \item{global_sim:}{
#'       d is same for all species. Calculated as standard deviation of all
#'       specificities calculated from permuted species. SES is much less sensitive
#'       to occupancy, and intuitively, species with greater occupancy have less
#'       extreme SES of specificity than rare species. However, using this d makes
#'       results from different abundance matrices NOT comparable, because d is a
#'       function of the distribution of species distributions unique to abunds_mat.
#'       Results for different variables and the same abundance matrix are still 
#'       comparable, though (obviously using same denom_type and other parameters).
#'     }
#'     \item{global_unif:}{
#'       d is same for all species. Calculated as variability in specificity under
#'       random uniform distribution (beta 1,1) of species abunds. This d is 
#'       comparable between different abundance matrices and between different
#'       variables. SES results using this denom_type are NOT comparable with results
#'       from any other. Due to its stability, this is the default option.
#'     }
#'     \item{raw:}{
#'       d is 1 for all species, so SES is not actually standardized (although the
#'       column label in output data will still be "SES"). This option only exists
#'       for diagnostic purposes, since units of this metric are incomprehensible.
#'     }
#'   }
#'
#'
#' @return data.frame where each row is an input species. First column is P-value
#'   ($Pval), second column is SES ($SES).
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
#' plot(ses_host$SES, ses_host_vaz$SES, ylab="bipartite::vaznull", xlab="naive")
#' abline(h=0);abline(v=0)
#' hist(ses_host_vaz$SES)
#' hist(ses_host$SES)
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
	p_adj="fdr", seed=1234567, tails=1, n_cores=2, verbose=TRUE, lowmem=FALSE, 
	vec_wvar=FALSE, phy_ent=FALSE, denom_type="global_unif"){

	# testing crap:
	# abunds_mat = data.frame(z_flat, z_elev_lo, z_elev_med, 
	#   z_elev_hi, z_elev_lo_all, z_elev_med_all, z_elev_hi_all)
	# env = metadata$Elevation
	# n_sim=1000
	# n_cores=10
	# hosts=NULL
	# hosts_phylo=NULL; n_sim=1000
	# sim_fun=function(m){ m[sample(1:nrow(m)), ] }
	# p_adj="fdr"; seed=1234567; tails=1; verbose=TRUE; lowmem=FALSE
	# vec_wvar=FALSE; phy_ent=FALSE; denom_type="global_unif"
	#   

	require(ape)
	require(geiger)

	# little message function, to clean up code a bit
	msg <- function(x){if(verbose){message(x)}}

	msg("Checking inputs.")
	# make sure p_adj is a real method
	if(! p_adj %in% p.adjust.methods){
		stop("p_adj not in p.adjust.methods.")
	}
	# decide inputs type
	env_prsnt <- !is.null(env)
	hos_prsnt <- !is.null(hosts)
	phy_prsnt <- !is.null(hosts_phylo)
	if(env_prsnt && (hos_prsnt || phy_prsnt)){
		stop("Too many inputs. Pick either env only, or hosts+hosts_phylo.")
	}else if(env_prsnt && (is.matrix(env) || is.data.frame(env) || class(env) == "dist")){
		# env present and it's 2-dimensional
		# make sure env is numeric
		if(!is.numeric(env)){
			stop("env not numeric.")
		}
		# if it's a data.frame or matrix, make sure it's square
		# and turn it into a dist
		if(is.matrix(env) || is.data.frame(env)){
			if(nrow(env) != ncol(env)){
				stop("env is 2-dimensional, but not square.")
			}else{
				env <- as.dist(env)
			}
		}
		data_type <- "dist"
	}else if(env_prsnt && is.vector(env)){
		# env present and it's a vector
		# make sure env is numeric
		if(!is.numeric(env)){
			stop("env not numeric.")
		}
		data_type <- "vec"
	}else if(hos_prsnt || phy_prsnt){
		if(sum(hos_prsnt, phy_prsnt) < 2){
			stop("both hosts AND hosts_phylo must be input for phylogenetic specificity.")
		}
		data_type <- "phy"
	}else{
		message(paste("Is env a vector?", is.vector(env)))
		message(paste("Is env a matrix?", is.matrix(env)))
		message(paste("Is env a data.frame?", is.data.frame(env)))
		message(paste("Is env a dist?", (class(env) == "dist")))
		message(paste("is hosts a vector?", is.vector(hosts)))
		message(paste("is hosts a factor?", is.factor(hosts)))
		message(paste("is hosts_phylo a tree?", is.phylo(hosts_phylo)))
		stop("Type error with env, hosts, or hosts_phylo. See diagnostics above.")
	}

	# check to make sure inputs are compatible
	predict_distlength <- function(x){((x^2)-x) /2}
	if(data_type == "phy"){
		if(! all(hosts_phylo$tip.label %in% hosts)){
			stop("Some hosts are missing from hosts_phylo.")
		}else if(length(hosts) != nrow(abunds_mat)){
			stop("hosts and abunds_mat have incompatible dimensions.")
		}else if(any(grepl(";", x=hosts_phylo$tip.label))){
			stop("hosts cannot contain semicolons.")
		}
	}else if (data_type == "vec"){
		if(! is.numeric(env)){
			stop("env is not numeric.")
		}else if(length(env) != nrow(abunds_mat)){
			stop("env and abunds_mat have incompatible dimensions.")
		}
	}else if(data_type == "dist"){
		if(! is.numeric(env)){
			stop("env is not numeric.")
		}else if(length(env) != predict_distlength(nrow(abunds_mat))){
			stop("env and abunds_mat have incompatible dimensions.")
		}
	}

	# turn vector or phy input into dist and set data_type to dist unless 
	# requested otherwise
	if(data_type == "phy" && phy_ent == FALSE){
		msg("Converting tree to dist.")
		env <- tree2mat(tree=hosts_phylo, x=hosts, n_cores=n_cores, delim=";")
		data_type <- "dist"
	}else if(data_type == "vec" && vec_wvar == FALSE){
		msg("Converting env vector to dist.")
		env <- dist(env)
		data_type <- "dist"
	}


	# define specificity function
	ns <- NULL # nested set for phy stuff
	if(data_type == "phy"){
		# if it's still phy, we're doing phylogenetic entropy. 
		# go ahead and prune tree and make nsested set, too.
		hosts_phylo <- keep.tip(hosts_phylo, hosts)
		ns <- make_nested_set(hosts_phylo, n_cores)
		# objects that spec_fun needs (other than abundance vector ab) 
		# are put into a list called extra_inputs
		extra_inputs <- list(hosts=hosts, hosts_phylo=hosts_phylo, ns=ns)
		spec_fun <- function(ab, extras){
			return(wpd(s=extras$hosts, s_phylo=extras$hosts_phylo, w=ab, nested_set=extras$ns, metric="Hp"))
		}
	}else if(data_type == "vec"){
		# if data_type is still "vec", we're doing weighted variance approach. 
		# need to define weighted variance function:
		w_var <- function(x, w){
			if(sum(w > 0) < 2){
				return(NA)
			}else{
				return(sum(w * (x - weighted.mean(x,w))^2) / (sum(w) - (sum(w^2 / sum(w)))) )
			}
		}
		# objects that spec_fun needs (other than abundance vector ab) 
		# are put into a list called extra_inputs
		extra_inputs <- list(env=env)
		spec_fun <- function(ab, extras){
			return(w_var(x=extras$env, w=ab))
		}
	}else if(data_type == "dist"){
		# data_type is "dist", the default. 
		# objects that spec_fun needs (other than abundance vector ab) 
		# are put into a list called extra_inputs
		# note that "env" on next line may be a phylo distmat or a euclidean dismat from a vector!!
		extra_inputs <- list(env=env)
		spec_fun <- function(ab, extras){
			return(rao_quad_ent(d=extras$env, w=ab, raw=FALSE))
		}
	}else{
		stop("mystery error - data type undefined.")
	}

	# make wrapper for lapply/mclapply so with cores=1 it just uses lapply
	if(n_cores > 1){
		require("parallel")
		lapply_fun <- function(X, FUN, ...){mclapply(X, FUN, mc.cores=n_cores, ...)}
	}else{
		lapply_fun <- function(X, FUN, ...){lapply(X, FUN, ...)}
	}

	# generate n_sim daughter seeds for generation of each permuted matrix
	msg("Generating daughter seeds.")
	seeds <- daughter_seeds(n=n_sim, s=seed)

	# function to generate matrix given seed:
	perm_mat <- function(s){
		set.seed(s)
		return((sim_fun(abunds_mat)))
	}
	if(lowmem == FALSE){
		# create n_sim permuted matrices, but convert them into a list of column vectors.
		msg(paste("Creating", n_sim, "permuted matrices."))
		perm_cols <- do.call("cbind", lapply_fun(X=seeds, FUN=perm_mat))
		# turn list of matrices into a list of column vectors
		perm_cols <- lapply(seq_len(ncol(perm_cols)), function(i){perm_cols[,i]})

		# apply spec_fun to each column vector in perm_cols
		msg("Calculating specificities for permuted matrices.")
		# bs = batch size = number of seqs to process at a time
		bs <- n_cores * 8
		batches <- rep(1:bs, each=ceiling(length(perm_cols) / bs))[1:length(perm_cols)]
		specs_sim <- rep(0, length(perm_cols))
		if(verbose){pb <- txtProgressBar(min=0, max=bs, style=3)}
		for(b in 1:bs){
			specs_sim[batches==b] <- unlist(lapply_fun(X=perm_cols[batches==b], 
				FUN=spec_fun, extras=extra_inputs))
			if(verbose){setTxtProgressBar(pb, b)}
		}
		if(verbose){close(pb)} # to kill of progress bar
		#specs_sim <- unlist(lapply_fun(X=perm_cols, FUN=spec_fun, extras=extra_inputs))
		# turn it into a matrix, where each row is a simulation, and each column is a species.
		specs_sim_mat <- matrix(specs_sim, ncol=ncol(abunds_mat), byrow = TRUE )
	}else{
		msg("Calculating specificities for permuted matrices (lomem).")
		specs_sim_mat <- matrix(0, nrow=n_sim, ncol=ncol(abunds_mat))
		if(verbose){pb <- txtProgressBar(min=0, max=n_sim, style=3)}
		nc <- ncol(abunds_mat)
		for(i in 1:n_sim){
			# permute, adjust abunds_mat, and split it into list of columns.
			pm <- perm_mat(seeds[i])
			perm_cols <- lapply(seq_len(nc), function(i){pm[,i]})
			specs_sim_mat[i,] <- unlist(lapply_fun(X=perm_cols, FUN=spec_fun, extras=extra_inputs))
			if(verbose){setTxtProgressBar(pb, i)}
		}
		if(verbose){close(pb)}# to kill of progress bar
	}

	# use spec_fun on empirical data
	msg("Calculating empirical specificities.")
	emp_cols <- lapply(seq_len(ncol(abunds_mat)), function(i){abunds_mat[,i]})
	specs_emp <- unlist(lapply_fun(X=emp_cols, FUN=spec_fun, extras=extra_inputs))
	# round specs_emp to 4 decimal places to match precision of cpp functions
	specs_emp <- round(specs_emp, 4)

	# check if permutation method INDUCED NAs. Like, if permutation method caused
	# some cols to be bad. This is a critical error.
	nas_emp <- is.na(specs_emp)
	nas_sim <- apply(X=specs_sim_mat, MAR=2, FUN=function(x){any(is.na(x))})
	if(any(nas_sim & !nas_emp)){
		stop("sim_fun induced errors, likely produced columns with < 3 nonzero values.")
	}

	# calculate p-values
	msg("Calculating P-values.")
	Pval <- rep(-1, length(specs_emp))
	for(i in 1:length(Pval)){
		Pval[i] <- pval_from_perms(emp=specs_emp[i], perm=specs_sim_mat[,i], tails=tails)
	}
	# adjust p-values
	Pval <- p.adjust(Pval, method=p_adj)

	# calculate SES
	msg("Calculating SES.")

	# get denominator
	if(denom_type == "global_sim"){
		d <- sd(as.vector(specs_sim_mat))
		denom <- rep(d, ncol(specs_sim_mat))
	}else if(denom_type == "species_sim"){
		denom <- apply(X=specs_sim_mat, MAR=2, FUN=sd)
	}else if(denom_type == "raw"){
		denom <- rep(1, ncol(specs_sim_mat))
	}else if(denom_type == "global_unif"){
		# make option?
		n_null <- 2000
		dseeds <- daughter_seeds(n=n_null, s=seed)
		make_null_col <- function(s, n_samp=nrow(abunds_mat)){
			set.seed(s)
			a <- runif(n_samp, min=0, max=1)
			return(a/sum(a))
		}
		null_cols <- lapply_fun(X=dseeds, FUN=make_null_col)
		null_specs <- unlist(lapply_fun(X=null_cols, FUN=spec_fun, extras=extra_inputs))
		denom <- rep(sd(null_specs), ncol(specs_sim_mat))
	}else{
		stop("Invalid denom_type.")
	}

	spec_sim_means <- apply(X=specs_sim_mat, MAR=2, FUN=mean)
	SES <- (specs_emp - spec_sim_means)/denom
	

	# format output object
	output <- data.frame(Pval, SES)
	# if col names were input, use them for output too.
	if( !is.null(colnames(abunds_mat))){
		rownames(output) <- colnames(abunds_mat)
	}

	msg("Done.")
	return(output)
}


#' pval_from_perms
#'
#' Calculates P-value for permutation tests. 
#'
#' @author John L. Darcy
#'
#' @param emp Numeric scalar. An empirical test statistic value. 
#' @param perm Numeric vector. Test statistic values similar to emp, but calculated
#'   from permuted data. 
#' @param tails integer. 
#'   \describe{
#'     \item{1:}{Left tail only.}
#'     \item{2:}{2-tailed test.}
#'     \item{3:}{Right tail only.}
#'     \item{0:}{No test, P=1.}
#'   }
#' @param threshold integer. Minimum number n of non-NA values in perm that are
#'   acceptable. If n < threshold, P=NA (DEFAULT: 100).
#'
#' @return a P-value.
#' 
#' @export
pval_from_perms <- function(emp, perm, tails, threshold = 100){
	n <- sum(!is.na(perm))
	if(is.na(emp) || n <= threshold){
		return(NA)
	}else{
		n_perm_blw <- sum(perm < emp)
		n_perm_abv <- sum(perm > emp)
		# if 0s, bump up to 1 because we don't actually know P=0
		if(n_perm_blw < 1){n_perm_blw <- 1}
		if(n_perm_abv < 1){n_perm_abv <- 1}
		# P-value for each tails case
		if(tails == 1){
			# left tail only
			return(n_perm_blw / n)
		}else if(tails == 3){
			# right tail only
			return(n_perm_abv / n)
		}else if(tails == 2){
			# 2-tailed
			return(min(c( n_perm_blw / n, n_perm_abv / n )))
		}else{
			# no test
			return(1)
		}
	}
}


#' daughter_seeds
#' 
#' Makes n daughter seeds from seed s. This is useful for processes one wishes to be
#' deterministic, but may not be executed in the same order every time.
#' 
#' @author John L. Darcy
#' 
#' @param n integer. Number of daughter seeds to make.
#' @param s integer. A seed (DEFAULT: 12345).
#' 
#' @return vector of length n containing integer seeds.
#' 
#' @export
daughter_seeds <- function(n,s=12345){
	set.seed(s)
	return(replicate(n, as.integer(paste(sample(0:9, nchar(s), replace=TRUE), collapse=""))))
}