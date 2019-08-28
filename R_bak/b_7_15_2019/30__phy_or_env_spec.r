#' phy_or_env_spec
#'
#' Calculates standardized effect size (SES) of host phylogenetic specificity (per 
#' Poulin et al. 2011) or standardized effect size of environmental specificity
#' for each species (column) in a matrix. 
#'
#' @author John L. Darcy
#' @references 
#' Poulin et al. (2011) Host specificity in phylogenetic and geographic
#'   space. Trends Parasitol 8:355-361. doi: 10.1016/j.pt.2011.05.003
#'
#' @param abunds_mat matrix or data frame of numeric values. Columns represent 
#'   species, rows are samples. For columns where the value is nonzero for two or
#'   fewer data points, environmental SES cannot be calculated, and NAs will be 
#'   returned. Negative values in abunds_mat are not allowed.
#' @param hosts character vector. Host identities corresponding to abunds. Only
#' required if calculating SES for phylogenetic specificity.
#' @param hosts_phylo phylo object. Tree containing all unique hosts as tips. Only
#' required if calculating SES for phylogenetic specificity.
#' @param env numeric vector. Environmental variable corresponding to abunds. For
#'   example, temperature. Only required if calculating SES for environmental
#'   specificity.
#' @param pd_metric character. Type of phylogenetic diversity to calculate for SPS.
#'   See ?wpd (DEFAULT: "Hp").
#' @param n_sim integer. Number of simulations of abunds_mat to do under the null 
#'   hypothesis that host or environmental association is random DEFAULT: 200).
#' @param sim_fun function. A function f(m) where f(abunds_mat) returns a matrix
#'   object with the same number of rows and columns as abunds_mat. Default is
#'   f=function(m){ m[sample(1:nrow(m)),]}, which just permutes the order of 
#'   rows in abunds_mat. Users may wish to use a null model that is able to
#'   preserve row and column totals such as the function permatswap() from the
#'   vegan package or the function vaznull() from the bipartite package. Either
#'   of these can be easily adapted to return only a single matrix (see examples). 
#'   However, neither can accomodate non-integer matrices. See adj_fun for using 
#'   non-integer abundance values. (DEFAULT: see above).
#' @param adj_fun function. A function f where f(abunds_mat) or 
#'   f(sim_fun(abunds_mat)) returns an object of class matrix with the same number
#'   of rows and columns as abunds_mat. This function can be some adjustment for
#'   integer values, or some rarefaction function, et cetera. It is applied to any 
#'   empirical or simulated matrix during analysis. Default behavior is to do no
#'   adjustment (DEFAULT: mat_passthru).
#' @param pair_weights_fun function. A function that computes pair-wise weight w12 given
#'   two weights w1 and w2. This is only used if env is a distance matrix. Since
#'   each column vector c in abunds_mat is 1-dimensional, but env is 2-dimensional,
#'   c must be transformed into a 2-dimensional object of the same dimensions as
#'   env. This is done via pairwise calculation, so a new weight is calculated for
#'   each unique pair-wise combination of c. Default is product of the two weights
#'   (DEFAULT: function(w1, w2){w1 * w2} ).
#' @param p_adj logical. Should output P-values be adjusted for multiple hypothesis
#'   testing? If true, will be done using FDR (Default = TRUE).
#' @param seed integer. Seed to use so that this is repeatable (DEFAULT: 1234557)
#' @param tails integer. 1 = 1-tailed, test for specificity only. 2 = 2-tailed.
#'   3 = 1-tailed, test for cosmopolitanism only. 4 = no test, P=1.0 (DEFAULT: 1).
#' @param n_cores integer. Number of CPU cores to use for parallel operations. If
#'   set to 1, lapply will be used instead of mclapply (DEFAULT: 2).
#' @param verbose logical. Should status messages be displayed? (DEFAULT: TRUE).
#' @param lowmem logical. Should this function be run in a way that saves memory?
#'   If TRUE, will run FAR slower but use almost no ram. If FALSE, lots of ram is
#'   used but it will run quickly. For example, on an 800x1600 matrix, with nperm=
#'   1000, n_cores=12 and lowmem=F, will use roughly 90 GB ram but finish in one
#'   minute. With the same settings but lowmem=T, will use less than 1 GB ram but 
#'   take over an hour. A progress bar is shown if TRUE (DEFAULT: FALSE).
#'
#' @return data.frame where each row is an input species. First column is P-value,
#'   second column is SES.
#'
#' @examples
#'	# phylogenetic specificity using phylocom data set
#'	# using default naive null model (permutation; Poulin et al. 2011).
#'	library(picante); data(phylocom)
#'	m <- t(phylocom$sample)
#'	ses_naive <- phy_or_env_spec(
#'		abunds_mat=m,
#'		hosts=rownames(m),
#'		hosts_phylo=phylocom$phylo,
#'		verbose=TRUE,
#'		n_cores=12
#'	)
#'
#'	# phylogenetic specificity using phylocom data set
#'	# using permatswap null model from vegan package
#'	library(vegan)
#'	ses_pms <- phy_or_env_spec(
#'		abunds_mat=m,
#'		hosts=rownames(m),
#'		hosts_phylo=phylocom$phylo,
#'		sim_fun=function(m){permatswap(m)$perm[[1]]},
#'		verbose=TRUE,
#'		n_cores=12
#'	)
#'
#'	# phylogenetic specificity using phylocom data set
#'	# using vazquez null model from bipartite package
#'	ses_vaz <- multi_phy_or_env_spec(
#'		abunds_mat=m,
#'		hosts=rownames(m),
#'		hosts_phylo=phylocom$phylo,
#'		sim_fun=function(m){bipartite::vaznull(1, m)[[1]]},
#'		verbose=TRUE,
#'		n_cores=12
#'	)
#'
#'	# compare results
#'  par(mfrow=c(2,1))
#'	plot(ses_naive$SES, ses_pms$SES, ylab="vegan::permatswap", xlab="naive")
#'  plot(ses_naive$SES, ses_vaz$SES, ylab="bipartite::vaznull", xlab="naive")
#'  
#'
#'
#' @export
	phy_or_env_spec <- function(abunds_mat, env=NULL, hosts=NULL, 
		hosts_phylo=NULL, pd_metric="Hp", n_sim=1000, 
		sim_fun=function(m){ m[sample(1:nrow(m)), ] }, adj_fun=mat_passthru, 
		seed=1234567, tails=1, n_cores=2, verbose=TRUE, 
		pair_weights_fun=function(x,y){sqrt(x*y)}, p_adj=TRUE, lowmem=FALSE){

		require(ape)
		require(geiger)

		# little message function, to clean up code a bit
		msg <- function(x){if(verbose){message(x)}}

		# decide which input type to do, and figure out if inputs look OK.
		# also define specificity function.
		msg("Checking inputs.")
		spec_type <- "error"
		# initialize spec_fun and ns just in case
		# ns is pre-computed nested set, which may or may not be calculated.
		spec_fun <- function(x){"ERROR: spec_fun not defined."}
		ns <- NULL
		if((is.null(hosts) || is.null(hosts_phylo)) && !is.null(env)){
			# either hosts or tree is null, and env is provided. Doing env spec.
			# define function for weighted sample variance
			w_var <- function(x, w){
				if(sum(w > 0) < 2){
					return(NA)
				}else{
					return(sum(w * (x - weighted.mean(x,w))^2) / (sum(w) - (sum(w^2 / sum(w)))) )
				}
			}

			# but is it a matrix or a vector?
			if(class(env) %in% c("matrix", "dist", "data.frame")){
				env <- as.dist(env)
				if( length(env) != (nrow(abunds_mat)^2 - nrow(abunds_mat))/2 ){
					stop("ERROR: abunds_mat and env have incompatible dimensions.")
				}
				spec_type <- "env_dist"
				spec_fun <- function(ab, env=NULL, hosts=NULL, hosts_phylo=NULL, pd_metric=NULL, pair_fun=pair_weights_fun, nested_set=NULL){
					ab_dist <-  as.dist(outer(ab, ab, pair_fun))
					return(weighted.mean(x=as.vector(env), w=as.vector(ab_dist)))
				}
			}else{
				if( length(env) != nrow(abunds_mat)){
					stop("ERROR: abunds_mat and env have incompatible dimensions.")
				}
				spec_type <- "env_vec"
				spec_fun <- function(ab, env=NULL, hosts=NULL, hosts_phylo=NULL, pd_metric=NULL, pair_fun=NULL, nested_set=NULL){
					return(w_var(x=env, w=ab))
				}
			}
		}else if(is.null(env) && (!is.null(hosts) && !is.null(hosts_phylo))){
			# env is null and hosts was given. doing phy spec.
			if(length(hosts) != nrow(abunds_mat)){
				stop("ERROR: abunds_mat and hosts have incompatible dimensions.")
			}else if( ! all(hosts %in% hosts_phylo$tip.label)){
				stop("ERROR: Not all hosts are in hosts_phylo as tips.")
			}
			spec_type <- "phy"
			# precompute nested set
			ns <- make_nested_set(phy=hosts_phylo, n_cores=n_cores)
			
			spec_fun <- function(ab, env=NULL, hosts=NULL, hosts_phylo=NULL, pd_metric=NULL, pair_fun=NULL, nested_set=NULL){
				return(wpd(s=hosts, s_phylo=hosts_phylo, w=ab, nested_set=ns, metric=pd_metric))
			}

		}else if(!is.null(hosts) && is.null(hosts_phylo)){
			stop("ERROR: hosts provided but hosts_phylo is missing.")
		}else if(is.null(env) && is.null(hosts)){
			stop("ERROR: hosts or env is missing.")
		}else{
			stop("ERROR: mystery input error.")
		}

		# make wrapper for lapply/mclapply so with cores=1 it just uses lapply
		if(n_cores > 1){
			require("parallel")
			lapply_fun <- function(X, FUN, ...){mclapply(X, FUN, mc.cores=n_cores, ...)}
		}else{
			lapply_fun <- function(X, FUN, ...){lapply(X, FUN, ...)}
		}

		msg("Adjusting input matrix with adj_fun.")
		abunds_mat <- adj_fun(abunds_mat)


		msg("Generating daughter seeds.")
		# generate n_sim daughter seeds for generation of each permuted matrix
		set.seed(seed)
		seeds <- replicate(n=n_sim, as.integer(paste(sample(0:9, nchar(seed), replace=TRUE), collapse="")))

		# function to generate matrix given seed:
		perm_mat <- function(seed){
			set.seed(seed)
			return((sim_fun(abunds_mat)))
		}

		if(lowmem == FALSE){
			# create n_sim permuted matrices, but convert them into a list of column vectors.
			msg(paste("Creating", n_sim, "permuted matrices."))
			perm_cols <- do.call("cbind", lapply_fun(FUN=perm_mat, X=seeds))
			perm_cols <- lapply(seq_len(ncol(perm_cols)), function(i){perm_cols[,i]})

			# apply spec_fun to each column vector in perm_cols
			msg("Calculating specificities for permuted matrices.")
			specs_sim <- unlist(lapply_fun(X=perm_cols, FUN=spec_fun, env=env, 
				hosts=hosts, hosts_phylo=hosts_phylo, pd_metric=pd_metric, 
				pair_fun=pair_weights_fun, nested_set=ns))
			# turn it into a matrix, where each row is a simulation, and each column is a species.
			specs_sim_mat <- matrix(specs_sim, nrow=n_sim, ncol=ncol(abunds_mat) )
		}else{
			msg("Calculating specificities for permuted matrices (lomem).")
			specs_sim_mat <- matrix(0, nrow=n_sim, ncol=ncol(abunds_mat))
			if(verbose){pb <- txtProgressBar(min=0, max=n_sim, style=3)}
			nc <- ncol(abunds_mat)
			for(i in 1:n_sim){
				# permute, adjust abunds_mat, and split it into list of columns.
				pm <- perm_mat(seed=seeds[i])
				perm_cols <- lapply(seq_len(nc), function(i){pm[,i]})
				specs_sim_mat[i,] <- unlist(lapply_fun(X=perm_cols, FUN=spec_fun, env=env, hosts=hosts, hosts_phylo=hosts_phylo, pd_metric=pd_metric, pair_fun=pair_weights_fun, nested_set=ns))
				if(verbose){setTxtProgressBar(pb, i)}
			}
		}

		# use spec_fun on empirical data
		msg("Calculating empirical specificities.")
		emp_cols <- lapply(seq_len(ncol(abunds_mat)), function(i){abunds_mat[,i]})
		specs_emp <- unlist(lapply_fun(X=emp_cols, FUN=spec_fun, env=env, hosts=hosts, hosts_phylo=hosts_phylo, pd_metric=pd_metric, pair_fun=pair_weights_fun, nested_set=ns))

		# check if permutation method INDUCED NAs. Like, if permutation method caused
		# some cols to be bad. This is a critical error.
		nas_emp <- is.na(specs_emp)
		nas_sim <- apply(X=specs_sim_mat, MAR=2, FUN=function(x){any(is.na(x))})
		if(any(nas_sim & !nas_emp)){
			stop("ERROR: sim_fun induced errors, likely produced columns with < 3 nonzero values.")
		}

		# calculate p-values
		msg("Calculating P-values.")
		get_pval <- function(i, tails, n){
			# handle NAs
			if(is.na(specs_emp[i])){
				return(NA)
			}else{
				n_sims_blw <- sum(specs_sim_mat[,i] < specs_emp[i])
				n_sims_abv <- sum(specs_sim_mat[,i] > specs_emp[i])

				# if n=0, bump it up to 1 because we don't actually know P=0
				if(n_sims_blw < 1){n_sims_blw <- 1}
				if(n_sims_abv < 1){n_sims_abv <- 1}
				# P-value for each tails case:
				if(tails == 1){ # left tail, underdispersion
					return(n_sims_blw / n)
				}else if(tails == 3){ # right tail, overdispersion
					return(n_sims_abv / n)
				}else if(tails == 2){ # 2-tailed
					return(min(c( (n_sims_blw / n), (n_sims_abv / n) )))
				}else{ # no test
					return(1)
				}
			}
		}
		Pval <- sapply(X=1:ncol(abunds_mat), FUN=get_pval, tails=tails, n=n_sim)

		# calculate SES
		msg("Calculating SES.")
		denom <- sd(as.vector(specs_sim_mat))

		spec_sim_means <- apply(X=specs_sim_mat, MAR=2, FUN=mean)
		SES <- (specs_emp - spec_sim_means)/denom

		# format output object
		output <- data.frame(Pval, SES)
		# if col names were input, use them for output too.
		if( !is.null(colnames(abunds_mat))){
			rownames(output) <- colnames(abunds_mat)
		}

		# adjust p-values
		if(p_adj == TRUE){
			output$Pval <- p.adjust(output$Pval, method="fdr")
		}

		msg("Done.")
		return(output)
	}
