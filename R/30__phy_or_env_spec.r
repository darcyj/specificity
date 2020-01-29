#' phy_or_env_spec
#'
#' Calculates species' specificities to either a 1-dimensional variable (vector), 
#' 2-dimensional variable (matrix), or to a phylogeny. Transforms all variable
#' input types into a matrix D, and calculates specificity by comparing empirical
#' RQE* (weighted mean of D's lower triangle and unique pairwise products of 
#' species abundances) to simulated RQE* (same but with permuted abundances). This
#' "raw" specificity is then divided by some denominator d in order to standardize
#' it such that specificities can be compared between different species and
#' different variables. Values closer to 0 indicate random assortment (null
#' hypothesis), and more negative values indicate stronger specificity.
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
#' @param p_method string. method argument to pval_from_perms (DEFAULT: "raw").
#' @param denom_type string. Type of denominator (d) to use (DEFAULT: "index"). Note
#'   that denominator type does NOT affect P-values.
#'   \describe{
#'     \item{"ses":}{
#'       d for species s is calculated as the standard deviation of RQE* values
#'       calculated from permuted species weights. This makes the output specificity
#'       a standardized effect size (SES). Unfortunately, this makes SES 
#'       counterintuitively sensitive to occupancy, where species with high occupancy
#'       have more extreme SES than rare species, due to their more deterministic sim
#'       specificities. Included for comparative purposes, not suggested.
#'     }
#'     \item{"global_unif":}{
#'       d is same for all species. Calculated as variability in RQE* under random
#'       uniform distribution (beta 1,1) of species abunds. This d is comparable between
#'       different abundance matrices and between different variables. Specificity is an
#'       SES using this denom_type, and insn't sensitive to occupancy like ses is. 
#'       Fairly sensitive to sample size (number of data points per species), so this is
#'       a better option than ses if you really want units of SDs, but is still not 
#'       suggested.
#'     }
#'     \item{"raw":}{
#'       d is 1 for all species, so output specificity has units of distance, i.e. the
#'       raw difference between empirical and simulated RQE*. This means that results
#'       from different variables are not comparable, since it is not scale-invariant to
#'       env or hosts_phylo. It IS still scale-invariant to the species weights in 
#'       aunds_mat. Not sensitive to number of samples. Not suggested because units are
#'       strange, and isn't comparable between variables. 
#'     }
#'     \item{"index":}{
#'       d is the mean of simulated (permuted) RQE* values for species that have stronger
#'       specificity than expected by chance, resulting in specificity values with range
#'       [-1, 0), with 0 as the null hypothesis. In this case, -1 indicates perfect
#'       specificity, where a species is associated with zero environmental variability.
#'       In the euclidean sense, this could be a species that is always found at the
#'       exact same elevation or the exact same pH. For species that have weaker specificity
#'       than expected by chance, d is x - the mean of simulated RQE* values, where x is
#'       the maximum possible dissimilarity. For vector inputs, this is the vector's range.
#'       For phylogenetic inputs, it's the greatest observed patristic distance. For matrix
#'       inputs, it's assumed to be 1, which can be overridden using matrix_tmax. In this
#'       way, overdispersed values have range (0, 1]. This d has other useful properties:
#'       scale invariance to env/hosts_phylo, insensitivity to the number of samples, 
#'       insensitivity to occupancy, and strong sensitivity to specificity (DEFAULT).
#'   }
#' @param matrix_tmax float. Theoretical maximum value of distances when a matrix or dist
#'   is supplied as env. Only used with denom_type = "index" (the default); if so an error
#'   will be thrown if there is a value inside a supplied matrix greater than matrix_tmax
#'   (DEFAULT: 1). 
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
	p_adj="fdr", seed=1234567, tails=1, n_cores=2, verbose=TRUE, lowmem=FALSE, 
	p_method="raw", denom_type="index", matrix_tmax=1, diagnostic=F){

	# testing stuff:
	# abunds_mat = data.frame(z_flat, z_elev_lo, z_elev_med, 
	#   z_elev_hi, z_elev_lo_all, z_elev_med_all, z_elev_hi_all)
	# env = metadata$Elevation
	# n_sim=100
	# n_cores=10
	# hosts=NULL
	# hosts_phylo=NULL
	# sim_fun=function(m){ m[sample(1:nrow(m)), ] }
	# p_adj="fdr"; seed=1234567; tails=1; verbose=TRUE; lowmem=FALSE
	# denom_type="index"
	# p_method <- "raw"


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

	# check matrix_tmax is ok
	if(data_type %in% c("mat", "dist") && denom_type == "index"){
		if(max(env) > matrix_tmax){
			stop("Values greater than matrix_tmax are present in env.")
		}
	}

	# turn it all into a dist, no matter what.
	if(data_type == "phy"){
		msg("Converting tree to dist.")
		env <- tree2mat(tree=hosts_phylo, x=hosts, n_cores=n_cores, delim=";")
		tmax <- max(env) # only used for "index"
	}else if(data_type == "vec"){
		msg("Converting env vector to dist.")
		tmax <- max(env) - min(env) # only used for "index"
		env <- dist(env)
	}else if(data_type == "mat"){
		env <- as.dist(env)
		tmax <- matrix_tmax # only used for "index"
	}else if(data_type == "dist"){
		tmax <- matrix_tmax # only used for "index"
	}

	# define specificity function, necessary because rao_quad_ent has args in wrong order.
	spec_fun <- function(ab, e){
		return(rao_quad_ent(d=e, w=ab, raw=FALSE))
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
			FUN=spec_fun, e=env))
		if(verbose){setTxtProgressBar(pb, b)}
	}
	if(verbose){close(pb)} # to kill of progress bar
	# turn it into a matrix, where each row is a simulation, and each column is a species.
	specs_sim_mat <- matrix(specs_sim, ncol=ncol(abunds_mat), byrow = TRUE )

	# use spec_fun on empirical data
	msg("Calculating empirical RQE*.")
	emp_cols <- lapply(seq_len(ncol(abunds_mat)), function(i){abunds_mat[,i]})
	specs_emp <- unlist(lapply_fun(X=emp_cols, FUN=spec_fun, e=env))
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
		Pval[i] <- pval_from_perms(emp=specs_emp[i], perm=specs_sim_mat[,i], 
			tails=tails, method=p_method, rounding=4)
	}
	# adjust p-values
	Pval <- p.adjust(Pval, method=p_adj)

	# calculate output specificity
	msg("Calculating specificities.")
	spec_sim_means <- apply(X=specs_sim_mat, MAR=2, FUN=mean)

	# get denominator
	if(denom_type == "ses"){
		denom <- apply(X=specs_sim_mat, MAR=2, FUN=sd)
	}else if(denom_type == "raw"){
		denom <- rep(1, ncol(specs_sim_mat))
	}else if(denom_type == "global_unif"){
		# todo: make option?
		n_null <- 2000
		dseeds <- daughter_seeds(n=n_null, s=seed)
		make_null_col <- function(s, n_samp=nrow(abunds_mat)){
			set.seed(s)
			a <- runif(n_samp, min=0, max=1)
			return(a)
		}
		null_cols <- lapply_fun(X=dseeds, FUN=make_null_col)
		null_specs <- unlist(lapply_fun(X=null_cols, FUN=spec_fun, e=env))
		denom <- rep(sd(null_specs), ncol(specs_sim_mat))
	#}else if(denom_type == "flat"){
	#	d <- spec_fun(rep(1, nrow(abunds_mat)), e=env)
	#	denom <- rep(d, ncol(specs_sim_mat))
	}else if(denom_type == "index"){
		denom <- spec_sim_means
		denom[specs_emp >  spec_sim_means] <- tmax
	}else{
		stop("Invalid denom_type.")
	}

	out_specs <- round(specs_emp - spec_sim_means, 4)/ round(denom, 4)
	
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
#' @param method string. Method by which P should be calculated from perms:
#'     \item{"raw":}{
#'       P is calculated as the sum of sim values more extreme than the empirical
#'       value plus one, divided by the number of sim values. 
#'     }
#'     \item{"dens_fit":}{
#'       P is calculated via kernel density estimation. No better than "raw".
#'     }
#'     \item{"gamma_fit":}{
#'       P is calculated by fitting a gamma distribution to sim values and calculating 
#'       area under the curve from (-inf,emp] or [emp,inf) depending on tailedness.
#'     }    
#' @param threshold integer. Minimum number n of non-NA values in perm that are
#'   acceptable. If n < threshold, P=NA (DEFAULT: 50).
#' @param rounding integer. Number of decimal places to round emp and perm This is only
#'   useful when emp and perm are expected to contain the exact same value, but the number
#'   of decimal places in that value is different between emp and perm Use a number less
#'   then zero to disable rounding (DEFAULT: -1).
#'
#' @return a P-value.
#' 
#' @export
pval_from_perms <- function(emp, perm, tails, method="raw", threshold=30, rounding=-1){
	require("fitdistrplus")
	# do rounding...?
	if(rounding >= 0){
		emp <- round(emp, rounding)
		perm <- round(perm, rounding)
	}
	# check if number of perms is below threshold
	n <- sum(!is.na(perm))
	if(n < threshold){
		LTP <- NA
		RTP <- NA
	}else if(method == "dens_fit"){
		# remove NAs
		perm <- perm[!is.na(perm)]
		frm <- min(c(mean(perm) - (mean(perm) - min(perm)) * 4, emp))
		too <- max(c(mean(perm) + (max(perm) - mean(perm)) * 4, emp))
		dns <- density(perm, from=frm, to=too, bw="sj")
		# Riemann sum
		dx <- dns$x[2] - dns$x[1]
		tot <- sum(dns$y * dx)
		LTP <- sum(dns$y[dns$x < emp] * dx) / tot
		RTP <- 1-LTP
	}else if(method == "raw"){
		# remove NAs
		perm <- perm[!is.na(perm)]
		# calculate P for right and left tails
		LTP <- (sum(perm <= emp) + 1)/n
		RTP <- (sum(perm >= emp) + 1)/n
	}else if(method == "gamma_fit"){
		gfit <- fitdistrplus::fitdist(data=perm, distr="gamma", lower=c(0,0))
		shp <- gfit$estimate[names(gfit$estimate) == "shape"]
		rte <- gfit$estimate[names(gfit$estimate) == "rate"]
		LTP <- pgamma(q=emp, shape=shp, rate=rte)
		RTP <- 1 - LTP
	}else{
		stop("Invalid method argument.")
	}
	# fix P > 1 in cases where all values of perm == emp
	if(LTP > 1){LTP <- 1}
	if(RTP > 1){RTP <- 1}
	# calculate P for tailed situations
	if(tails == 1){
		return(LTP)
	}else if(tails == 3){
		return(RTP)
	}else if(tails == 2){
		return(min(c(LTP, RTP)))
	}else{
		stop("Invalid tails argument.")
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


#' check_pes_inputs
#' 
#' Function used by phy_or_env_spec. 
#' checks abunds_mat, env, hosts, and hosts_phylo inputs to phy_or_env_spec to
#' make sure there are no problems. This could include missing species in trees,
#' incompatible dimensions, non-numeric inputs, etc. Returns an input type, which
#' is just a string that can be "mat", "dist", "vec", "phy", or "error".
#' 
#' @param abunds_mat (required, see phy_or_env_spec)
#' @param env (required, can be NULL, see phy_or_env_spec)
#' @param hosts (required, can be NULL, see phy_or_env_spec)
#' @param hosts_phylo (required, can be NULL, see phy_or_env_spec)
#' @param verbose logical. Should status messages be displayed? (DEFAULT: TRUE).
#' 
#' @return string. either "mat", "dist", "vec", "phy", or "error".
#' 
#' @export
check_pes_inputs <- function(abunds_mat, env, hosts, hosts_phylo, verbose=TRUE){
	msg <- function(x){if(verbose){message(x)}}
	# decide inputs type
	env_prsnt <- !is.null(env)
	hos_prsnt <- !is.null(hosts)
	phy_prsnt <- !is.null(hosts_phylo)
	if(env_prsnt && (hos_prsnt || phy_prsnt)){
		# both env AND some phylogenetic stuff present, error!
		data_type <- "error"
		msg("Error: Too many inputs. Pick either env only, or hosts+hosts_phylo.")
	}else if(env_prsnt && (is.matrix(env) || is.data.frame(env))){
		# env present and it's 2-dimensional
		data_type <- "mat"
	}else if(env_prsnt && class(env) == "dist"){
		# env present and it's a dist (LT from 2-dimensional)
		data_type <- "dist"
	}else if(env_prsnt && is.vector(env)){
		# env present and it's a vector
		data_type <- "vec"
	}else if(env_prsnt && is.factor(env)){
		data_type <- "error"
		msg("Error: env is a factor.")
	}else if(hos_prsnt && phy_prsnt){
		# hosts and hosts_phylo present, it's phy
		data_type <- "phy"
	}else if(xor(hos_prsnt, phy_prsnt)){
		# only one of hosts or hosts_phylo are present, but no env. error!
		data_type <- "error"
		msg("Error: Both hosts and hosts_phylo required for phylogenetic analysis.")
	}else if(! env_prsnt || hos_prsnt){
		# none of them are present!
		data_type <- "error"
		msg("Error: Either env or hosts is required.")
	}else{
		# mystery error!
		data_type <- "error"
		msg("Error: Mystery error.")
	}

	# check for errors specific to 2D data:
	if(data_type %in% c("mat", "dist")){
		# if it's a dist, just convert to mat. all same checks apply.
		if(data_type == "dist"){
			env <- as.matrix(env)
		}
		# is env numeric?
		if(!is.numeric(env)){
			data_type <- "error"
			msg("Error: env not numeric.")
		# is env square?
		}else if(nrow(env) != ncol(env)){
			data_type <- "error"
			msg("Error: env is 2-dimensional, but not square.")
		# is env same dim as abunds_mat?
		}else if(nrow(env) != nrow(abunds_mat)){
			data_type <- "error"
			msg("Error: env and abunds_mat have incompatible dimensions.")
		}
	}

	# check for errors specific to 1D env data:
	if(data_type == "vec"){
		# is env numeric?
		if(!is.numeric(env)){
			data_type <- "error"
			msg("Error: env not numeric.")
		# is env same dim as abunds_mat?
		}else if(length(env) != nrow(abunds_mat)){
			data_type <- "error"
			msg("Error: env and abunds_mat have incompatible dimensions.")
		}
	}

	# check for errors specific to phylogenetic data:
	if(data_type == "phy"){
		# are all the right tips in hosts_phylo?
		if(! all(hosts %in% hosts_phylo$tip.label)){
			data_type <- "error"
			msg("Some hosts are missing from hosts_phylo.")
		# is hosts the right length?
		}else if(length(hosts) != nrow(abunds_mat)){
			data_type <- "error"
			msg("hosts and abunds_mat have incompatible dimensions.")
		# do any hosts contain semicolons?
		}else if(any(grepl(";", x=hosts_phylo$tip.label))){
			data_type <- "error"
			msg("hosts cannot contain semicolons.")
		}
	}

	return(data_type)
}


