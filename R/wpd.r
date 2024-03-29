#' wpd
#'
#' Calculates weighted Phylogenetic Diversity for a vector s of species
#' observations, weighted by the frequency of each species within s. For example,
#' if s={a, a, b, a, b, c, a}, then species a will have weight 4, species b will
#' have weight 2, and species c will have weight 1. Unobserved species have weight
#' zero. However, one may wish to exclude observations that do not meet some criterion,
#' such as co-observation of a symbiote or parasite. For this reason, a second set
#' of weights w can be provided as a vector of numeric values that are paired with
#' s. These weights are then implicitely combined with the weights discussed above 
#' depending on which weighted metric is chosen. In the case of Phylogenetic 
#' Entropy (Hw), per-tip weights are calculated as the sums of w. In the case of
#' Weighted Faith (WF), per-tip weights are averages of w. 
#'
#' @author John L. Darcy
#' @references
#' \itemize{
#'   \item Allen B, Kon M, Bar-Yam Y (2009) A new phylogenetic diversity measure
#'     generalizing the Shannon index and its application to Phyllostomid bats.
#'     American Naturalist 174(2).
#'   \item Swenson NG (2014) Functional and Phylogenetic Ecology in R. 
#'     Springer UseR! Series, Springer, New York, New York, U.S.A.
#'   \item Faith DP (1992) Conservation evaluation and phylogenetic diversity. 
#'     Biological Conservation 61.
#' }
#'
#' @param s character vector. One species name per observation. If no species was
#'   observed for a given datum, use NA. s can also be provided as a vector of
#'   unique species identities, in which case counts of those species can be
#'   given as w.
#' @param s_phylo phylo object. Tree containing all unique names in s as tips.
#'   Must not contain duplicate tip labels.
#' @param w numeric vector. Optional weights for s, e.g. number of parasites 
#'   observed in each sample, or boolean weights corresponding to presence or 
#'   absence of parasite species, or confidence species was observed, etc. If w is
#'   not provided but a weighted metric is specified, w will be set to 1 for each
#'   value of s. Thus, weights for each unique species in s would be equal to the
#'   number of times that species appears in s. w is not used for unweighted
#'   metrics (PD). Any NA values in w will be pairwise removed from w and s
#'   (default: NULL).
#' @param nested_set matrix. The output of make_nested_set(s_phylo). If not
#'   provided, will be calculated on the fly. Precalculation only provides speedup
#'   with very large trees (default: NULL). 
#' @param metric character. Abbreviated name of desired tree-based phylogenetic
#'   diversity metric. Available metrics are:
#'   \describe{
#'     \item{Hp:}{
#'       Phylogenetic Entropy. Insensitive to 0 weights, cannot increase with removal
#'       of taxa. Allen et al. 2009.
#'     }
#'     \item{WF:}{
#'       Weighted Faith's PD. Sensitive to 0 weights, i.e. a clade that was heavily
#'       sampled but has lots of zeroes will cause its sister clades to be
#'       underrepresented. Swenson 2014.
#'     }
#'     \item{PD:}{
#'       Original Faith's Phylogenetic Diversity. Unweighted. Simply a sum of branch-
#'       lengths in your tree (but only for taxa in s). Faith 1992.
#'     }
#'   }
#'
#' @return Single WPD or PD value.
#'
#' @examples
#' # library(specificity)
#' # set.seed(12345)
#' # s_phylo <- get(data(endophyte))$supertree
#' # w <- sample(c(0, 1), replace=TRUE, size=10)
#' # s <- sample(s_phylo$tip.label, replace=TRUE, size=10)
#' # wpd(s, s_phylo, w, metric="Hp")
#'
#' @export
wpd <- function(s, s_phylo, w=NULL, nested_set=NULL, metric="Hp"){

	# if no w provided, make fake w (all 1s) 
	if(is.null(w)){
		w <- rep(1, length(s))
	}

	# check that s and w are same length
	if( length(s) != length(w) ){
		stop("ERROR: Vectors s and w are not same length.")
	}

	# check that all unique values of s are in s_phylo
	if( ! all(unique(s) %in% s_phylo$tip.label)){
		stop("ERROR: Not all values of s are in s_phylo as tips.")
	}

	# remove NA data
	keepers <- (!is.na(w)) & (!is.na(s))
	w <- w[keepers]
	s <- s[keepers]

	# make sure tip names in s_phylo are unique
	if(any(duplicated(s_phylo$tip.label))){
		stop("ERROR: Duplicate tip labels in s_phylo.")
	}

	# if no nested set provided make one. necessary even for regular PD.
	if(is.null(nested_set)){
		nested_set <- make_nested_set(s_phylo, n_cores=1)
	}

	# aggregate raw weights for each tip in s_phylo. returns NA for tips that
	# are not in s, so that they can be ignored later.
	get_raw_tip_weights <- function(tipname){ 
		b <- s == tipname
		if(sum(b) == 0){
			return(NA)
		}else{
			return( sum(w[b] )) # there was an error here where sum was ommitted...?
		}
	}
	w_tips <- sapply(X=s_phylo$tip.label, FUN=get_raw_tip_weights)

	# naignore function, different from na.rm. this allows downstream functions f(x)
	# to return NA if all values of x are NA, but do the na.rm bit if only some
	# values are NA. This is necessary so that unused branches can be ignored safely.
	# the alternative is pruning the tree every time this function is used, and that's
	# a) hard with a nested set, and b) computationally expensive. So we don't want
	# f(c(NA,NA,NA)) to return 0, since that counts toward an upstream average. 
	# Instead, we want f(x) to return NA so it can be ignored. In other words, 
	# NA is being used as a second type of 0. Type 1 is a true zero, where the user 
	# wanted weight 0 for the branch. Type 2 is a fake 0, an NA, where we don't want to
	# consider that zero for functions like mean or sum. I had to write all this for
	# when I inevitably forget why I wrote the "naig" function later on.
	naig <- function(x){
		b <- is.na(x)
		if(all(b)){
			return(NA)
		}else{
			return(x[!b])
		}
	}
	sum_naig <- function(x){ sum(naig(x)) }
	mean_naig <- function(x){ mean(naig(x)) }

	# similarly, for phylogenetic entropy (secretly Shannon's index), we need a log
	# function that returns 0 for 0 * log(0). Google it. 
	plogp <- function(p){
		plogp_scalar <- function(p_i){
			if(p_i==0){
				return(0)
			}else{
				return(p_i*log(p_i))
			}
		}
		return(sapply(X=p, FUN=plogp_scalar))
	}

	# calculate final weight for each branch in nested_set matrix
	# and return result
	if(metric == "Hp"){
		# wfun returns pb, proportional weight of branch b.
		wfun <- function(nsb, rw){
			# nsb = nested set branch, a row of nested_set df
			# rw = raw weights
			return(sum_naig(rw[nsb[2]:nsb[3]]) / sum_naig(rw) )
		}

		# from Allen et al. 2009, Equation 1
		pb <- apply(X=nested_set, MARGIN=1, FUN=wfun, rw=w_tips)
		lb <- s_phylo$edge.length
		# since Hp is insensitive to zeroes (like Shannon), just turn all NAs
		# into zeroes.
		pb[is.na(pb)] <- 0
		# calculate Hp (remember, plogp(pb) means pb * log(pb) )
		return( -1 * sum( lb * plogp(pb) ))

	}else if (metric == "WF"){
		wfun <- function(nsb, rw){
			return(mean_naig(unlist(rw[ nsb[2]:nsb[3] ])))
		}
		# from Swenson 2014, pp. 36
		a <- apply(X=nested_set, MARGIN=1, FUN=wfun, rw=w_tips)
		n <- sum(!is.na(a))
		l <- s_phylo$edge.length
		return(n * sum_naig(a * l) / sum_naig(a))
	}else if (metric == "PD"){
		# even though PD is unweighted, wfun just figures out which branches
		# are to be ignored because they don't have descendents in S.
		wfun <- function(nsb, rw){
			# get descendent weights
			desc <- rw[ nsb[2]:nsb[3] ]
			if(all(is.na(desc) | desc <= 0)){
				# if all descendents of edge are NA or 0 (not in s):
				return(0)
			}else{
				# if any descendents of edge aren't NA:
				return(1)
			}
		}
		branch01 <- apply(X=nested_set, MARGIN=1, FUN=wfun, rw=w_tips)
		# from Faith 1992
		return(sum(s_phylo$edge.length * branch01))
	}else{
		stop(paste("ERROR: metric \"", metric, "\" not defined.", sep=""))
	}
}

