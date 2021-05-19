#' wpd_table
#'
#' Calculates phylogenetic entropy (Hp) for each column vector s of species
#' observations within matrix m, weighted by the frequency of each species within
#' s. Can also calculate Faith's PD.
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
#' @param m numeric matrix or data.frame of weights, where columns are species and
#'   rows are samples.
#' @param s_names species names for m if not colnames(m). NULL will use colnames (default: NULL)
#' @param s_phylo phylo object. Tree containing all unique names in s as tips.
#'   Must not contain duplicate tip labels.
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
#' @param ncores integer. Number of CPU cores to use for parallel operations (default: 4).
#' 
#' @return multiple WPD or PD values, one for each column of m.
#'
#' @examples
#' # library(specificity)
#' # set.seed(12345)
#' # s_phylo <- get(data(endophyte))$supertree
#' # w <- sample(c(0, 1), replace=TRUE, size=10)
#' # nspec <- 12
#' # m <- t(as.matrix(data.frame(
#' #   a=runif(nspec, 0, 100),
#' #   b=runif(nspec, 0, 100),
#' #   c=runif(nspec, 0, 100)
#' # )))
#' # colnames(m) <- sample(s_phylo$tip.label, ncol(m))
#' # wpd_table(m, s_phylo)
#'
#' @export
wpd_table <- function(m, s_phylo, s_names=NULL, nested_set=NULL, metric="Hp", ncores=4){
	# handle s_names == NULL
	if(is.null(s_names)){s_names <- colnames(m)}
	# check all names occur in tree
	if( ! all(s_names %in% s_phylo$tip.label) ){
		stop("Not all s_names (or colnames(m) if s_names was NULL) are in s_phylo")
	}
	# make nested set?
	if(is.null(nested_set)){
		s_phylo <- ape::keep.tip(phy=s_phylo, tip=s_names)
		nested_set <- make_nested_set(s_phylo, ncores)
	}
	samp_list <- as.list(as.data.frame(t(m)))

	if(ncores <= 1){
		output <- sapply(
			X=samp_list, 
			FUN=function(w){
				wpd(s=s_names, s_phylo=s_phylo, w=w, nested_set=nested_set, metric=metric) 
			}
		)
	}else{
		output <- simplify2array(parallel::mclapply(
			X=samp_list, 
			FUN=function(w){
				wpd(s=s_names, s_phylo=s_phylo, w=w, nested_set=nested_set, metric=metric) 
			},
			mc.cores=ncores
		))
	}
	return(output)
}