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
#' @param m matrix of species observation vectors (s). See s argument of wpd().
#' @param s_phylo phylo object. Tree containing all unique names in s as tips.
#'   Must not contain duplicate tip labels.
#' @param nested_set matrix. The output of make_nested_set(s_phylo). If not
#'   provided, will be calculated on the fly. Precalculation only provides speedup
#'   with very large trees (DEFAULT: NULL). 
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
#' @return multiple WPD or PD values, one for each column of m.
#'
#' @examples
#'   none yet written.
#'
#' @export
wpd_table <- function(m, s_phylo, nested_set, metric="Hp", ncores=4){
	wpd_vec <- function(x, names, s_phylo, nested_set, metric){
		wpd(s=names, s_phylo=s_phylo, w=x, nested_set=nested_set, metric=metric)
	}

	samp_list <- as.list(as.data.frame(m))

	output <- simplify2array(mclapply2(X=samp_list, FUN=wpd_vec, names=rownames(m), s_phylo=s_phylo, nested_set=nested_set, metric=metric, mc.cores=ncores))
	return(output)
}

