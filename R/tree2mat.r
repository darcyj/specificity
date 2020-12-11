
#' tree2mat
#'
#' Transforms a phylogenetic tree into a dist object containing patristic distances
#' between tips. Dists are just lower triangles of matrices, and the rows and
#' columns of that matrix are defined by a user-supplied vector of tip labels,
#' which can include duplicate values. Contrast with ape::cophenetic.phylo, which
#' produces a distance matrix containing only unique pairwise patristic distances
#' within the phylogeny. 
#'
#' @author John L. Darcy
#'
#' @param tree phylo object. Tree containing all unique species in x as tips. May
#'   contain tips that are not in x.
#' @param x character vector. Vector of species identities, each of which must
#'   be in tree as a tip label. May contain any given species identity more than
#'   once.
#' @param n_cores integer. Number of cores to use for parallel computation. No
#'   parallelization will be done if n_cores = 1. Multithreading should only be used
#'   for large trees where x has low redundancy (default: 1).
#' @param delim string. Delimiter character or string for internal use. Must not be
#'   present in tree$tip.label. This is checked by the function and will return an
#'   error otherwise (default: ";").
#' @return dist object, of vector length equal to (l^2-l)/2 where l is length(x); i.e. 
#'   values are the lower triangle of a patristic distance matrix with rows=x and cols=x. 
#'
#' @examples
#'   library(specificity)
#'   library(ape)
#'   example_tree <- ape::read.tree(text=" ((((a:1,b:1):1,c:2):1,d:3):1,(e:1,f:1):3);")
#'   example_x <- c("a", "a", "a", "b", "c", "d", "c", "a", "f")
#'   # unique patristic distance matrix:
#'   ape::cophenetic.phylo(example_tree)
#'   # dist object for example_x:
#'   tree2mat(tree=example_tree, x=example_x)
#' 
#'   # examples with other delimiters
#'   tree2mat(tree=example_tree, x=example_x, delim="@")
#'   tree2mat(tree=example_tree, x=example_x, delim="i love cats")
#'   # should fail since "a" is in a tip name:
#'   # tree2mat(tree=example_tree, x=example_x, delim="a")
#' 
#' @export
tree2mat <- function(tree, x, n_cores=1, delim=";"){
	# check that each species in x is in tree exactly once
	times_in_tree <- sapply(x, FUN=function(x){sum(tree$tip.label %in% x)})
	if(0 %in% times_in_tree){ stop("ERROR: some species in x not in tree.") }
	if(max(times_in_tree > 1)){ stop("ERROR: some species in x in tree multiple times.")}

	# check that delim character isn't in tree tip names
	if(any(grepl(delim, x=tree$tip.label))){
		stop("ERROR: delim character is present in tip labels.")
	}

	# initialize output dist object
	output <- as.dist(matrix(nrow=length(x), ncol=length(x), dimnames=list(x,x), data=0))

	# make vector of which comparison types we need for output
	lt <- function(x){x[lower.tri(x)]}
	comp_pairs <- lt(outer(X=x, Y=x, FUN=function(x, y){paste(x, y, sep=delim)}))
	# make order (a-b vs b-a) irrelevant
	unorder <- function(x){
		xv <- unlist(strsplit(x, delim))
		xv <- xv[order(xv)]
		return(paste(xv, collapse=delim))
	}
	comp_pairs <- sapply(X=comp_pairs, FUN=unorder)

	# drop tips in tree not in x
	tree <- ape::keep.tip(tree, x)

	# make nested set object
	ns <- make_nested_set(tree)

	# wrapper function for bl_distance_ns that takes an item in comp_pairs as 1st arg
	d_ab <- function(ab, treei=tree, nsi=ns, delimi=delim){
		a <- unlist(strsplit(ab, delimi))[1]
		b <- unlist(strsplit(ab, delimi))[2]
		return(bl_distance_ns(tipa=a, tipb=b, treei, nsi))
	}

	# get branch-length distance for each unique pair
	u_comp_pairs <- unique(comp_pairs)
	if(n_cores <= 1){
		u_dists <- as.vector(sapply(X=u_comp_pairs, FUN=d_ab))
	}else{
		u_dists <- simplify2array(parallel::mclapply(X=u_comp_pairs, FUN=d_ab, mc.cores=n_cores))
	}

	# put u_dists into the right places in output
	for(i in 1:length(u_comp_pairs)){
		output[comp_pairs == u_comp_pairs[i]] <- u_dists[i]
	}

	# all done
	return(output)
}

