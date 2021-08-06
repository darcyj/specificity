
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
#' @return dist object, of vector length equal to (l^2-l)/2 where l is length(x); i.e. 
#'   values are the lower triangle of a patristic distance matrix with rows=x and cols=x. 
#'
#' @examples
#' # library(specificity)
#' # library(ape)
#' # example_tree <- ape::read.tree(text=" ((((a:1,b:1):1,c:2):1,d:3):1,(e:1,f:1):3);")
#' # example_x <- c("a", "a", "a", "b", "c", "d", "c", "a", "f")
#' # # unique patristic distance matrix:
#' # ape::cophenetic.phylo(example_tree)
#' # # dist object for example_x:
#' # tree2mat(tree=example_tree, x=example_x)
#' 
#' @export
tree2mat <- function(tree, x, n_cores=1 ){
	# check that each species in x is in tree exactly once
	times_in_tree <- sapply(x, FUN=function(x){sum(tree$tip.label %in% x)})
	if(0 %in% times_in_tree){ stop("ERROR: some species in x not in tree.") }
	if(max(times_in_tree > 1)){ stop("ERROR: some species in x in tree multiple times.")}

	# make vector of which comparison types we need for output
	lt <- function(x){x[lower.tri(x)]}

	# initialize output dist object
	output <- as.dist(matrix(nrow=length(x), ncol=length(x), dimnames=list(x,x), data=0))

	# two columns of names: a and b
	# same order as output lower-triangle dist output
	comps2do <- data.frame(
		a=lt(outer(X=x, Y=x, FUN=function(a, b){a})),
		b=lt(outer(X=x, Y=x, FUN=function(a, b){b})),
		stringsAsFactors=FALSE
	)

	# order so that names a and b are in alphabetical order, without messing up
	# ordering per dist output
	comps2do <- data.frame(
		a=pmin(comps2do$a, comps2do$b), 
		b=pmax(comps2do$a, comps2do$b),
		stringsAsFactors=FALSE
	)

	# unique comparisons. same data structure as above but unique.
	# as.matrix is required because of funkines with rows of data.frames being
	# data.frames themselves. 
	uniquecomps <- as.matrix(unique(comps2do))

	# transform uniquecomps to list, to enable lapply functions without messy
	# index-as-argument approach
	uniquecomps <- lapply(X=1:nrow(uniquecomps), FUN=function(i){uniquecomps[i,]} )

	# drop tips in tree not in x
	tree <- ape::keep.tip(tree, x)

	# make nested set object
	ns <- make_nested_set(tree)

	# wrapper function for bl_distance_ns that takes an item in uniquecomps as 1st arg
	d_ab <- function(ab, treei=tree, nsi=ns){
		return(bl_distance_ns(tipa=ab[1], tipb=ab[2], treei, nsi))
	}

	# get branch-length distance for each unique comparison
	if(n_cores <= 1){
		uniquedists <- as.vector(sapply(X=uniquecomps, FUN=d_ab))
	}else{
		uniquedists <- simplify2array(parallel::mclapply(X=uniquecomps, FUN=d_ab, mc.cores=n_cores))
	}

	# put values where they belong
	for(i in 1:length(uniquecomps)){
		ab <- uniquecomps[[i]]
		d <- uniquedists[[i]]
		if(ab[1] != ab[2]){
			output[(comps2do$a == ab[1]) & (comps2do$b == ab[2])] <- d
		}
	}

	# all done
	return(output)
}

