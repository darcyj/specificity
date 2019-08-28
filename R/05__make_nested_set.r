#' make_nested_set
#'
#' Makes a nested set table for a phylo object. Phylo objects made by the ape
#' package store phylogenies as an "adjacency  list", which in R is a table within
#' which any given edge is represented by the two node numbers it connects. With
#' this data structure, it is very computationally expensive to figure out which
#' tips are the descendents of a given node. Instead, using a "nested set" data 
#' structure, this operation is trivial. A nested set stores the minimum and 
#' maximum tip index for each node, such that the descendents of that node are
#' given by the inclusive range between those values. 
#'
#' @author John L. Darcy
#' @references
#'   https://en.wikipedia.org/wiki/Nested_set_model
#'   https://en.wikipedia.org/wiki/Adjacency_list
#'  
#' @seealso 
#'   ape::phylo
#'   geiger::tips
#'
#' @param phy phylo object. Must be rooted, and sorted such that tip indices are
#'   ordered. This is the default for rooted trees read in using ape's read.tree
#'   function.
#' @param n_cores integer. Number of CPU cores to use (DEFAULT: 2). lapply will
#'   be used instead of mclapply if ncores is 1.
#'
#' @return Matrix object representing a nested set of nodes. Each row matches
#'   rows of the "edges" object within phy. Object has the following columns:
#'   \describe{
#'     \item{1 (node)}{Node value in the original phylo object.}
#'     \item{2 (min)}{minimum tip index subtended by node.}
#'     \item{3 (max)}{maximum tip index subtended by node.}
#'     \item{4 (contig)}{Is min:max congiguous? 1 (true) or 0 (false).}
#'   }
#'
#' @examples
#'   library(geiger)
#'   library(ape)
#'   library(parallel)
#'   phy <- get(data(geospiza))$phy
#'   # check if tree is rooted:
#'   is.rooted(phy)
#'   # make nested set table:
#'   phy_ns <- make_nested_set(phy)
#'   # show that nested set table matches up with edges table in phy:
#'   all(phy$edge[,2] == phy_ns[,1])
#'
#' @export
make_nested_set <- function(phy, n_cores=2){
	require(geiger)		
	require(ape)
	require(parallel)
	if( ! ape::is.rooted(phy)){
		stop("ERROR: phylogeny is not rooted. Maybe try a midpoint root?")
	}
	# this function is just a wrapper for geiger's tips function, so that
	# it returns only the node, min, max, and errors
	get_node_tip_range <- function(node){
		# find indices of descendent tips
		descs <- which(phy$tip.label %in% geiger::tips(phy=phy, node=node))
		min <- min(descs)
		max <- max(descs)
		# check to make sure descs are contiguous throughout range
		# this is an assumption of the nested set approach.
		# initialize contig ("conriguous") to fail state:
		contig <- FALSE
		# if tips are contiguous, set contig to true.
		if( all(descs %in% min:max)){ contig <- TRUE }
		# return a data.frame object
		return(data.frame(node=node, min=min, max=max, contig=contig))
	}
	# parallelize get_node_tip_range across each node
	if(n_cores <= 1){
		dfs <- lapply(X=phy$edge[,2], FUN=get_node_tip_range)
	}else{
		dfs <- mclapply(X=phy$edge[,2], FUN=get_node_tip_range, mc.cores=n_cores)
	}
	# combine data frames. if data.table is available, use that since it's FAST
	# otherwise just use do.call
	if(requireNamespace("data.table", quietly=TRUE)){
		dfs <- as.data.frame(data.table::rbindlist(dfs))
	}else{
		dfs <- do.call("rbind", dfs)
	}
	# if there was an contiguity error, do a warning
	if( ! all(dfs$contig) ){
		warning("One or more nodes did not contain contiguous tip ranges.")
	}

	return(as.matrix(dfs))
}

#' bl_distance_ns
#'
#' Calculates branch-length distance between tipa and tipb in a phylogenetic tree
#' using nested-set optomization. Requires a pre-calculated nested-set.
#'
#' @author John L. Darcy
#'
#' @param tipa string. Name of a tip in tree.
#' @param tipb string. Name of another tip in tree.
#' @param tree phylo object. Tree containing all unique species in x as tips. May
#'   contain tips that are not in x.
#' @param ns matrix. Nested-set matrix for tree; use make_nested_set(tree).
#' 
#' @return Distance between tipa and tipb.
#'
#' @examples
#'   library(ape)
#'   example_tree <- ape::read.tree(text=" ((((a:1,b:1):1,c:2):1,d:3):1,(e:1,f:1):3);")
#'   plot(example_tree); axis(side=1)
#'   example_ns <- make_nested_set(example_tree)
#'   bl_distance_ns("a", "c", example_tree, example_ns) # should be 4
#'   bl_distance_ns("a", "f", example_tree, example_ns) # should be 8
#'   bl_distance_ns("d", "c", example_tree, example_ns) # should be 6
#'   
#' 
#' @export
bl_distance_ns <- function(tipa, tipb, tree, ns){
	# numbers for tipa and tipb
	nab <- which(tree$tip.label %in% c(tipa, tipb))
	if(tipa == tipb){
		return(0)
	}else if(length(nab) < 2){
		stop("ERROR: Either tipa or tipb is not in tree.")
	}else if(length(nab) > 2){
		stop("ERROR: tipa or tipb are in tree multiple times.")
	}
	# let b be true only for nodes that descend in tipa or tipb, but not both.
	# ns[,3] is maxes, ns[,2] is mins.
	b <- xor((ns[,3] >= nab[1] & ns[,2] <= nab[1]), (ns[,3] >= nab[2] & ns[,2] <= nab[2]))
	# sum up edges immediately descendent from nodes where b is true
	return(sum(tree$edge.length[b]))
}
