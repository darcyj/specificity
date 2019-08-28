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
#'   \item{node}{Node value in the original phylo object.}
#'   \item{min}{minimum tip index subtended by node.}
#'   \item{max}{maximum tip index subtended by node.}
#'   \item{contig}{Is min:max congiguous? T(1) or F(0).}
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
#'   all(phy$edge[,2] == phy_ns$node)
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
	# if there was an conriguity error, do a warning
	if( ! all(dfs$contig) ){
		warning("One or more nodes did not contain contiguous tip ranges.")
	}

	return(as.matrix(dfs))
}

