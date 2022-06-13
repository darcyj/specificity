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
#'
#' @param phy phylo object. Must be rooted, and sorted such that tip indices are
#'   ordered. This is the default for rooted trees read in using ape's read.tree
#'   function.
#' @param n_cores integer. Number of CPU cores to use for parallel operations. If
#'   set to 1, lapply will be used instead of mclapply. A warning will be shown if
#'   n_cores > 1 on Windows, which does not support forked parallelism (default: 2).
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
#'   # library(specificity)
#'   # library(ape)
#'   # phy <- get(data(endophyte))$supertree
#'   # # check if tree is rooted:
#'   # ape::is.rooted(phy)
#'   # # make nested set table:
#'   # phy_ns <- make_nested_set(phy)
#'   # # show that nested set table matches up with edges table in phy:
#'   # all(phy$edge[,2] == phy_ns[,1])
#'
#' @export
make_nested_set <- function(phy, n_cores=2){
	# check if tree is rooted
	if( ! ape::is.rooted(phy)){
		stop("ERROR: phylogeny is not rooted. Maybe try a midpoint root?")
	}
    # warn if ncores > 1 and platform isn't  "unix" (windows can't do forked parallelism)
    if(n_cores > 1 && .Platform$OS.type != "unix"){
        warning("Windows is incompatible with n_cores > 1.")
    }
	# this function is just a wrapper for tips_from_node, so that
	# it returns only the node, min, max, and errors
	.get_node_tip_range <- function(node){
		# find indices of descendent tips
		descs <- tips_from_node(node, anc=phy$edge[,1], des=phy$edge[,2])
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
		dfs <- lapply(X=phy$edge[,2], FUN=.get_node_tip_range)
	}else{
		dfs <- parallel::mclapply(X=phy$edge[,2], FUN=.get_node_tip_range, mc.cores=n_cores)
	}
	# combine data frames.
	dfs <- do.call("rbind", dfs)
	# if there was an contiguity error, do a warning
	if( ! all(dfs$contig) ){
		warning("One or more nodes did not contain contiguous tip ranges.")
	}

	return(as.matrix(dfs))
}



