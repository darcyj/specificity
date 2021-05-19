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
#' #  library(specificity)
#' #  library(ape)
#' #  example_tree <- ape::read.tree(text=" ((((a:1,b:1):1,c:2):1,d:3):1,(e:1,f:1):3);")
#' #  plot(example_tree); axis(side=1)
#' #  example_ns <- make_nested_set(example_tree)
#' #  bl_distance_ns("a", "c", example_tree, example_ns) # should be 4
#' #  bl_distance_ns("a", "f", example_tree, example_ns) # should be 8
#' #  bl_distance_ns("d", "c", example_tree, example_ns) # should be 6
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
