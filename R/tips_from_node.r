#' tips_from_node
#'
#' Determines which tip indices in a phylogeny descend from a given node. Called
#' by make_nested_set(), not intended for use otherwise, but some may find it handy.
#' Data should come from a rooted phylogeny, but this function doesn't check that so
#' be careful. 
#'
#' @author John L. Darcy
#'
#' @seealso 
#'   ape::phylo
#' 
#' @param nodes integer vector or scalar. The node index or indices for which
#'   tip indices are desired. 
#' @param anc integer vector. "ancestor" column vector from an adjacency matrix.
#'   For an ape::phylo object phy, anc=phy$edge[,1].
#' @param des integer vector. "descendant" column vector from an adjacency matrix.
#'   For an ape::phylo object phy, des=phy$edge[,2].
#' 
#' 
#' @return integer vector of tip indices, in no particular order.
#'
#' @examples
#'   library(specificity)
#'   library(ape)
#'   phy <- get(data(endophyte))$supertree
#'   # check if tree is rooted:
#'   ape::is.rooted(phy)
#'   # which tips are in the Cucurbitales?
#'   plot(phy) # need to stretch out the plot to see...
#'   nodelabels(adj=c(0,-1), bg="yellow") # node numbers
#'   nodelabels(phy$node.label, adj=c(0,1), bg="lightblue") # node names
#'   # we can see that Cucurbitales is node 107  
#'   cuc_tips <- tips_from_node( nodes=107, anc=phy$edge[,1], des=phy$edge[,2] )
#'   cuc_tips
#'   phy$tip.label[cuc_tips]
#'
#' @export
tips_from_node <- function(nodes, anc, des){
	# make sure no nodes are self-refferential! NOT ALLOWED!
	if(any(anc == des)){ stop("Bad phylogeny format - self refferential edges.") }
	# catch special case where input nodes contains root - return all tips!
	rootnode <- unique(anc[which(! anc %in% des)])
	if(rootnode %in% nodes){
		return(des[!des %in% anc]) # all tips
	}else{
		# while not all nodes are terminal, i.e. some nodes are in anc
		while( sum(nodes %in% anc) > 0 ){ 
			descs <- des[anc %in% nodes]
			terminals <- nodes[ ! nodes %in% anc]
			nodes <- c(descs, terminals)
		}
		return(nodes)
	}
}

