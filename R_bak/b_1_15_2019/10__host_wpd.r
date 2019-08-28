#' host_wpd
#'
#' Calculates Weighted Phylogenetic Diversity (WPD) of host species, given 
#' occurrences of some species i that is hosted by those species. 
#' WPD is not weighted by the abundances of species
#' i, but rather by the number of times species i was observed for a given host.
#' For example, if the vector of abundances for species i is {3, 40, 0, 12} and
#' the vector of host identities is {a, a, b, b}, then the weights will be 2 for
#' a and 1 for b. This function can randomly permute host identities so that a 
#' null WPD can be calculated. Note that all data must be included for this null
#' option, i.e. observations where species i was NOT found (0 abundance, see 
#' example above) must be included. This function is adapted from a function in
#' the package lefse, which can be found at 'github.com/NGSwenson/lefse_0.5'. 
#'
#' @author John L. Darcy
#' @references
#' Swenson NG (2014) Functional and Phylogenetic Ecology in R. 
#'   Springer UseR! Series, Springer, New York, New York, U.S.A.
#'  
#' @seealso github.com/NGSwenson/lefse_0.5
#'
#' @param abunds numeric vector. Abundances for species i.
#' @param hosts character vector. Host identities corresponding to abunds.
#' @param hosts_phylo phylo object. Tree containing all names in hosts as tips.
#' @param null logical. If TRUE, hosts will be randomly permuted (DEFAULT: FALSE).
#' @param null_occ integer. If > 0, null model will be run at specified
#'   occupancy, i.e. number of observations will be forced to null_occ. If given,
#'   abunds argument is arbitrary but must still be provided (DEFAULT: 0).
#' @param unweighted logical. If TRUE, will use PD instead of WPD (DEFAULT: FALSE).
#'
#' @return Single WPD or PD value.
#'
#' @examples
#' # none yet written.
#'
#' @export
	host_wpd <- function(abunds, hosts, hosts_phylo, null=FALSE, null_occ=0, unweighted=FALSE){
		devtools::use_package("ape")
		devtools::use_package("geiger")
		# if null model is desired, scramble host identities.
		if(null==TRUE && null_occ > 0){ 
			abunds <- sample(c(rep(1, null_occ), rep(0, length(hosts) - null_occ)))
		}else if(null==TRUE && null_occ <= 0){
			abunds <- sample(abunds)
		}
		# tabulate host IDs where abunds > 0, set col names for cdoe readability
		hosts_table <- data.frame(table(hosts[abunds > 0]))
		colnames(hosts_table) <- c("host", "count")
		# prune tree to only contain tips in hosts_table
		hosts_phylo <- ape::drop.tip(hosts_phylo, hosts_phylo$tip.label[! hosts_phylo$tip.label %in% hosts_table$host])
		if(unweighted==TRUE){
			# regular PD is SO SIMPLE to calculate - this works because tree is pruned.
			return(sum(hosts_phylo$edge.length))
		}else{
			# make edges data frame, with col names for readability
			# note: this function could be simplified by omitting this df, but I think
			#   that it is more legible this way.
			edges <- as.data.frame(matrix(0, nrow=nrow(hosts_phylo$edge), ncol=3))
			colnames(edges) <- c("node", "edge_length", "mean_abund")
			# Fill first two columns with node numbers defining each edge
			edges$node <- hosts_phylo$edge[,2]
			# Fill the third column with the length of that edge
			edges$edge_length <- hosts_phylo$edge.length
			# for each edge, calculate the average abundance of its descendent tips
			node_mean_abund <- function(x){
				tip_names <- geiger::tips(phy=hosts_phylo, node=x)
				return(mean(hosts_table$count[ hosts_table$host %in% tip_names ]))
			}
			edges$mean_abund <- sapply(X=edges$node, FUN=node_mean_abund)
			# calculate Weighted Phylogenetic Diversity
			return(nrow(edges) * ((sum(edges$edge_length * edges$mean_abund)) / sum(edges$mean_abund)))
		}
	}
