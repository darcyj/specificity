
#' phy_spec_simulate
#'
#' Simulates inputs for env_spec_multiple, by 
#' creating a species distribution over an artificial (or real) 
#' environmental variable. That distribution has a mean at the "ideal" 
#' environmental value for the simulated species, and the standard
#' deviation of that distribution controls the extent to which the species 
#' is specific to the variable. A high SD means less specificity, and a 
#' low SD means more specificity. 
#'
#' @author John L. Darcy
#'
#' @param sdev numeric vector. Standard deviation of the probability 
#'   distribution P(species|phy), in units of branch-length distance. 
#'   Low values mean that the species is found with a narrow phylogenetic 
#'   grouping of hosts, i.e. specificity. High values mean that the species 
#'   is found across a wider group of hosts, i.e. cosmopolitanism. Multiple 
#'   values can be input in order to simulate a range of specificities 
#'   simultaneously. To get a handle on this somewhat obscure variable, 
#'   consider plotting a histogram of cophenetic distances within hosts_phylo. 
#'   Can be length 1 or n.
#' @param ideal character vector. Tip label of hosts_phylo that is ideal (or
#'   closest to ideal) for the simulated species. Does not have to be in hosts.
#'   Can be length 1 or n. 
#' @param hosts character vector. Real of fake host identities. All must be 
#'   tips within hosts_phylo. 
#' @param hosts_phylo phylo object. Tree containing all unique hosts as tips.
#' @param n_obs integer vector. Number of positive observations to make, i.e. 
#'   occupancy of simulated species. Can be length 1 or n. 
#' @param sim_type character vector. "Fixed" means that values of hosts will
#'   be sampled using R's sample() function, resulting in exactly n_obs 
#'   observations. "Probabilistic" means that the total sum of probabilities 
#'   will be multiplied such that its sum is equal to n_obs, thus the most
#'   likely number of positive observations is n_obs. Whole word or first
#'   letter can be used ("p" or "f") ' (DEFAULT: "p").
#' @param n_cores integer. Number of CPU cores for parallel computation 
#'   (DEFAULT: 2).
#' @return List object containing "matrix" and "params" objects:
#'   \item{matrix}{
#'     matrix where each column is a vector of simulated boolean observations
#'     corresponding to a value of env; each row represents a simulated species.
#'     Matrix where a given column j is a vector of simulated boolean observations
#'     corresponding to the jth value of env. Each row 
#'   }
#'   \item{params}{
#'     data.frame of parameters (columns) used to simulate each species (rows).
#'   }
#'
#' @examples
#'   none yet written.
#'
#' @export
	phy_spec_simulate <- function(sdev, ideal, hosts, hosts_phylo, n_obs, sim_type="p", n_cores=2){
		require("ape")
		require("parallel")
		# make sure hosts is character and not factor
		hosts <- as.character(hosts)
		# deal with variable inputs by constructing table for each simulated species
		# each row is the 4 parameters for a species.
		var_lens <- c(length(sdev), length(ideal), length(n_obs), length(sim_type))
		if( ! all(var_lens) %in% c(1, max(var_lens))){
			stop("All input variables (sdev, ideal, n_obs, sim_type) must be either the length of the longest input variable or length 1.")
		}
		# mini 1-row df for each species that will be simulated.
		var_df <- data.frame(sdev=as.numeric(sdev), ideal=as.character(ideal), n_obs=as.integer(n_obs), sim_type=as.character(sim_type), stringsAsFactors=FALSE)
		var_list <- split(var_df, seq(nrow(var_df)))

		# function for calculating patristic distance between only two tips
		bl_distance <- function(tipa, tipb, tre){
			# just prune the tree so only tipa and tipb are in it, then sum branch lengths
			tre_simp <- ape::drop.tip(tre, tre$tip.label[! tre$tip.label %in% c(tipa, tipb)])
			return(sum(tre_simp$edge.length))
		}

		# function to apply for each simulated species
		# returns a logical vector describing whether species was observed for each
		# value of hosts
		phy_spec_sim_1species <- function(params){
			# vector of differences between each host in hosts and ideal
			# NOT calculated by taking a column out of a full cophenetic distance matrix
			# because some trees are much too large for that.
			dists2ideal <- sapply(X=hosts, FUN=function(x){ bl_distance(tipa=x, tipb=params$ideal, tre=hosts_phylo) } )
			# get probability of observing species with distance to ideal dists2ideal[i],
			# given that the species is observed exactly once.
			probs_1obs <- dnorm(x=dists2ideal, mean=0, sd=params$sdev) * 2

			# handle different ways of extrapolating 
			if(any(startsWith(params$sim_type, c("F", "f")))){
				# fixed n_obs, just use sample()
				# make vector of just falses for each host, trues added later
				output <- rep(FALSE, length(hosts))
				s_inds <- sample(1:length(hosts), size=params$n_obs, prob=probs_1obs)
				output[s_inds] <- TRUE
			}else if(any(startsWith(params$sim_type, c("P", "p")))){
				# probabilistic n_obs
				# adjust probs such that sum(probs) = n_obs
				probs_adj <- (probs_1obs/sum(probs_1obs)) * params$n_obs
				# clip any probs over 1 down to 1
				probs_adj[probs_adj > 1] <- 1
				# use weighted coin-flips to get observations
				# could also use sample(c(TRUE,FALSE), size=1, prob=c(p, 1-p))
				return(sapply(X=probs_adj, FUN=function(x){rbinom(n=1, size=1, prob=x)}))
			}else{
				stop(paste("ERROR: sim_type", params$sim_type, "not defined."))
			}
		}

		out_matrix <- simplify2array(mclapply(X=var_list, FUN=phy_spec_sim_1species, mc.cores=n_cores))
		return(list(
			matrix=out_matrix, 
			params=var_df
		))
	}












