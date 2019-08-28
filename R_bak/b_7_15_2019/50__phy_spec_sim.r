
#' phy_spec_sim
#'
#' Simulates inputs for phy_or_env_spec, by creating a species distribution over
#' an artificial (or real) host phylogenetic tree. For a phylogeny, the species
#' probability distribution P(s) is based on patristic distances within the
#' tree, such that P(s) is maximized at zero patristic distance between a tip
#' in the tree and the ideal host species for s. This distribution is given by a
#' a truncated normal distribution centered on zero, using only positive values.
#' A uniform proportion (up) to that distribution may be added as well, to add a
#' baseline probability to P(s). The standard deviation of P(s) can be raised or
#' lowered to simulate cosmopolitanism or specificity. 
#'
#' @author John L. Darcy
#'
#' @param sdev numeric vector. Standard deviation of the probability distribution
#'   P(s), in units of patristic distance in hosts_phylo. Low values  mean that
#'   species s is found with a narrow grouping of hosts, i.e. specificity. High
#'   values mean that s is found across a wider group of hosts, i.e.
#'   cosmopolitanism. Multiple values can be input in order to simulate a range
#'   of specificities, simultaneously. To get a handle on this somewhat opaque 
#'   variable, consider plotting a histogram of patristic distances within 
#'   hosts_phylo (see: ape::cophenetic.phylo). Can be length 1 or n.
#' @param ideal character vector. Tip label of hosts_phylo that is ideal (or
#'   closest to ideal) for the simulated species. Does not have to be in hosts.
#'   Can be length 1 or n. 
#' @param hosts character vector. Real of fake host identities. All must be 
#'   tips within hosts_phylo. Analogous to env argument to env_spec_sim.
#' @param hosts_phylo phylo object. Tree containing all unique hosts as tips.
#' @param n_obs integer vector. Number of positive observations to make, i.e.
#'   occupancy of simulated species. Can be length 1 or n.
#' @param up numeric vector. up=uniform proportion. This is the proportion of
#'   the probability distribution P(species) that is composed of a uniform
#'   distribution, if desired. If set to a value above zero (and blow 1), 
#'   P(species) will be a weighted sum of the normal distribution described above,
#'   and a uniform distribution. The weight for the uniform distribution will be
#'   up, and the weight for the normal distribution will be 1-up (default: 0).
#' @param sim_type character vector. Determines what type of data are created.
#'   Weighted data are integer values <= 0, presence-absence data are boolean.
#'   \item{unw_f}{
#'     Presence-absence simulation, with a fixed number of observations.
#'   }
#'   \item{unw_p}{
#'     Presence-absence simulation, with a probabilistic number of observations.
#'   }
#'   \item{wtd_f}{
#'     Weighted simulation, with a fixed number of observations (DEFAULT).
#'   }
#' @param fail_rm logical. Should failed species simulations be removed T/F? This
#'   can happen if hosts does not contain many species that are phylogenetically
#'   close to ideal. Output objects will be SMALLER than inputs if TRUE. If FALSE,
#'   some columns of matrix output will be filled with NAs, but number of output
#'   species will match inputs.
#' @param n_cores integer. Number of CPU cores for parallel computation (DEFAULT: 2).
#' @param seed integer. Seed for randomization. Daughter seeds will be generated for
#'   parallel computations, each with the same number of digits as seed 
#'   (DEFAULT: 1234567).
#'
#' @return List object containing "matrix" and "params" objects:
#'   \item{matrix}{
#'     matrix where each column is a vector of simulated observations corresponding
#'     to a value of hosts; each row represents a simulated species.
#'   }
#'   \item{params}{
#'     data.frame of parameters (columns) used to simulate each species (rows).
#'     A column called "index" is included so that simulated species can be mapped
#'     back onto original data structures when some species are ommitted due to 
#'     simulation failure (see fail_rm).
#'   }
#'
#' @examples
#'   none yet written.
#'
#' @export
	prob_phy_1species <- function(d, sdev, up){
		unif_part <- dunif(x=d, min=min(d), max=max(d))
		norm_part <- dnorm(x=d, mean=0, sd=sdev)
		# make sure each part sums to 1, and multiply by weights. then add.
		unif_part <- (unif_part / sum(unif_part)) * up
		norm_part <- (norm_part / sum(norm_part)) * (1 - up)
		# combine the two using add
		return(unif_part + norm_part)
	}

	phy_spec_sim <- function(sdev, ideal, hosts, hosts_phylo, n_obs, up, sim_type="wtd_f", 
		fail_rm=TRUE, n_cores=2, seed=1234567){

		require("parallel")
		# make sure hosts is character and not factor
		hosts <- as.character(hosts)
		# deal with variable inputs by constructing table for each simulated species
		# each row is the 4 parameters for a species.
		var_lens <- c(length(sdev), length(ideal), length(n_obs), length(up), length(sim_type))
		if( ! all(var_lens) %in% c(1, max(var_lens))){
			stop("All input variables (sdev, ideal, n_obs, sim_type) must be either the length of the longest input variable or length 1.")
		}
		# generate n_sim daughter seeds
		set.seed(seed)
		seeds <- replicate(n=max(var_lens), as.integer(paste(sample(0:9, nchar(seed), replace=TRUE), collapse="")))
		seeds <- formatC(seeds, width=nchar(seed), format="d", flag="0")
		# mini 1-row df for each species that will be simulated.
		var_df <- data.frame(index=1:max(var_lens), sdev, ideal, n_obs, up, seed=seeds, sim_type, stringsAsFactors=F)
		var_list <- split(var_df, seq(nrow(var_df)))

		# make nested set for hosts_phylo
		nested_set <- make_nested_set(phy=hosts_phylo, n_cores=n_cores)

		# function for calculating patristic distance between only two tips
		# requires nested set because that makes it about 12x faster
		bl_distance_ns <- function(tipa, tipb, tre=hosts_phylo, ns=nested_set){
			# numbers for tipa and tipb
			nab <- which(tre$tip.label %in% c(tipa, tipb))
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
			return(sum(tre$edge.length[b]))
		}

		# function to apply for each simulated species
		phy_spec_sim_1species <- function(params){
			# vector of differences between each host in hosts and ideal
			# NOT calculated by taking a column out of a full cophenetic distance matrix
			# because some trees are much too large for that.
			dists2ideal <- sapply(X=hosts, FUN=function(x){ bl_distance_ns(tipa=x, tipb=params$ideal) } )
			# get probability of observing species with distance to ideal dists2ideal[i],
			# given that the species is observed exactly once.
			probs_1obs <- prob_phy_1species(d=dists2ideal, sdev=params$sdev, up=params$up)


			# handle different ways of extrapolating 
			if(params$sim_type == "unw_f"){
				# check that there are enough positive probabilities to even simulate a species.
				# if not, return NAs instead.
				if(sum(probs_1obs > 0) < params$n_obs ){
					return(rep(NA, length(env)))
				}else{
					# fixed n_obs, just use sample()
					# make vector of just falses for each ENV, trues added later
					output <- rep(FALSE, length(env))
					s_inds <- sample(1:length(env), size=params$n_obs, prob=probs_1obs)
					output[s_inds] <- TRUE
					return(output)
				}
			}else if(params$sim_type == "unw_p"){
				# check that there are enough positive probabilities to even simulate a species.
				# if not, return NAs instead.
				if(sum(probs_1obs > 0) < params$n_obs ){
					return(rep(NA, length(env)))
				}else{
					# probabilistic n_obs
					# adjust probs such that sum(probs) = n_obs
					probs_adj <- (probs_1obs/sum(probs_1obs)) * params$n_obs
					# clip any probs over 1 down to 1
					probs_adj[probs_adj > 1] <- 1
					# use weighted coin-flips to get observations
					# could also use sample(c(TRUE,FALSE), size=1, prob=c(p, 1-p))
					return(sapply(X=probs_adj, FUN=function(x){rbinom(n=1, size=1, prob=x)}))
				}
			}else if(params$sim_type == "wtd_f"){
				# weighted n_obs, fixed.
				return(as.vector(rmultinom(n=1, size=params$n_obs, prob=probs_1obs)))
			}else{
				stop(paste("ERROR: sim_type", params$sim_type, "not defined."))
			}
		}

		out_matrix <- simplify2array(mclapply(X=var_list, FUN=phy_spec_sim_1species, mc.cores=n_cores))

		# remove failures if requested
		if(fail_rm){
			goods <- apply(X=out_matrix, MAR=2, FUN=function(x){!all(is.na(x))})
			out_matrix <- out_matrix[,goods]
			var_df <- var_df[goods,]
		}

		return(list(
			matrix=out_matrix, 
			params=var_df
		))
	}

