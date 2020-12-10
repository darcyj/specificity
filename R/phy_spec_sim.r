
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
#'   closest to ideal) for the simulated species. Does not have to be in hosts,
#'   but MUST be in hosts_phylo. Can be length 1 or n. 
#' @param ideal2 character vector. Tip label of hosts_phylo that is secondary
#'   ideal host for the simulated species. Does not have to be in hosts, but MUST
#'   be in hosts_phylo. Can be blank ("") if corresponding n_ideal < 2. Can be
#'   length 1 or n (default: "").
#' @param ideal3 character vector. Tip label of hosts_phylo that is tertiary
#'   ideal host for the simulated species. Does not have to be in hosts, but MUST
#'   be in hosts_phylo. Can be blank ("") if corresponding n_ideal < 3. Can be
#'   length 1 or n (default: "").
#' @param n_ideal integer vector. number of ideal hosts to use. Must be 1, 2, or
#'   3 (default: 1).
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
#' @param oceanp numeric vector. See ?env_spec_sim for help.
#' @param n_cores integer. Number of CPU cores for parallel computation (default: 2).
#' @param seed integer. Seed for randomization. Daughter seeds will be generated for
#'   parallel computations, each with the same number of digits as seed 
#'   (default: 1234567).
#'
#' @return List object containing "matrix" and "params" objects:
#'   \describe{
#'     \item{matrix:}{
#'       matrix where each column is a vector of simulated observations corresponding
#'       to a value of hosts; each row represents a simulated species.
#'     }
#'     \item{params:}{
#'       data.frame of parameters (columns) used to simulate each species (rows).
#'       A column called "index" is included so that simulated species can be mapped
#'       back onto original data structures when some species are ommitted due to 
#'       simulation failure (see fail_rm).
#'     }
#'   }
#'
#' @examples
#'   none yet written.
#'
#' @export
phy_spec_sim <- function(sdev, ideal, ideal2="", ideal3="", n_ideal=1, hosts, hosts_phylo, 
	n_obs, up=0, oceanp=0, n_cores=2, seed=1234567){

	require("parallel")
	# make sure hosts is character and not factor
	hosts <- as.character(hosts)
	# deal with variable inputs by constructing table for each simulated species
	# each row is the 4 parameters for a species.
	var_lens <- c(length(sdev), length(ideal), length(n_obs), length(up))
	if( ! all(var_lens) %in% c(1, max(var_lens))){
		stop("All input variables (sdev, ideals, n_ideal, n_obs, up, oceanp) must be either
			the length of the longest input variable or length 1.")
	}
	# check to make sure certain variables are within expected ranges:
	checkrange <- function(x, r=c(0,1)){
		if(min(x) < r[1] || max(x) > r[2]){	return(FALSE) }else{ return(TRUE) }
	}
	if(!checkrange(up)){stop("Range of up not in [0,1]")}
	if(!checkrange(oceanp)){stop("Range of oceanp not in [0,1]")}


	# generate n_sim daughter seeds
	set.seed(seed)
	seeds <- replicate(n=max(var_lens), as.integer(paste(sample(0:9, nchar(seed), replace=TRUE), collapse="")))
	seeds <- formatC(seeds, width=nchar(seed), format="d", flag="0")
	# mini 1-row df for each species that will be simulated.
	var_df <- data.frame(index=1:max(var_lens), sdev, ideal, ideal2, ideal3, n_ideal, n_obs, up, 
		oceanp, seed=seeds, stringsAsFactors=F)
	var_list <- split(var_df, seq(nrow(var_df)))

	# validate ideal, ideal2, ideal3
	if(! all(var_df$n_ideal %in% 1:3)){
		stop("Not all values of n_ideal were in {1,2,3}.")
	}
	ideal_v <- data.frame(
		i1v = var_df$ideal  %in% hosts_phylo$tip.label, # all must be true for validity
		i2v = var_df$ideal2 %in% hosts_phylo$tip.label, # can be FALSE if n_ideal < 2
		i3v = var_df$ideal3 %in% hosts_phylo$tip.label # can be FALSE if n_ideal < 3
	)
	ideal_v$i2v[var_df$n_ideal < 2] <- TRUE # OK to have invalid ideal2 if n_ideal < 2
	ideal_v$i3v[var_df$n_ideal < 3] <- TRUE # OK to have invalid ideal3 if n_ideal < 3
	badspecies <- apply(X=!ideal_v, MAR=1, FUN=any)
	if(any(badspecies)){
		stop("One or more species had invalid ideal/ideal2/ideal3 values.")
	}

	# make nested set for hosts_phylo
	nested_set <- make_nested_set(phy=hosts_phylo, n_cores=n_cores)

	# d = distances from each tip in tree to ideal tip
	prob_phy_1species <- function(d, params){
		oceansamps <- sample(1:length(d), round(params$oceanp * length(d)))
		unif_part <- dunif(x=d, min=min(d), max=max(d))
		norm_part <- dnorm(x=d, mean=0, sd=params$sdev)
		# make sure each part sums to 1, and multiply by weights. then add.
		unif_part <- (unif_part / sum(unif_part)) * params$up
		norm_part <- (norm_part / sum(norm_part)) * (1 - params$up)
		# combine the two using add
		output <- unif_part + norm_part
		# downgrade oceansamps to unif_part
		output[oceansamps] <- unif_part[oceansamps]
		return(output)
	}

	# function to apply for each simulated species
	phy_spec_sim_1species <- function(params){
		# vector of differences between each host in hosts and ideal
		# NOT calculated by taking a column out of a full cophenetic distance matrix
		# because some trees are much too large for that.
		dists2ideal <- sapply(X=hosts, FUN=function(x){ bl_distance_ns(tipa=x, tipb=params$ideal, tree=hosts_phylo, ns=nested_set) } )
		if(params$n_ideal == 2){
			dists2ideal2 <- sapply(X=hosts, FUN=function(x){ bl_distance_ns(tipa=x, tipb=params$ideal2, tree=hosts_phylo, ns=nested_set) } )
			dists2ideal <- apply(X=cbind(dists2ideal, dists2ideal2), MAR=1, FUN=min)
		}else if(params$n_ideal == 3){
			dists2ideal2 <- sapply(X=hosts, FUN=function(x){ bl_distance_ns(tipa=x, tipb=params$ideal2, tree=hosts_phylo, ns=nested_set) } )
			dists2ideal3 <- sapply(X=hosts, FUN=function(x){ bl_distance_ns(tipa=x, tipb=params$ideal3, tree=hosts_phylo, ns=nested_set) } )
			dists2ideal <- apply(X=cbind(dists2ideal, dists2ideal2, dists2ideal3), MAR=1, FUN=min)
		}
		# get probability of observing species with distance to ideal dists2ideal[i],
		# given that the species is observed exactly once.
		probs_1obs <- prob_phy_1species(d=dists2ideal, params)
		return(tryCatch({
			return(as.vector(rmultinom(n=1, size=params$n_obs, prob=probs_1obs)))
		}, error = function(e){
			rep(NA, length(probs_1obs))
		}))
		#return(as.vector(rmultinom(n=1, size=params$n_obs, prob=probs_1obs)))
	}

	out_matrix <- simplify2array(mclapply(X=var_list, FUN=phy_spec_sim_1species, mc.cores=n_cores))

	return(list(
		matrix=out_matrix, 
		params=var_df
	))
}

