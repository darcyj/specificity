#' env_spec_sim
#'
#' Simulates inputs for phy_or_env_spec, by creating a species distribution over
#' an artificial (or real) environmental variable. That distribution has a mean
#' at the "ideal" environmental value for the simulated species, and the standard
#' deviation of that distribution controls the extent to which the species is
#' specific to the variable. A high SD means weaker specificity, and a low SD 
#' means stronger specificity. 
#'
#' Since this process can result in failures (if a species is requested that's
#' highly specific to a region of env that isn't samples), some output species
#' will be failures. default operation is to remove those failures from output
#' matrix and output params data frame, but this can be changed.
#'
#' @author John L. Darcy
#'
#' @param sdev numeric vector. Standard deviation of the probability distribution
#'   P(species), in the same units of env. Low values mean that the species
#'   is found acrosss only a narrow range of env, i.e. specificity. High values 
#'   mean that the species is found across a wide range of env, i.e. 
#'   cosmopolitanism. Multiple values can be input in order to simulate a range
#'   of specificities simultaneously. Can be length 1 or n.
#' @param ideal numeric vector. Value of env that is ideal for the simulated
#'   species. This is the mode of the probability distribution P(species).
#'   Can be length 1 or n.
#' @param ideal2 numeric vector. Value of env that is the second ideal for the 
#'   simulated species. Only used if n_ideal >= 2. This is the second mode of the
#'   probability distribution P(species). Can be length 1 or n.
#' @param ideal3 numeric vector. Value of env that is the third ideal for the 
#'   simulated species. Only used if n_ideal = 3. This is the third mode of the
#'   probability distribution P(species). Can be length 1 or n.
#' @param n_ideal integer vector. Number of ideal values for the simulated species,
#'   i.e. modality of that species' distribution across env; 1 for unimodal, etc.
#'   Only can ue values 1, 2, or 3, which correspond to ideal, ideal2, and ideal3.
#'   Can be length 1 or n (default: 1).
#' @param env numeric vector. Real or fake environmental variable.
#' @param n_obs integer vector. Number of positive observations to make, i.e.
#'   occupancy of simulated species. Can be length 1 or n (default: 1).
#' @param up numeric vector. up=uniform proportion. This is the proportion of
#'   the probability distribution P(species) that is composed of a uniform
#'   distribution, if desired. If set to a value above zero (and blow 1), 
#'   P(species) will be a weighted sum of the normal distribution described above,
#'   and a uniform distribution. The weight for the uniform distribution will be
#'   up, and the weight for the normal distribution will be 1-up (default: 0).
#' @param oceanp numeric vector. oceanp=ocean proportion. This is the proportion of 
#'   samples in env that are "in the ocean", i.e. samples where the species would
#'   not expect to be found even if env is permissive. If aliens were calculating
#'   specificity of cows to temperature, they might look in the ocean at sites where
#'   the temperature is 17C (great for cows). But cows are not found in the ocean.
#'   This proportion is used to randomly select ocean sites within env, and then
#'   p(s|env|ocean) = up. Can be length 1 or n (default: 0).
#' @param n_cores integer. Number of CPU cores for parallel computation (default: 2).
#' @param seed integer. Seed for randomization. Daughter seeds will be generated for
#'   parallel computations, each with the same number of digits as seed 
#'   (default: 1234567).
#'
#' @return List object containing "matrix" and "params" objects:
#'   \describe{
#'     \item{matrix:}{
#'       matrix where each column is a vector of simulated observation frequencies
#'       (counts) corresponding to a value of env; each row represents a simulated 
#'       species.
#'     }
#'     \item{params:}{
#'       data.frame of parameters (columns) used to simulate each species (rows).
##'    }
#'   }
#' @examples
#'   none yet written.
#'
#' @export
env_spec_sim <- function(sdev, ideal, ideal2=0, ideal3=0, n_ideal=1, env,  n_obs, up=0, oceanp=0, n_cores=2, seed=1234567){
	require("parallel")
	# force n_obs to be integer
	n_obs <- round(n_obs)
	# deal with variable inputs by constructing table for each simulated species
	# each row is the 4 parameters for a species.
	var_lens <- c(length(sdev), length(ideal), length(ideal2), length(ideal3), length(n_ideal),
		length(up), length(oceanp), length(n_obs))
	if( ! all(var_lens) %in% c(1, max(var_lens))){
		stop("All input variables (sdev, ideal, n_obs) must be either the length of the longest input variable or length 1.")
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
	# mini 1-row df for each species that will be simulated. ="params" below.
	var_df <- data.frame(index=1:max(var_lens), sdev, ideal, ideal2, ideal3, n_ideal, n_obs, up, oceanp, seed=seeds)
	var_list <- lapply(X=as.list(1:nrow(var_df)), FUN=function(x){ var_df[x,] })

	# make sure all n_ideal values are in 1:3
	if(! all(var_df$n_ideal %in% 1:3)){
		stop("ERROR: Not all values of n_ideal were in {1,2,3}.")
	}

	prob_env_1species <- function(x, params){
		# this function is called from within env_spec_sim_1species()
		ideals <- c(params$ideal, params$ideal2, params$ideal3)[1:params$n_ideal]
		d2i <- sapply(X=x, FUN=function(x){min(abs(x-ideals))})
		# figure out which samples are in "ocean" (indices)
		oceansamps <- sample(1:length(x), round(params$oceanp * length(x)))
		# turn distances into probabilities
		unif_part <- rep(1/length(d2i), length(d2i)) * params$up
		norm_part <- dnorm(x=d2i, mean=0, sd=params$sdev)
		norm_part <- (norm_part / sum(norm_part)) * (1-params$up)
		output <- unif_part + norm_part
		# downgrade oceansamps to unif_part
		output[oceansamps] <- unif_part[oceansamps]
		return(output)
	}

	# function to apply for each simulated species
	# returns a logical vector describing whether species was observed for each
	# value of env.
	env_spec_sim_1species <- function(params){
		set.seed(params$seed)
		# get probability of observing species with env[i]
		# if that the species is observed exactly once.
		probs_1obs <- prob_env_1species(x=env, params)
		return(as.vector(rmultinom(n=1, size=params$n_obs, prob=probs_1obs)))
	}

	out_matrix <- simplify2array(mclapply(X=var_list, FUN=env_spec_sim_1species, mc.cores=n_cores))

	return(list(
		matrix=out_matrix, 
		params=var_df
	))
}