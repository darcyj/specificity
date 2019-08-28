
#' env_spec_sim
#'
#' Simulates inputs for phy_or_env_spec, by creating a species distribution over
#' an artificial (or real) environmental variable. That distribution has a mean
#' at the "ideal" environmental value for the simulated species, and the standard
#' deviation of that distribution controls the extent to which the species is
#' specific to the variable. A high SD means less specificity, and a low SD means
#' more specificity. 
#'
#' Since this process can result in failures (if a species is requested that's
#' highly specific to a region of env that isn't samples), some output species
#' will be failures. Default operation is to remove those failures from output
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
#'   species. This is the mean of the probability distribution P(species).
#'   This value does NOT have to actually be within env, or even within the 
#'   range of env. Can be length 1 or n.
#' @param env numeric vector. Real or fake environmental variable.
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
#'   can happen if env has a large gap (no values) and you specify an ideal and sd
#'   such that p(s|env) is zero for all input values of env. Output objects will be
#'   SMALLER than inputs if TRUE. If FALSE, some columns of matrix output will be
#'   filled with NAs, but number of output species will match inputs (DEFAULT: TRUE).
#' @param n_cores integer. Number of CPU cores for parallel computation (DEFAULT: 2).
#' @param seed integer. Seed for randomization. Daughter seeds will be generated for
#'   parallel computations, each with the same number of digits as seed 
#'   (DEFAULT: 1234567).
#'
#' @return List object containing "matrix" and "params" objects:
#'   \item{matrix}{
#'     matrix where each column is a vector of simulated boolean observations
#'     corresponding to a value of env; each row represents a simulated species.
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
	env_spec_sim <- function(sdev, ideal, env,  n_obs, up=0, sim_type="wtd_f", fail_rm=TRUE, n_cores=2, seed=1234567){
		require("parallel")
		# force n_obs to be integer
		n_obs <- round(n_obs)
		# deal with variable inputs by constructing table for each simulated species
		# each row is the 4 parameters for a species.
		var_lens <- c(length(sdev), length(ideal), length(up), length(n_obs), length(sim_type))
		if( ! all(var_lens) %in% c(1, max(var_lens))){
			stop("All input variables (sdev, ideal, n_obs, sim_type) must be either the length of the longest input variable or length 1.")
		}
		# generate n_sim daughter seeds
		set.seed(seed)
		seeds <- replicate(n=max(var_lens), as.integer(paste(sample(0:9, nchar(seed), replace=TRUE), collapse="")))
		seeds <- formatC(seeds, width=nchar(seed), format="d", flag="0")
		# mini 1-row df for each species that will be simulated.
		var_df <- data.frame(index=1:max(var_lens), sdev, ideal, n_obs, up, seed=seeds, sim_type)
		var_list <- lapply(X=as.list(1:nrow(var_df)), FUN=function(x){ var_df[x,] })

		prob_env_1species <- function(x, sdev, ideal, up){
			mean <- ideal
			unif_part <- dunif(x=x, min=min(x), max=max(x))
			norm_part <- dnorm(x=x, mean=ideal, sd=sdev)
			# make sure each part sums to 1, and multiply by weights. then add.
			unif_part <- (unif_part / sum(unif_part)) * up
			norm_part <- (norm_part / sum(norm_part)) * (1 - up)
			# combine the two using add
			return(unif_part + norm_part)
		}

		# function to apply for each simulated species
		# returns a logical vector describing whether species was observed for each
		# value of env.
		env_spec_sim_1species <- function(params){
			set.seed(params$seed)
			# get probability of observing species with env[i]
			# if that the species is observed exactly once.
			probs_1obs <- prob_env_1species(x=env, sdev=params$sdev, ideal=params$ideal, up=params$up)

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

		out_matrix <- simplify2array(mclapply(X=var_list, FUN=env_spec_sim_1species, mc.cores=n_cores))

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


