
#' env_spec_simulate
#'
#' Simulates inputs for env_spec_multiple, by creating a species distribution over
#' an artificial (or real) environmental variable. That distribution has a mean
#' at the "ideal" environmental value for the simulated species, and the standard
#' deviation of that distribution controls the extent to which the species is
#' specific to the variable. A high SD means less specificity, and a low SD means
#' more specificity. 
#'
#' @author John L. Darcy
#'
#' @param sdev numeric vector. Standard deviation of the probability distribution
#'   P(species|env), in the same units of env. Low values mean that the species
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
#' @param sim_type character vector. "Fixed" means that values of env will be
#'   sampled using R's sample() function, resulting in exactly n_obs observations.
#'   "Probabilistic" means that the total sum of probabilities will be multiplied
#'   such that its sum is equal to n_obs, thus the most likely number of positive
#'   observations is n_obs. Whole word or first letter can be used ("p" or "f")
#'   (DEFAULT: "p").
#' @param n_cores integer. Number of CPU cores for parallel computation (DEFAULT: 2).
#' @return List object containing "matrix" and "params" objects:
#'   \item{matrix}{
#'     matrix where each column is a vector of simulated boolean observations
#'     corresponding to a value of env; each row represents a simulated species.
#'   }
#'   \item{params}{
#'     data.frame of parameters (columns) used to simulate each species (rows).
#'   }
#'
#' @examples
#'   none yet written.
#'
#' @export
	env_spec_simulate <- function(sdev, ideal, env, n_obs, sim_type="p", n_cores=2){
		require("parallel")
		# deal with variable inputs by constructing table for each simulated species
		# each row is the 4 parameters for a species.
		var_lens <- c(length(sdev), length(ideal), length(n_obs), length(sim_type))
		if( ! all(var_lens) %in% c(1, max(var_lens))){
			stop("All input variables (sdev, ideal, n_obs, sim_type) must be either the length of the longest input variable or length 1.")
		}
		# mini 1-row df for each species that will be simulated.
		var_df <- data.frame(sdev, ideal, n_obs, sim_type=as.character(sim_type))
		var_list <- lapply(X=as.list(1:nrow(var_df)), FUN=function(x){ var_df[x,] })

		# function to apply for each simulated species
		# returns a logical vector describing whether species was observed for each
		# value of env.
		env_spec_sim_1species <- function(params){
			# first, get probability of observing species with env[i]
			# given that the species is observed exactly once.
			probs_1obs <- dnorm(x=env, mean=params$ideal, sd=params$sdev)
			# handle different ways of extrapolating 
			if(any(startsWith(as.character(params$sim_type), c("F", "f")))){
				# fixed n_obs, just use sample()
				# make vector of just falses for each ENV, trues added later
				output <- rep(FALSE, length(env))
				s_inds <- sample(1:length(env), size=params$n_obs, prob=probs_1obs)
				output[s_inds] <- TRUE
				return(output)
			}else if(any(startsWith(as.character(params$sim_type), c("P", "p")))){
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

		out_matrix <- simplify2array(mclapply(X=var_list, FUN=env_spec_sim_1species, mc.cores=n_cores))
		return(list(
			matrix=out_matrix, 
			params=var_df
		))
	}




