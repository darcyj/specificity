#' check_pes_inputs
#' 
#' Function used by phy_or_env_spec. 
#' checks abunds_mat, env, hosts, and hosts_phylo inputs to phy_or_env_spec to
#' make sure there are no problems. This could include missing species in trees,
#' incompatible dimensions, non-numeric inputs, etc. Returns an input type, which
#' is just a string that can be "mat", "dist", "vec", "phy", or "error".
#' 
#' @param abunds_mat (required, see phy_or_env_spec)
#' @param env (required, can be NULL, see phy_or_env_spec)
#' @param hosts (required, can be NULL, see phy_or_env_spec)
#' @param hosts_phylo (required, can be NULL, see phy_or_env_spec)
#' @param verbose logical. Should status messages be displayed? (default: TRUE).
#' 
#' @return string. either "mat", "dist", "vec", "phy", or "error".
#' 
#' @examples
#' library(specificity)
#' attach(endophyte)
#' m <- occ_threshold(prop_abund(otutable), threshold=10)
#' check_pes_inputs(m, env=metadata$Elevation, hosts=NULL, hosts_phylo=NULL)
#' check_pes_inputs(m, env=NULL, hosts=metadata$PlantGenus, hosts_phylo=supertree)
#' aspect_dis <- circularize2dist(metadata$Aspect, 360)
#' check_pes_inputs(m, env=aspect_dis, hosts=NULL, hosts_phylo=NULL)
#' 
#' @export
check_pes_inputs <- function(abunds_mat, env, hosts, hosts_phylo, verbose=TRUE){
	msg <- function(x){if(verbose){message(x)}}
	# decide inputs type
	env_prsnt <- !is.null(env)
	hos_prsnt <- !is.null(hosts)
	phy_prsnt <- !is.null(hosts_phylo)
	if(env_prsnt && (hos_prsnt || phy_prsnt)){
		# both env AND some phylogenetic stuff present, error!
		data_type <- "error"
		msg("Error: Too many inputs. Pick either env only, or hosts+hosts_phylo.")
	}else if(env_prsnt && (is.matrix(env) || is.data.frame(env))){
		# env present and it's 2-dimensional
		data_type <- "mat"
	}else if(env_prsnt && class(env) == "dist"){
		# env present and it's a dist (LT from 2-dimensional)
		data_type <- "dist"
	}else if(env_prsnt && is.vector(env)){
		# env present and it's a vector
		data_type <- "vec"
	}else if(env_prsnt && is.factor(env)){
		data_type <- "error"
		msg("Error: env is a factor.")
	}else if(hos_prsnt && phy_prsnt){
		# hosts and hosts_phylo present, it's phy
		data_type <- "phy"
	}else if(xor(hos_prsnt, phy_prsnt)){
		# only one of hosts or hosts_phylo are present, but no env. error!
		data_type <- "error"
		msg("Error: Both hosts and hosts_phylo required for phylogenetic analysis.")
	}else if(! env_prsnt || hos_prsnt){
		# none of them are present!
		data_type <- "error"
		msg("Error: Either env or hosts is required.")
	}else{
		# mystery error!
		data_type <- "error"
		msg("Error: Mystery error.")
	}

	# check for NA or NAN values in abunds_mat
	if(anyNA(abunds_mat)){
		data_type <- "error"
		msg("Error: NA or NaN values in abunds_mat.")
	}

	# check for errors specific to 2D data:
	if(data_type %in% c("mat", "dist")){
		# if it's a dist, just convert to mat. all same checks apply.
		if(data_type == "dist"){
			env <- as.matrix(env)
		}
		# is env numeric?
		if(!is.numeric(env)){
			data_type <- "error"
			msg("Error: env not numeric.")
		# is env square?
		}else if(nrow(env) != ncol(env)){
			data_type <- "error"
			msg("Error: env is 2-dimensional, but not square.")
		# is env same dim as abunds_mat?
		}else if(nrow(env) != nrow(abunds_mat)){
			data_type <- "error"
			msg("Error: env and abunds_mat have incompatible dimensions.")
		# are there NA/NaN values in env?
		}else if(anyNA(env)){
			data_type <- "error"
			msg("Error: env contains NA or NaN values.")
		}
	}

	# check for errors specific to 1D env data:
	if(data_type == "vec"){
		# is env numeric?
		if(!is.numeric(env)){
			data_type <- "error"
			msg("Error: env not numeric.")
		# is env same dim as abunds_mat?
		}else if(length(env) != nrow(abunds_mat)){
			data_type <- "error"
			msg("Error: env and abunds_mat have incompatible dimensions.")
		# are there NA/NaN values in env?
		}else if(anyNA(env)){
			data_type <- "error"
			msg("Error: env contains NA or NaN values.")
		}
	}

	# check for errors specific to phylogenetic data:
	if(data_type == "phy"){
		# are there NA/NaN values in hosts?
		if(anyNA(hosts)){
			data_type <- "error"
			msg("Error: hosts contains NA or NaN values.")
		# are all the right tips in hosts_phylo?
		}else if(! all(hosts %in% hosts_phylo$tip.label)){
			data_type <- "error"
			msg("Some hosts are missing from hosts_phylo.")
		# is hosts the right length?
		}else if(length(hosts) != nrow(abunds_mat)){
			data_type <- "error"
			msg("hosts and abunds_mat have incompatible dimensions.")
		# do any hosts contain semicolons?
		}else if(any(grepl(";", x=hosts_phylo$tip.label))){
			data_type <- "error"
			msg("hosts cannot contain semicolons.")
		}
	}

	return(data_type)
}

