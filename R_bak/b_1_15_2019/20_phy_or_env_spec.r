#' phy_or_env_spec
#'
#' Calculates standardized effect size (SES) of host phylogenetic specificity (per 
#' Poulin et al. 2011) or standardized effect size of environmental specificity
#' for a species i. Also returns a P-value calculated via simulation, which has a
#' minimum possible value of 1/n_sim. Note that the calculation of SPS/SES is '
#' meaningless unless all observations where species i MAY HAVE been found are 
#' included, including observations where species i was ABSENT. 
#' 
#' For phylogenetic specificity, lower SES values indicate host phylogenetic 
#' underdispersion (i.e. strong host specificity) of species i, values close to 
#' zero indicate random host association, and values above zero indicate host 
#' phylogenetic overdispersion (more cosmopolitanism) than expected by chance.
#'
#' For environmental specificity, lower SES values indicate that when species i
#' is present, environmental variability is lower than expected by chance
#' (environmental specificity). SES values close to zero indicate random
#' environmental association, and SES values above zero indicate more environmental
#' variability than expected by chance when species i is present(cosmopolitanism).
#' 
#' Because SES can be either negative (specificity), 0 (null), or positive
#' (cosmopolitanism), the statistical test can be 1- or 2-tailed. Default is to
#' only test for specificity.
#'
#' @author John L. Darcy
#' @seealso phy_or_env_spec_multiple() paralellizes this function across multiple 
#' species that were found across the same hosts or environmental variable.
#' @references 
#' Poulin et al. (2011) Host specificity in phylogenetic and geographic
#'   space. Trends Parasitol 8:355-361. doi: 10.1016/j.pt.2011.05.003
#'
#' @param abunds numeric vector. Abundances for species i, corresponding to either
#' hosts or env.
#' @param hosts character vector. Host identities corresponding to abunds. Only
#' required if calculating SES for phylogenetic specificity.
#' @param hosts_phylo phylo object. Tree containing all unique hosts as tips. Only
#' required if calculating SES for phylogenetic specificity.
#' @param env numeric vector or dist or square matrix. Environmental variable. In 
#'   the case of a vector, must be same length as abunds. For a dist or square
#'   matrix, the matrix's dimensions must both be equal to the length of abunds. 
#'   A data.frame can be used, too, as long as as.matrix(env) yields the appropriate
#'   square matrix input. This argument is only required if calculating SES for 
#'   environmental specificity. 
#' @param unweighted logical. If TRUE, will use PD instead of WPD for phylogenetic
#'   specificity (DEFAULT: FALSE).
#' @param n_sim integer. number of simulations to do under the null hypothesis that
#'   host or environmental association is random (DEFAULT: 1000)
#' @param seed integer. Seed to use so that this is repeatable (DEFAULT: 12345)
#' @param tails integer. 1 = 1-tailed, test for specificity only. 2 = 2-tailed.
#'   3 = 1-tailed, test for cosmopolitanism only. 4 = no test, P=1.0 (DEFAULT: 1).
#' @param occ_stab logical. If TRUE, output statistic is stabilized for occupancy.
#'   if FALSE, not stabilized as in Poulin et al. 2011 (DEFAULT=TRUE).
#' @param occ_stab_frac float. Fractional occupancy between 0 and 1 for occupancy-
#'   stable calculation of denominator. Values between 0.10 and 0.40 work best
#'   (DEFAULT: 0.20).
#' @param force_denom float. If nonzero and occ_stab==TRUE, occupancy stabilization
#'   is done by forcing SES denominator to be a given value. Use 1 for non-
#'   standardized effect size (units of env or pd) or use 0 for denominator to be
#'   simulated using occ_stab_frac. Useful to make the code run faster by pre-
#'   computing the denominator if it's constant (DEFAULT: 0).
#' @param denom_only logical. If TRUE, only SES denominator will be returned. This
#'   option should only be used to calculate the denominator for the force_denom
#'   option if this is run multiple times. Changes return to be a single denominator
#'   value (DEFAULT: FALSE).
#'
#' @return vector of P-value and SPS.
#'
#' @examples
#' # none yet written.
#'
#' @export
	phy_or_env_spec <- function(abunds, hosts=NULL, hosts_phylo=NULL, env=NULL, 
		unweighted=FALSE, n_sim=1000, seed=12345, tails=1, occ_stab=TRUE, 
		occ_stab_frac=0.20, force_denom=0, denom_only=FALSE){
		set.seed(seed)
		# make sure there are enough abunds to actually calculate variance or pd:
		if(sum(abunds > 0) < 2 && denom_only==FALSE){
			pval <- NA; ses <- NA; denom <- NA
		}else{
			# decide whether to calculate phy or env specificity:
			if((is.null(hosts) || is.null(hosts_phylo)) && !is.null(env)){
				# either hosts or tree is null, and env is provided. Doing env spec.
				# but is it a matrix or a vector?
				if(class(env) %in% c("matrix", "dist", "data.frame")){
					spec_type <- "env_dist"
					env <- as.dist(env)
				}else{
					spec_type <- "env_vec"
				}
			}else if(is.null(env) && (!is.null(hosts) && !is.null(hosts_phylo))){
				# env is null and hosts was given. doing phy spec.
				spec_type <- "phy"
			}else if(!is.null(hosts) && is.null(hosts_phylo)){
				stop("ERROR: hosts provided but hosts_phylo is missing.")
			}else if(is.null(env) && is.null(hosts)){
				stop("ERROR: hosts or env is missing.")
			}else{
				stop("ERROR: mystery input error.")
			}
			# check occ_stab_frac
			if(occ_stab_frac > 1 || occ_stab_frac < 0){
				stop("ERROR: occ_stab_frac not within [0,1].")
			}

			# calculate empirical specificity and simulate n_sim specificities under null hypothesis
			if(spec_type == "phy"){
				spec <- host_wpd(abunds, hosts, hosts_phylo, null=FALSE, unweighted)
				specs_sim <- replicate(n_sim, host_wpd(abunds, hosts, hosts_phylo, null=TRUE, unweighted))
				if(occ_stab == TRUE && force_denom <= 0){
					# occ_stab is true, but no forced denominator. use occ_stab_frac.
					occ_not <- round(length(abunds) * occ_stab_frac)
					denom <- sd(replicate(n_sim, host_wpd(abunds, hosts, hosts_phylo, null=TRUE, null_occ=occ_not, unweighted)))
				}else if(occ_stab == TRUE && force_denom > 0){
					denom <- force_denom
				}else{
					denom <- sd(specs_sim)
				}
			}else if(spec_type == "env_vec"){
				spec <- var(env[abunds > 0])
				specs_sim <- replicate(n_sim, var(env[sample(abunds > 0)]))
				if(occ_stab == TRUE && force_denom <= 0){
					# occ_stab is true, but no forced denominator. use occ_stab_frac.
					occ_not <- round(length(abunds) * occ_stab_frac)
					denom <- sd(replicate(n_sim, var(sample(env, size=occ_not))))
				}else if(occ_stab == TRUE && force_denom > 0){
					denom <- force_denom
				}else{
					denom <- sd(specs_sim)
				}
			}else if(spec_type == "env_dist"){
				subdist <- function(d, v){
					d[unlist(sapply(X=1:(length(v)-1), FUN=function(x){v[1:x]}))]
				}
				emp_pres <- (abunds > 0)
				spec <- var(subdist(env, emp_pres))
				specs_sim <- replicate(n_sim, subdist(env, sample(emp_pres)])))
				if(occ_stab == TRUE && force_denom <= 0){
					# occ_stab is true, but no forced denominator. use occ_stab_frac.
					occ_not <- round(length(abunds) * occ_stab_frac)
					occ_pres <- rep(F, length(abunds))
					occ_pres[sample(1:length(occ_pres), size=occ_not)] <- TRUE
					denom <- sd(replicate(n_sim, subdist(env, sample(occ_pres)]))))
				}
			}else{ 
				stop("ERROR: mystery error")
			}


			# calculate SES and P-val only if denom_only is FALSE
			if(denom_only == FALSE){
				# calculate SES
				ses <- (spec - mean(specs_sim)) / denom

				# calculate p-value
				num_sims_below <- sum(specs_sim < spec)
				num_sims_above <- sum(specs_sim > spec)
				# if 0, bump it up to 1 because we don't actually KNOW if p == 0
				if(num_sims_below < 1){num_sims_below <- 1}
				if(num_sims_above < 1){num_sims_above <- 1}
				# p-value for each tails case:
				if(tails == 1){
					pval <- (num_sims_below / n_sim)
				}else if(tails == 3){
					pval <- (num_sims_above / n_sim)
				}else if(tails == 2){
					pval <- min(c( (num_sims_below / n_sim), (num_sims_above / n_sim) ))
				}else{
					pval <- 1
				}
				return(c(Pval=pval, SES=ses))
			}else{
				return(denom)
			}

		}

	}
