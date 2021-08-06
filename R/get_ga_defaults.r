#' get_ga_defaults
#'
#' Simply returns default parameters for the genetic algorithm in rao_genetic_max().
#' This function has no arguments.
#' 
#' @author John L. Darcy
#'
#' @return named list of parameters.
#'
#' @examples
#' # get_ga_defaults()
#'
#' @export
get_ga_defaults <- function(){
	list(
		swap_freq=c(1,1,2,3),
		term_cycles=10,
		maxiters=400,
		popsize_swap=150,
		popsize_perm=150,
		keep=5,
		prc=0.0001
	)
}



