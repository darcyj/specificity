#' get_ga_defaults
#'
#' Simply returns default parameters for the genetic algorithm in rao_genetic_max().
#' This function has no parameters.
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
		term_cycles=10,
		maxiters=400,
		popsize=300,
		keep=5,
		prc=0.01
	)
}