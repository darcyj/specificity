
#' onto2nwk
#'
#' Converts an ontology (higherarchical categories) into a nwk phylogeny.
#'
#' @author John L. Darcy
#'
#' @param df a data.frame object where columns represent ontology levels, which
#'   are assumed to be nested hierarchically. this function does not check for 
#'   proper hierarchical nestedness - it is the user's job to check that each 
#'   node and tip name is monophyletic. Lower levels (e.g. tips) should be the 
#'   rightmost column of df, and higher levels (e.g. roots) should be leftmost
#'   column, with intermediate columns ordered between.
#'
#' @return A newick (nwk) format string.
#'
#' @examples
#' df <- data.frame(
#'	l1 = c( "a", "a", "a", "a", "a", "a", "a", "b", "b", "b", "b", "b", "b", "c", "d"),
#'	l2 = c( "e", "e", "e", "e", "f", "f", "g", "h", "h", "h", "i", "j", "j", "k", "l"),
#'	l3 = c( "m", "n", "o", "o", "p", "p", "q", "r", "r", "s", "t", "u", "v", "w", "x")
#' )
#' nwk_str <- onto2nwk(df)
#' a <- ape::read.tree(text=nwk_str)
#' plot(a, show.node.label=TRUE)
#'
#' @export
onto2nwk <- function(df){
	# get rid of identical rows
	df <- df[!duplicated(df),]
	# add root if not rooted
	if(all(df[,1] == df[1,1])){
		# do nothing, it's rooted (maybe do something here later)
	}else{
		# root
		df <- data.frame(root=rep("root", nrow(df)), df)
	}
	# nwkize function
	nwkize <- function(nodes){paste0("(", paste(paste0(nodes, ":1"), collapse=","), ")")}
	for(col in ncol(df):2){
		# agg col by col-1
		ag <- aggregate(df[,col], by=list(df[,col-1]), nwkize)$x
		# agg df by col-1
		rows2keep <- aggregate(1:nrow(df), by=list(df[,col-1]), function(x){x[1]})$x
		df <- df[rows2keep,]
		df[,col-1] <- paste0(ag, df[,col-1])

	}
	return(paste0(df[1,1], ";"))
}