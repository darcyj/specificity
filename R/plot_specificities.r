#' plot_specs_stacks
#'
#' Visualizes results from phy_or_env_spec as stacked histograms. Aliased to 
#' plot_specificities() for backward compatibility. 
#'
#' @author John L. Darcy
#'
#' @param specs_list list of data.frames. Each data.frame must be an output from
#'   phy_or_env_spec; must have columns "Spec" and "Pval".
#' @param n_bins integer. Number of bins for stacked violins (default: 20).
#' @param col_sig string. Color name or hex code for species where Pval <= alpha
#'   (default: "black").
#' @param col_nsig string. Color name or hex code for species where Pval > alpha
#'   (default: "gray").
#' @param col_bord string. Color name or hex code for border color. Use NA for no
#'   border (default: NA). 
#' @param alpha float. alpha value for determining statistical significance; see
#'   col_sig and col_nsig above (default: 0.05).
#' @param label_cex float. Used to change size of x-axis labels (default: 1).
#' @return returns nothing (a plot is made).
#'
#' @examples
#' data(endophyte); attach(endophyte)
#' otutable <- occ_threshold(prop_abund(otutable), 10)
#' elev_spec <- phy_or_env_spec(otutable, env=metadata$Elevation, n_sim=100, p_method="gamma_fit")
#' ndvi_spec <- phy_or_env_spec(otutable, env=metadata$NDVI, n_sim=100, p_method="gamma_fit")
#' plot_specs_stacks(list(Elevation=elev_spec, NDVI=ndvi_spec))
#'
#' @export
plot_specs_stacks <- plot_specificities <- function(specs_list, n_bins=20, 
	col_sig="black", col_nsig="gray", col_bord=NA, alpha=0.05, label_cex=1){

	# find min and max specificities
	min_spec <- min( sapply(X=specs_list, FUN=function(x){min(x$Spec)}) )
	max_spec <- max( sapply(X=specs_list, FUN=function(x){max(x$Spec)}) )

	# make binning table
	bin_table <- data.frame(
		bin = 1:n_bins,
		mins = seq(from=min_spec, to=max_spec, length.out=n_bins+1)[1:n_bins],
		maxs = seq(from=min_spec, to=max_spec, length.out=n_bins+1)[2:(n_bins+1)]
	)
	# add a little bit to top of final bin so we can use <
	bin_table$maxs[n_bins] <- bin_table$maxs[n_bins] + ((max_spec - min_spec) / 100)

	# function to do binning
	bin_spec_data <- function(spec, bin_table){

		spec_sig <- spec$Spec[spec$Pval <= alpha]
		spec_nsig<- spec$Spec[spec$Pval > alpha]
		# count how many are in each bin
		count_sig <- count_nsig <- rep(0, n_bins)
		for(i in 1:n_bins){
			count_sig[i] <- sum(spec_sig >= bin_table$mins[i] & spec_sig < bin_table$maxs[i])
			count_nsig[i] <- sum(spec_nsig >= bin_table$mins[i] & spec_nsig < bin_table$maxs[i])
		}
		# convert to proportions and return
		return(data.frame(
			bin=bin_table$bin,
			sig=count_sig / nrow(spec),
			nsig=count_nsig / nrow(spec)
		))

	}

	# make a big df of binned proportions
	bin_props_list <- lapply(X=specs_list, FUN=bin_spec_data, bin_table=bin_table)
	# get variable names in there and transform to df
	for(i in 1:length(bin_props_list)){
		bin_props_list[[i]] <- data.frame(var= rep(names(bin_props_list)[i], nrow(bin_props_list[[i]])), bin_props_list[[i]])
	}
	plotdf <- do.call("rbind", bin_props_list)
	rownames(plotdf) <- 1:nrow(plotdf)

	# scale plotdf so max bin width is 1
	max_bin_width <- max(plotdf$sig + plotdf$nsig)
	plotdf$sig <- plotdf$sig * (1/max_bin_width)
	plotdf$nsig<- plotdf$nsig* (1/max_bin_width)

	# make plot area
	n_vars <- length(specs_list)
	plot(0, type="n", xlab="", ylab="Spec", xaxt='n', 
		ylim=c(min_spec, max_spec), xlim=c(0.5, n_vars+0.5))

	# function to plot boxes in a row
	# x is a row from plotdf, m is the middle of the stack (x-axis)
	plotboxes <- function(x, m){
		width <- sum(x$sig, x$nsig)
		ymin <- bin_table$min[x$bin]
		ymax <- bin_table$max[x$bin]
		xmin <- m - (width/2)
		xmax <- m + (width/2)
		split <- xmin + x$sig
		# rectangle for sig, from xmin to split
		if(x$sig > 0.0001){ # prevents width 0 rect, which is drawn as a line in some pdf viewers!
			rect(xmin, ymin, split, ymax, col=col_sig, border=col_bord, lwd=0.5)
		}
		# rectangle for nsig, from split to xmax
		if(x$nsig > 0.0001){ # prevents width 0 rect, which is drawn as a line in some pdf viewers!
			rect(split, ymin, xmax, ymax, col=col_nsig, border=col_bord, lwd=0.5)
		}
	}

	# for each variable, plotboxes
	var_names <- unique(plotdf$var)
	for(i in 1:length(var_names)){
		var_i <- var_names[i]
		subdf <- plotdf[plotdf$var == var_i,]
		sublist <- split(subdf, seq(nrow(subdf)))
		lapply(X=sublist, FUN=plotboxes, m=i)
	}

	# labels for variables
	axis(side=1, at=1:length(var_names), labels=var_names, cex.axis=label_cex)
}



#' plot_specs_violin
#'
#' Visualizes results from phy_or_env_spec as violins. Violin area is proportionally 
#' divided such that lighter colors represent density of non-significant features, 
#' and darker colors represent statistically significant features.
#'
#' @author John L. Darcy
#'
#' @param specs_list list of data.frames. Each data.frame must be an output from
#'   phy_or_env_spec; must have columns "Spec" and "Pval".
#' @param cols character vector of color names or hex codes. If only one value is
#'   given, all violins will be that color. Otherwise, one value may be given per
#'   item in specs_list, corresponding to its order (default: "black").
#' @param cols_bord character vector of color names or hex codes. Color name or hex
#'   code for borders drawn around and within violins. Length 1 or length n, just
#'   like cols. For no borders, use cols_bord=NA (default: "white").
#' @param alpha float. alpha value for determining statistical significance (default: 0.05).
#' @param label_cex float. Used to change size of x-axis labels (default: 1).
#' @param nsig_trans float between 0 and 1 (inclusive). Determines how transparent
#'   violin area will be for nonsignificant features, with 0 meaning totally
#'   transparent and 1 meaning totally opaque (default: 0.4).
#' @param minval minimum possible value for specificity statistic (default: -1).
#' @param maxval maximum possible value for specificity statistic (default: 1).
#' @param ylab y-axis label for plot (default:"Spec").
#' @param ... additional arguments to be passed to polygon().
#' @return returns nothing (a plot is made).
#'
#' @examples
#' data(endophyte); attach(endophyte)
#' otutable <- occ_threshold(prop_abund(otutable), 10)
#' elev_spec <- phy_or_env_spec(otutable, env=metadata$Elevation, 
#'   n_sim=100, p_method="gamma_fit")
#' ndvi_spec <- phy_or_env_spec(otutable, env=metadata$NDVI, 
#'   n_sim=100, p_method="gamma_fit")
#' # default black
#' plot_specs_violin(list(Elevation=elev_spec, NDVI=ndvi_spec))
#' # with colors
#' plot_specs_violin(list(Elevation=elev_spec, NDVI=ndvi_spec), 
#'   cols=c("orange", "forestgreen"))
#' # with border colors
#' plot_specs_violin(list(Elevation=elev_spec, NDVI=ndvi_spec), 
#'   cols=c("orange", "forestgreen"), 
#'   cols_bord=c("blue", "red"))
#' # with thicker borders (arg "lwd" is passed to polygon())
#' plot_specs_violin(list(Elevation=elev_spec, NDVI=ndvi_spec), 
#'   cols=c("orange", "forestgreen"), cols_bord="black", lwd=3)
#' @export
plot_specs_violin <- function(specs_list, cols="black", cols_bord="white", 
	alpha=0.05, label_cex=1, nsig_trans=0.30, minval=-1, maxval=1, ylab="Spec", ...){
	n <- length(specs_list)
	# for each spec result in specs_list, make a density fit using
	# the exact same "x" (soon to be y) scale.
	densities_all <- lapply(X=specs_list, FUN=function(x){
		density(x$Spec, bw="SJ", n=500, from=minval, to=maxval)})
	densities_sig <- lapply(X=specs_list, FUN=function(x){
		density(x$Spec[x$Pval < alpha], bw="SJ", n=500, from=minval, to=maxval)})
	# simplify objects
	densities_all_mat <- sapply(X=densities_all, FUN=function(x){x$y}, simplify="array")
	densities_sig_mat <- sapply(X=densities_sig, FUN=function(x){x$y}, simplify="array")
	yvals <- densities_all[[1]]$x
	# convert densities mats to proportion, just in case density has some variation
	# in total densities due to missing data etc
	densities_all_mat <- apply(X=densities_all_mat, MAR=2, FUN=function(x){x/sum(x)})
	densities_sig_mat <- apply(X=densities_sig_mat, MAR=2, FUN=function(x){x/sum(x)})
	# proportionalize densities_sig
	sig_props <- sapply(X=specs_list, FUN=function(x){sum(x$Pval < alpha) / nrow(x)})
	for(i in 1:ncol(densities_sig_mat)){
		densities_sig_mat[,i] <- densities_sig_mat[,i] * sig_props[i]
	}
	# reconcile larger sig densities than all
	sig_higher <- densities_sig_mat > densities_all_mat
	densities_sig_mat[sig_higher] <- densities_all_mat[sig_higher]
	# handle colors
	if(length(cols) == 1){
		cols <- rep(cols, n)
	}else if(length(cols) == n){
		# do nothing...?
	}else{
		stop("cols must be length 1, or same length as specs_list")
	}
	if(length(cols_bord) == 1){
		cols_bord <- rep(cols_bord, n)
	}else if(length(cols_bord) == n){
		# do nothing...?
	}else{
		stop("cols_bord must be length 1, or same length as specs_list")
	}
	max_dens <- max(densities_all_mat)
	# calculate plot yaxis limits
	min_spec <- min(sapply(X=specs_list, FUN=function(x){min(x$Spec)}))
	max_spec <- min(sapply(X=specs_list, FUN=function(x){max(x$Spec)}))
	# function to make color transparent
	col2trans <- function(col, trans){
		rgb( col2rgb(col)[1,1]/255, col2rgb(col)[2,1]/255, col2rgb(col)[3,1]/255, trans)
	}
	# each subplot is 1 unit wide, centered on n-0.5. 
	plot(0, type="n", xlim=c(0, n), ylim=c(min_spec, max_spec), xaxt="n", xlab="", ylab=ylab,
		bty="n")
	# for each variable, draw a violin
	for(i in 1:n){
		x_all <- (densities_all_mat[,i] / max_dens) / 2
		x_sig <- (densities_sig_mat[,i] / max_dens) / 2

		x_mid <- i - 0.5
		# "all" polygon

		polygon(
			x=c(x_mid-x_all + 2*x_sig, rev(x_mid + x_all )), 
			y=c(yvals, rev(yvals)),
			col=col2trans(cols[i], nsig_trans),
			border=cols_bord[i], ...
		)

		polygon(
			x=c(x_mid-x_all, rev(x_mid - x_all + 2*x_sig)), 
			y=c(yvals, rev(yvals)),
			col=cols[i],
			border=cols_bord[i], ...
		)
	}
	# labels for variables
	var_names <- names(specs_list)
	mtext(var_names, side=1, at=((1:length(var_names))-0.5), cex=label_cex )
}




