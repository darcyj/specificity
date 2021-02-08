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
#' @aliases plot_specificities
#'
#' @examples
#' library(specificity)
#' attach(endophyte)
#' # only analyze species with occupancy >= 20
#' m <- occ_threshold(prop_abund(otutable), 20)
#' # create list to hold phy_or_env_spec outputs
#' specs_list <- list()
#' specs_list$NDVI <- phy_or_env_spec(m, env=metadata$NDVI, 
#'   n_cores=10, n_sim=50, p_method="gamma_fit")
#' specs_list$Evapotranspiration <- phy_or_env_spec(m,
#'   env=metadata$Evapotranspiration, n_cores=10, 
#'   n_sim=100, p_method="gamma_fit")
#' specs_list$Rainfall <- phy_or_env_spec(m, env=metadata$Rainfall,
#'   n_cores=10, n_sim=50, p_method="gamma_fit")
#' plot_specs_stacks(specs_list)
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
#' library(specificity)
#' attach(endophyte)
#' # only analyze species with occupancy >= 20
#' m <- occ_threshold(prop_abund(otutable), 20)
#' # create list to hold phy_or_env_spec outputs
#' specs_list <- list()
#' specs_list$NDVI <- phy_or_env_spec(m, env=metadata$NDVI, 
#'   n_cores=10, n_sim=50, p_method="gamma_fit")
#' specs_list$Evapotranspiration <- phy_or_env_spec(m,
#'   env=metadata$Evapotranspiration, n_cores=10, 
#'   n_sim=100, p_method="gamma_fit")
#' specs_list$Rainfall <- phy_or_env_spec(m, env=metadata$Rainfall,
#'   n_cores=10, n_sim=50, p_method="gamma_fit")
#' # default black
#' plot_specs_violin(specs_list)
#' # with colors
#' plot_specs_violin(specs_list, cols=c("forestgreen", "gold", "blue"))
#' # with border colors
#' plot_specs_violin(specs_list, cols=c("forestgreen", "gold", "blue"),
#'   cols_bord=c("red", "blue", "black"))
#' # with thicker borders (arg "lwd" is passed to polygon())
#' plot_specs_violin(specs_list, cols=c("forestgreen", "gold", "blue"),
#'   cols_bord="black", lwd=3)
#' # with NO borders
#' plot_specs_violin(specs_list, cols=c("forestgreen", "gold", "blue"),
#'   cols_bord=NA)
#' @export
plot_specs_violin <- function(specs_list, cols="black", cols_bord="white", 
	alpha=0.05, label_cex=1, nsig_trans=0.30, minval=-1, maxval=1, ylab="Spec", ...){
	n <- length(specs_list)
	res <- 500 # should be a user parameter; number of points to use for density fits

	# check specs_list is a list
	if(!is.list(specs_list)){
		stop("specs_list needs to be a list")
	}
	# check for names in list
	if(any(is.null(names(specs_list))) | any("" %in% names(specs_list))){
		stop("All items in specs_list must be named.")
	}
	# check for cols "Spec" and "Pval" in each item in list
	for(i in 1:length(specs_list)){
		if(!is.data.frame(specs_list[[i]])){
			stop(paste0("Item '", names(specs_list)[i], "' is not a data.frame."))
		}else if(!all(c("Spec", "Pval") %in% colnames(specs_list[[i]]))){
			stop(paste0("Item '", names(specs_list)[i], "' is missing Spec and/or Pval columns."))
		}
	}
	# check for species with NA or NaN P-values or specs
	# treat P=NA as P=1, remove Spec=NA
	for(i in 1:length(specs_list)){
		npnas <- sum(is.na(specs_list[[i]]$Pval))
		if(npnas > 0){
			warning(paste0("Item '", names(specs_list)[i], "' contains ", npnas, " NA P-values, they will be treated as P=1."))
			specs_list[[i]]$Pval[is.na(specs_list[[i]]$Pval)] <- 1
		}
		nspecnas <- sum(is.na(specs_list[[i]]$Spec))
		if(nspecnas > 0){
			warning(paste0("Item '", names(specs_list)[i], "' contains ", nspecnas, " NA Spec values, they will be ignored."))
			specs_list[[i]] <- specs_list[[i]][!is.na(specs_list[[i]]$Spec), ]
		}

	}

	# function to make dens and sigprop for each df in specs_list
	makedens <- function(x, xname){
		# check if enough points for density fit
		if(nrow(x) < 10){stop(paste0("Insufficient observations (", 
			nrow(x), "<", "10) in item '", xname, "'."))}
		fit <- density(x$Spec, bw="SJ", n=res, from=minval, to=maxval)
		dens_all <- fit$y
		yvals <- fit$x
		n_sig <- sum(x$Pval <= alpha)
		n_nsig<- sum(x$Pval > alpha)
		if(n_sig == 0){
			sigprop <- rep(0, length(dens_all))
		}else if(n_nsig==0){
			sigprop <- rep(1, length(dens_all))
		}else if(n_sig>10 & n_nsig>10){
			# scale curves to n_sig and n_nsig
			dens_sig <- density(x$Spec[x$Pval<=alpha], bw="SJ", n=res, from=minval, to=maxval)$y
			dens_sig <- (dens_sig / sum(dens_sig)) * n_sig
			dens_nsig <- density(x$Spec[x$Pval>alpha], bw="SJ", n=res, from=minval, to=maxval)$y
			dens_nsig <- (dens_nsig / sum(dens_nsig)) * n_nsig
			# add a tiny bit to avoid dividing by zero
			tinybit <- max(c(dens_sig, dens_nsig)) / 1000
			dens_sig <- dens_sig + tinybit
			dens_nsig <- dens_nsig + tinybit
			# proportion significant density
			sigprop <- dens_sig / (dens_sig + dens_nsig)
		# the following 2 conditions exist for data sets with low numbers of either sig 
		# or nsig species, which couldn't be fit to a density curve. 
		}else if(n_sig > n_nsig){
			# scale curves to n_sig and n_nsig
			dens_sig <- density(x$Spec[x$Pval<=alpha], bw="SJ", n=res, from=minval, to=maxval)$y
			dens_sig <- (dens_sig / sum(dens_sig)) * n_sig
			dens_all <- (dens_all / sum(dens_all)) * (n_sig + n_nsig)
			# add a tiny bit to avoid dividing by zero
			tinybit <- max(dens_all) / 1000
			dens_sig <- dens_sig + tinybit
			# proportion significant density
			sigprop <- dens_nsig / (dens_all + tinybit)
			# get rid of negative values
			sigprop[sigprop < 0] <- 0
		}else if(n_sig < n_nsig){
			# scale curves to n_sig and n_nsig
			dens_nsig <- density(x$Spec[x$Pval>alpha], bw="SJ", n=res, from=minval, to=maxval)$y
			dens_nsig <- (dens_nsig / sum(dens_nsig)) * n_nsig
			dens_all <- (dens_all / sum(dens_all)) * (n_sig + n_nsig)
			# add a tiny bit to avoid dividing by zero
			tinybit <- max(dens_all) / 1000
			# proportion significant density
			sigprop <- (dens_all - dens_nsig + tinybit) / (dens_all + tinybit)
			# get rid of negative values
			sigprop[sigprop < 0] <- 0
		}else{
			stop("unknown error")
		}
		return(list(dens=dens_all, sigprop=sigprop, yvals=yvals, n=nrow(x)))
	}
	# use it to make stuff
	densities <- mapply(makedens, specs_list, names(specs_list), SIMPLIFY=FALSE)

	# handle colors
	if(length(cols) == 1){
		cols <- rep(cols, n)
	}else if(length(cols) != n){
		stop("cols must be length 1, or same length as specs_list")
	}
	if(length(cols_bord) == 1){
		cols_bord <- rep(cols_bord, n)
	}else if(length(cols_bord) != n){
		stop("cols_bord must be length 1, or same length as specs_list")
	}

	# get ready for plotting
	# re-scale curves so they all add up to 1
	for(i in 1:length(densities)){
		densities[[i]]$dens <- densities[[i]]$dens / sum(densities[[i]]$dens)
	}
	# calculate max observed density
	max_dens <- max(sapply(X=densities, FUN=function(x){max(x$dens)}))
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
		x_all <- (densities[[i]]$dens / max_dens) / 2
		x_sig <- x_all * densities[[i]]$sigprop

		x_mid <- i - 0.5
		# "all" polygon
		polygon(
			x=c(x_mid-x_all + 2*x_sig, rev(x_mid + x_all )), 
			y=c(densities[[i]]$yvals, rev(densities[[i]]$yvals)),
			col=col2trans(cols[i], nsig_trans),
			border=cols_bord[i], ...
		)

		polygon(
			x=c(x_mid-x_all, rev(x_mid - x_all + 2*x_sig)), 
			y=c(densities[[i]]$yvals, rev(densities[[i]]$yvals)),
			col=cols[i],
			border=cols_bord[i], ...
		)
	}
	# labels for variables
	var_names <- names(specs_list)
	mtext(var_names, side=1, at=((1:length(var_names))-0.5), cex=label_cex )
}

