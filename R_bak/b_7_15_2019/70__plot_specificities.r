#' plot_specificities
#'
#' Visualizes results from phy_or_env_spec
#'
#' @author John L. Darcy
#'
#' @param specs_list list of data.frames. Each data.frame must be an output from
#'   phy_or_env_spec; must have columns "SES" and "Pval".
#' @param n_bins integer. Number of bins for stacked violins (DEFAULT: 20).
#' @param col_sig string. Color name or hex code for species where Pval <= alpha
#'   (DEFAULT = "black").
#' @param col_nsig string. Color name or hex code for species where Pval > alpha
#'   (DEFAULT = "gray").
#' @param col_bord string. Color name or hex code for border color. Use NA for no
#'   border (DEFAULT = NA). 
#' @param alpha float. alpha value for determining statistical significance; see
#'   col_sig and col_nsig above (DEFAULT = 0.05).
#' @return returns nothing (a plot is made).
#'
#' @examples
#'   none yet written.
#'
#' @export
	plot_specificities <- function(specs_list, n_bins=20, col_sig="black", col_nsig="gray", 
		col_bord=NA, alpha=0.05, label_cex=0.60){

		# find min and max specificities
		min_spec <- min( sapply(X=specs_list, FUN=function(x){min(x$SES)}) )
		max_spec <- max( sapply(X=specs_list, FUN=function(x){max(x$SES)}) )

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

			spec_sig <- spec$SES[spec$Pval <= alpha]
			spec_nsig<- spec$SES[spec$Pval > alpha]
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
		plot(0, type="n", xlab="", ylab="SES", xaxt='n', ylim=c(min_spec, max_spec), xlim=c(0.5, n_vars+0.5))

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
			rect(xmin, ymin, split, ymax, col=col_sig, border=col_bord, lwd=0.5)
			# rectangle for nsig, from split to xmax
			rect(split, ymin, xmax, ymax, col=col_nsig, border=col_bord, lwd=0.5)
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
