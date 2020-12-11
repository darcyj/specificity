#' plot_pairwise_spec
#'
#' Plots pairwise correlations between specificity to multiple variables. Specificity
#' results are supplied to this function as a list of specificity tables, i.e. a list
#' where each object within the list is an output of phy_or_env_spec, and all were
#' created using the same abunds_mat object (see: ?phy_or_env_spec). 
#'
#' @author John L. Darcy
#' 
#' @param sl "specs list" list of outputs from phy_or_env_spec as described above.
#' @param label_cex float. Size of variable labels, which will be displayed along the
#'   plot's diagonal. Use cex units; see ?par (default: 1).
#' @param point_cex float. Size of points in the plot's lower triangle. Useful to reduce
#'   this if you are plotting lots of species. Use cex units; see ?par (default: 1).
#' @param cor_cex float. Size of text for correlations displayed in plot's upper
#'   triangle. Use cex units; see ?par (default: 1).
#' @param cor_red_lim float. Correlation coefficients will be shown in red if they are
#'   equal to or more extreme than this value (default: 0.70). 
#' @param method string. Preferred correlation method. see ?cor for options (default: 
#'   "pearson").
#' @return Returns nothing. Plots correlations in a square matrix of subplots, where
#'   variable names are shown in the diagonal, pairwise specificities are plotted in
#'   the lower triangle, and correlation coefficients are displayed in the upper
#'   triangle. For plots in the lower triangle, each point represents a species.
#'
#' @examples
#' # # example commented out since they are computationally intense
#' # # this is so they don't cause testing the package to take forever.
#' # attach(endophyte)
#' # otutable_over10 <- occ_threshold(otutable, threshold = 10)
#' # specs_list <- list()
#' # specs_list$NDVI <- phy_or_env_spec(otutable_over10, env=metadata$NDVI, 
#' #   n_cores=10, n_sim=100, p_method="gamma_fit")
#' # specs_list$Evapotranspiration <- phy_or_env_spec(otutable_over10,
#' #   env=metadata$Evapotranspiration, n_cores=10, n_sim=100, p_method="gamma_fit")
#' # specs_list$Rainfall <- phy_or_env_spec(otutable_over10, env=metadata$Rainfall,
#' #   n_cores=10, n_sim=100, p_method="gamma_fit")
#' # plot_pairwise_spec(specs_list)
#'
#' @export
plot_pairwise_spec <- function(sl, label_cex=1, point_cex=1, cor_cex=2, cor_red_lim=0.70, method="pearson"){
    # convert sl to a data.frame for easier plotting
    df <- simplify2array(lapply(X=sl, FUN=function(x){x$Spec}))

    # store old par so it can be resetted
    old.par <- par(no.readonly = TRUE)
    # set new par
    n <- ncol(df)
    par(mfrow = c(n,n), oma = c(2,2,2,2), mar = c(0,0,0,0) )

    # make a matrix to figure out which type of plot to do at position i,j
    # lower tri = scatterplots, diag=names, upper tri = correlation coefficients
    typemat <- matrix("D", nrow=n, ncol=n)
    typemat[lower.tri(typemat)] <- "L"
    typemat[upper.tri(typemat)] <- "U"
    for(i in 1:n){for(j in 1:n){
        if(typemat[i,j] == "L"){
            # lower tri - do scaterplot
            plot(x=df[,j], y=df[,i], axes = FALSE, xlab="", ylab="", pch=20, cex=point_cex)
            abline(v=0)
            abline(h=0)
            box()
        }else if(typemat[i,j] == "D"){
            # diag - write variable name
            plot(1, type="n", xlim=c(-1, 1), ylim=c(-1, 1), axes = FALSE, xlab="", ylab="", pch=20)
            text(x=0, y=0, labels=colnames(df)[i], cex=label_cex, srt=-45)
            box()
        }else if(typemat[i,j] == "U"){
            # upper tri - nicely display correlation coefficient (r)
            cor_ij <- cor(df[,j], df[,i], use="complete.obs", method=method) 
            if(cor_ij >= cor_red_lim || cor_ij <= (-1 * cor_red_lim)){
                col_ij <- "red"
            }else{
                col_ij <- "black"
            }
            cor_ij <- sprintf("%.2f", round(cor_ij,2))

            plot(1, type="n", xlim=c(-1, 1), ylim=c(-1, 1), axes = FALSE, xlab="", ylab="", pch=20)
            text(x=0, y=0, labels=cor_ij, cex=cor_cex, col=col_ij)
            box()
        }
    }}

    # reset par back to normal
    par(old.par)

}
