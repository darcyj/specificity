#' aggregate_specs_list
#'
#' Aggregates a list of outputs from phy_or_env_spec() into a single data.frame object.
#' Can also include feature data (e.g. species taxonomy) into that output. Output can 
#' also be byFeature, with one row per feature, and multiple columns for different 
#' variables.
#'
#' @author John L. Darcy
#'
#' @param sl specs_list. A named list of outputs from phy_or_env_spec. See examples.
#' @param byFeature bool. If true, each feature will occupy only one row, with multiple
#'   columns to represent the different variables in specs_list (default: FALSE)
#' @param fd data.frame. Optional feature data - a data.frame object with one row per 
#'   feature, including some column with feature IDs that includes feature IDs in sl as
#'   rownames (default:NULL)
#' @param fd_id integer or string. If integer, specifies the column index of fd that 
#'   contains feature ids. If character, specifies the column name (default: 1).
#' 
#' @return a data.frame object.
#'
#' @examples
#' # attach(endophyte)
#' # otutable <- occ_threshold(prop_abund(otutable), 20)
#' # specs_list <- list()
#' # # note: "index_rough" is only being used here to save time for demonstration purposes.
#' # specs_list$elevation <- phy_or_env_spec(otutable, env=metadata$Elevation, 
#' #   n_cores=20, n_sim=100, denom_type="sim_center")
#' # specs_list$rainfall <- phy_or_env_spec(otutable, env=metadata$Rainfall, 
#' #   n_cores=20, n_sim=100, denom_type="sim_center")
#' # # aggregate, long mode, like for ggplot:
#' # specs_df_long <- aggregate_specs_list(specs_list, byFeature=FALSE, fd=taxonomy, fd_id=1)
#' # # aggregate, wide mode:
#' # specs_df_wide <- aggregate_specs_list(specs_list, byFeature=TRUE, fd=taxonomy, fd_id=1)
#' # # example plot with ggplot:
#' # library(ggplot2)
#' # ggplot(specs_df_long, aes(x=Variable, y=Spec)) + geom_violin() + geom_jitter(width=0.3)
#'
#' @export
  aggregate_specs_list <- function(sl, byFeature=FALSE, fd=NULL, fd_id=1){
    # validate sl as a specs list
    if(!is.list(sl)){
      stop("sl is not a list")
    }else if( !all(sapply(X=sl, FUN=class) %in% c("data.frame", "matrix"))){
      stop("not every object in sl is data.frame or matrix")
    }else if( !all(sapply(X=sl, FUN=function(x){all(c("Pval", "Spec") %in% colnames(x))})) ){
      stop("not every object in sl has Spec and Pval columns")
    }else if(any(is.null(names(sl))) | any(names(sl)=="")){
      stop("names missing from sl")
    }else if(length(names(sl)) != length(unique(names(sl)))){
      stop("sl names not unique")
    }

    # convert to "long" format
    # add "Variable"
    sl <- mapply(FUN=function(df, nm){df$Variable<-nm; df}, df=sl, nm=names(sl), SIMPLIFY=FALSE)
    # add "FeatureID"
    sl <- lapply(X=sl, FUN=function(df){df$FeatureID <- rownames(df); df})
    # remove other columns (may be present if diagnostic mode was used in phy_or_env_spec()) and sort
    sl <- lapply(X=sl, FUN=function(df){df[, colnames(df) %in% c("Pval", "Spec", "Variable", "FeatureID")]})
    sl <- lapply(X=sl, FUN=function(df){rownames(df) <- NULL; df[, order(colnames(df))]})
    # rbind it all together
    output <- do.call("rbind", sl)
    rownames(output) <- NULL

    # if byFeature==TRUE, squish output
    if(byFeature){
      features <- unique(output$FeatureID)
      variables <- names(sl)
      spec_colnames <- paste0(variables, "_Spec")
      pval_colnames <- paste0(variables, "_Pval")
      # check those colnames are safe
      if(any(spec_colnames %in% colnames(output))){
        stop("combined variable_Spec column name already in use")
      }
      if(any(pval_colnames %in% colnames(output))){
        stop("combined variable_Pval column name already in use")
      }
      rowFromFeature <- function(f){
        specs <- sapply(X=variables, FUN=function(v){ output$Spec[output$FeatureID==f & output$Variable==v][1] })
        pvals <- sapply(X=variables, FUN=function(v){ output$Pval[output$FeatureID==f & output$Variable==v][1] })
        names(specs) <- spec_colnames
        names(pvals) <- pval_colnames
        cbind( data.frame(FeatureID=f, stringsAsFactors=FALSE), t(as.matrix( c(specs, pvals))) )
      }
      output <- do.call("rbind", lapply(X=features, FUN=rowFromFeature))
    }

    # if fd is present, check it out
    if(!is.null(fd)){
      # validate
      if(! class(fd) %in% c("matrix", "data.frame")){
        stop("fd is not a matrix or data.frame")
      }
      # get id column of fd
      if(is.numeric(fd_id) && fd_id %% 1 == 0){
        fd_id_col <- fd_id
      }else if(is.character(fd_id)){
        if(fd_id %in% colnames(fd)){
          fd_id_col <- which(colnames(fd) == fd_id)[1]
        }else{
          stop("fd_id string input not in colnames of fd")
        }
      }else{
        stop("invalid fd_id")
      }
      if(fd_id_col > ncol(fd)){
        stop("invalid fd_id")
      }
      # errors for reserved variables
      if(any(c("Spec", "Pval", "Variable") %in% colnames(fd))){
        stop("one or more reserved variable names (Spec, Pval, Variable) in fd colnames")
      }
      # remove fd[,fd_id_col] entries that are NA or ""
      fd <- fd[!is.na(fd[,fd_id_col]),]
      fd <- fd[fd[,fd_id_col] != "",]
      # check that fd[,fd_id_col] is unique
      if(length(fd[,fd_id_col]) != length(unique(fd[,fd_id_col]))){
        stop("fd_id column of fd contains duplicate entries")
      }

      # figure out order of fd[,fd_id_col] that matches output$FeatureID
      fd2out_order <- sapply(X=output$FeatureID, FUN=function(x){
        if(x %in% fd[,fd_id_col]){
          return(which(fd[,fd_id_col] == x))
        }else{
          return(NA)
        }
      })

      # for each column in fd except fd_id_col, move it to output using fd2out_order
      for(j in 1:ncol(fd)){if(j != fd_id_col){
        output <- data.frame(output, fd[,j][fd2out_order], stringsAsFactors=FALSE)
        colnames(output)[ncol(output)] <- colnames(fd)[j]
      }}

    }

    return(output)
  }

