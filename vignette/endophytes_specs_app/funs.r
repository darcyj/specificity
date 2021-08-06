#' aggregate_specs_list
#'
#' Aggregates a list of outputs from specificity::phy_or_env_spec() into a single 
#' data.frame object. Can also include feature data (e.g. species taxonomy) into 
#' that output. Output can also be byFeature, with one row per feature, and multiple
#' columns for different variables.
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
#' # attach(antarctica)
#' # # to see details of the antarctica data set, see ?antarctica
#' # # aggregate data long format:
#' # ant_long <- aggregate_specs_list(antarctica_specs_list, byFeature=FALSE,
#' #   fd=antarctica_taxonomy, fd_id=1)
#' # head(ant_long)
#' # ant_wide <- aggregate_specs_list(antarctica_specs_list, byFeature=TRUE,
#' #   fd=antarctica_taxonomy, fd_id=1)
#' # head(ant_wide)
#' # # aggregation can also be done without feature data (fd):
#' # ant_long <- aggregate_specs_list(antarctica_specs_list, byFeature=FALSE)
#' # head(ant_long)
#' # ant_wide <- aggregate_specs_list(antarctica_specs_list, byFeature=TRUE)
#' # head(ant_wide)
#' 
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

#' plot_specs_shiny
#' 
#' Runs an interactive shiny visualization. Data can be re-downloaded from that
#' visualization, in either "wide" or "long" format. See ?aggregate_specs_list 
#' for explanation of those formats.
#' 
#' @author John L. Darcy
#'
#' @param sl specs_list. A named list of outputs from phy_or_env_spec. See examples.
#' @param fd data.frame. Optional feature data - a data.frame object with one row per 
#'   feature, including some column with feature IDs that includes feature IDs in sl as
#'   rownames. If NULL, no feature data will be used (default:NULL).
#' @param fd_id_col integer or string. Either the column number of fd containing feature
#'   IDs (see above), or the column name of that column (default: 1).
#' @return Returns output of shiny::shinyApp(). If run locally, this will manifest
#'   as a browser window opening with an interactive visualization. 
#'
#' @examples
#' # attach(antarctica)
#' # plot_specs_shiny(antarctica_specs_list, antarctica_taxonomy, 1)
#' 
#' 
#' @export
plot_specs_shiny <- function(sl, fd=NULL, fd_id_col=1){
  # handle fd_id_col, but only if fd is provided
  if(!is.null(fd)){
    if(is.numeric(fd_id_col) && (fd_id_col %% 1 == 0) && (fd_id_col <= ncol(fd)) && (fd_id_col > 0)){
      fd_id_col_name <- colnames(fd)[fd_id_col]
    }else if(is.character(fd_id_col) && fd_id_col %in% colnames(fd)){
      fd_id_col_name <- fd_id_col
    }else{
      stop("invalid fd_id_col")
    }
  }

  # aggregate data into long and wide data types
  # long - used for plotting with ggplot, does NOT have user metadata
  # fat - used for displaying selected data to user
  # baselong - only used for download outputs; contains everything
  longdata <- aggregate_specs_list(sl, byFeature=FALSE)
  widedata <- aggregate_specs_list(sl, byFeature=TRUE, fd=fd, 
    fd_id=fd_id_col_name)
  baselongdata <- aggregate_specs_list(sl, byFeature=FALSE, fd=fd,
    fd_id=fd_id_col_name)

  # figure out what columns are available for searching
  fieldChoices <- colnames(widedata)
  searchdefault <- "FeatureID"

  # add plotting helper variables to longdata
  # jittercenter: the x-axis center of violins
  # jitterx : the x-axis position of jittered points
  # This manual jittering has to be done so that the base ggplot type can be
  # a scatterplot instead of violin, since violin won't work with brushing.
  # said another way, scales for x and y must be continuous.
  longdata$jitterx <- longdata$jittercenter <- rep(0, nrow(longdata))
  jitterwidth <- 0.6 # (1 means jittered columns are touching)
  vars <- unique(longdata$Variable)
  for(i in 1:length(vars)){
    var_i <- vars[i]
    longdata$jitterx[longdata$Variable==var_i] <- runif(
      n=sum(longdata$Variable==var_i), 
      min=i-(jitterwidth/2), 
      max=i+(jitterwidth/2)
    )
    longdata$jittercenter[longdata$Variable==var_i] <- i
  }
  rm(var_i)
  # hilight: whether or not a point is hilighted. 0=don't hilight, 1=hilight.
  # to be interpreted as alpha (opacity) value by ggplot.
  longdata$hilight01 <- rep(0, nrow(longdata))
  # show01: show this point, true or false? actually interpreted in ggplot as
  # alpha, i.e. opacity. 1=show, 0=hide. points will be hidden if not found by
  # search function.
  longdata$show01 <- rep(1, nrow(longdata))

  # define download choices
  downloadChoices <- c(
    "All data (long)", 
    "All data (wide)",
    "Shown/searched data (long)",
    "Shown/searched data (wide)",
    "Brushed data (long)",
    "Brushed data (wide)"
  )

  # page setup stuff
  ui <- shiny::fluidPage(
    shiny::fluidRow(
      shiny::column(width = 12,
        shiny::plotOutput("plot1", height = 300,
          brush = shiny::brushOpts(id = "plot1_brush")
        )
      )
    ),
    shiny::fluidRow(
      shiny::column(width=4, 
        shiny::textInput("searchQuery", label = shiny::h4("Search string:"), value = ""),
      ),
      shiny::column(width=4, 
        shiny::selectInput("searchField", label = shiny::h4("Search field:"), 
          choices=fieldChoices, selected=searchdefault)
      ),
      shiny::column(width=2, 
        shiny::selectInput("searchIgnoreCase", label = shiny::h4("Ignore Case?"),
          choices=c(TRUE,FALSE), selected=TRUE)
      ),
      shiny::tags$style(type='text/css', "#searchButtonClick { width:100%; margin-top: 45px;}")

    ),
    # these are formatted as such because I had difficulty getting the parens right
    shiny::fluidRow(
      shiny::column(width=2,
        colourpicker::colourInput("sigColor", 
          value="black",
          label=shiny::h4(
            shiny::HTML(
              paste0(
                "<em>P</em>", 
                intToUtf8(8804), 
                intToUtf8(945), 
                " color:"
              )
            )
          )
        )
      ),
      shiny::column(width=2, 
        colourpicker::colourInput("nsigColor", 
          value="gray",
          label=shiny::h4(
            shiny::HTML(
              paste0(
                "<em>P</em>", 
                ">", 
                intToUtf8(945), 
                " color:"
              )
            )
          )
        )
      ),
      shiny::column(width=2, 
        shiny::numericInput("alpha", label=shiny::h4(paste0(intToUtf8(945), ":")), value=0.05, step=0.01),
      ),
      shiny::column(width=2, 
        shiny::numericInput("pointSize", label = shiny::h4("Point size:"), value=0.5, step=0.1),
      ),
      shiny::column(width=2, 
        shiny::selectInput("pointPCH", label = shiny::h4("Point type:"),
          choices=as.character(0:20), selected="19")
      )
    ),

    shiny::fluidRow(
      shiny::column(width=2,
        colourpicker::colourInput("hilightColor", label=shiny::h4("Brush hilight color:"),
          value="red")
      ),
      shiny::column(width=2, 
        shiny::numericInput("hilightSize", label=shiny::h4("Brush hilight size:"),
          value=2, step=0.1),
      ),
      shiny::column(width=4, 
        shiny::selectInput("downloadType", label = shiny::h4("Download type:"),
          choices=downloadChoices, selected=1)
      ),
      shiny::column(width=2,
        shiny::downloadButton("downloadData", "Download")
      ),
      shiny::tags$style(type='text/css', "#downloadData { width:100%; margin-top: 45px;}")
    ),

    shiny::fluidRow(
      shiny::column(width=12,
        shiny::checkboxGroupInput(inputId="displayFields",
          label=shiny::h4("Fields to show in brushed points (below):"),
          choices = fieldChoices, selected=fieldChoices, inline=TRUE)
      )
    ),

    shiny::fluidRow(
      shiny::column(width = 12,
        shiny::h4("Brushed points:"),
        # shiny::verbatimTextOutput("brush_info") # OLD
        DT::dataTableOutput("brush_info") # NEW
      )
    )
  )
  # function to build points_df, the data.frame used for plotting geom_point via ggplot
  make_points_df <- function(
    longdf,             # complete data object, long format (longdata)
    widedf,             # complete data object, wide format (widedata)
    searchField,        # which field to search
    searchQuery,        # text, what to search for
    searchIgnoreCase,   # if true, ignore case
    sigColor,           # color for significant features
    nsigColor,          # color for non-significant features
    alpha){             # alpha value to determine significance (e.g. 0.05)
    
    # longdf <- longdata
    # widedf <- widedata
    # searchField <- "tax"
    # searchQuery <- "cyano"
    # searchIgnoreCase <- TRUE
    # sigColor <- "red"
    # nsigColor <- "black"
    # alpha <- 0.05

    # determine significance, assign colors
    longdf$pointColor <- rep(nsigColor, nrow(longdf))
    longdf$pointColor[longdf$Pval <= alpha] <- sigColor

    # do search, unless blank.
    searchQuery <- gsub(" ", "", searchQuery)
    if(searchQuery != ""){
      j <- which(colnames(widedf) == searchField)[1]
      hits <- grepl(pattern=searchQuery, x=as.character(widedf[,j]), 
        ignore.case=searchIgnoreCase)
      hitIDs <- widedf$FeatureID[hits]
      longdf$show01 <- rep(0, nrow(longdf))
      longdf$show01[longdf$FeatureID %in% hitIDs] <- 1
    }else{
      longdf$show01 <- rep(1, nrow(longdf))
    }

    # significance
    longdf$sigTF <- longdata$Pval <= alpha

    # subset to exclude unshown points
    return( longdf[longdf$show01 == 1, ] )
  }

  # function to make hilight_df, which is plotted UNDER points_df for hilight effect
  make_hilight_df <- function(
    points_df,     # from make_points_df
    brush_df){     # from brushedPoints()

    return( points_df[points_df$FeatureID %in% brush_df$FeatureID, ] )
  }

  # function to build display_df (rows shown for selected points)
  make_display_df <- function(
    hilight_df,        # see above
    widedf,            # widedata full version
    displayFields){    # character vector of which fields to show

    # subset widedf to only include visible+hilighted features in plotdf
    features2keep <- unique(hilight_df$FeatureID)
    output <- widedf[widedf$FeatureID %in% features2keep, ]

    # subset widedf to only include columns selected by user
    output <- output[, colnames(output) %in% displayFields]

    return(output)
  }

  # function to extract colors from points_df since ggplot is STUPID and can't take 
  # colors as a column inside aes()
  get_colors_from_points_df <- function(points_df){points_df$pointColor}
  get_FeatureIDs <- function(df){unique(df$FeatureID)}

  # server function
  server <- function(input, output, session) {

    # build points_df
    points_df <- shiny::reactive({make_points_df(
      longdf=longdata,
      widedf=widedata,
      searchField=input$searchField,
      searchQuery=input$searchQuery,
      searchIgnoreCase=input$searchIgnoreCase,
      sigColor=input$sigColor,
      nsigColor=input$nsigColor,
      alpha=input$alpha
    )})

    # build hilight_df
    hilight_df <- shiny::reactive({make_hilight_df(
      points_df = points_df(),
      brush_df = shiny::brushedPoints(longdata, input$plot1_brush)
    )})

    # build display_df
    display_df <- shiny::reactive({make_display_df(
      hilight_df=hilight_df(),
      widedf=widedata,
      displayFields=input$displayFields
    )})

    # output brushdata for rendering
    # output$brush_info <- shiny::renderPrint({ display_df() }) # OLD
    output$brush_info <- DT::renderDataTable(display_df(), 
      options = list(scrollX=TRUE, sDom='<"top">lrt<"bottom">ip')) # NEW

    # draw plot
    # using aes_string() instead of aes() because R CMD check thinks that the unquoted
    # variable names used by aes() are undeclared variables. Easy enough to use aes_string()
    # instead, which allows the variable names to be quoted.
    output$plot1 <- shiny::renderPlot({
      ggplot2::ggplot(data=longdata, mapping=ggplot2::aes_string(x="jitterx", y="Spec")) + 
        ggplot2::scale_x_continuous(limits=c(0.5, length(vars)+0.5), breaks=1:length(vars), 
          labels=vars, expand = c(0, 0)) + ggplot2::theme(axis.title.x = ggplot2::element_blank()) +
        ggplot2::geom_violin(ggplot2::aes_string(group="jittercenter", y="Spec"), scale="count" ) + 

        ggplot2::geom_point(data=hilight_df(), ggplot2::aes_string(x="jitterx", y="Spec"), 
          color=input$hilightColor, size=input$hilightSize, shape=as.integer(input$pointPCH)) +

        ggplot2::geom_point(data=points_df(), ggplot2::aes_string(x="jitterx", y="Spec"), 
          color=get_colors_from_points_df(points_df()), 
          size=input$pointSize, shape=as.integer(input$pointPCH))

    })

    # make data to download
    data2download <- shiny::reactive({
      if(input$downloadType == downloadChoices[1]){ #"All data (long)"
        baselongdata
      }else if(input$downloadType == downloadChoices[2]){ #"All data (wide)"
        widedata
      }else if(input$downloadType == downloadChoices[3]){ #"Shown/searched data (long)"
        baselongdata[baselongdata$FeatureID %in% get_FeatureIDs(points_df()),]
      }else if(input$downloadType == downloadChoices[4]){ #"Shown/searched data (wide)"
        widedata[widedata$FeatureID %in% get_FeatureIDs(points_df()),]
      }else if(input$downloadType == downloadChoices[5]){ #"Brushed data (long)"
        baselongdata[baselongdata$FeatureID %in% get_FeatureIDs(display_df()),]
      }else if(input$downloadType == downloadChoices[6]){ #"Brushed data (wide)"
        widedata[widedata$FeatureID %in% get_FeatureIDs(display_df()),]
      }else{
        # just do full long data
        baselongdata
      }
    })

    output$downloadData <- shiny::downloadHandler(
      filename = function() {
        filename <- gsub(pattern="[ \\(\\)\\/]+", replacement="_", x=input$downloadType)
        filename <- sub(pattern="_+$", replacement="", x=filename)
        paste(filename, "_", Sys.Date(), ".tsv", sep="")
      },
      content = function(file) {
        write.table(x=data2download(), file=file, sep="\t", row.names=FALSE, quote=FALSE)
      }
    )

  }

  # run program
  shiny::shinyApp(ui, server)
}


