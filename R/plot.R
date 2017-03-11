#' @title Correlation Plot
#'
#' @description This function easily draws the correlation plot of
#' surrogate variables (sv)
#' and variables.
#'
#' @param sv The matrix of sv identified by \code{\link{DaMiR.SV}}
#' function
#' @param df A data frame with class and known variables; at least
#' one column with
#' 'class' label must be included
#' @param type Type of correlation metric to be applied; default is
#' "pearson"
#' @param sig.level The significance level of the correlation; default
#'  is 0.0001
#'
#' @details Factorial variables are allowed. They will be tranformed as
#'  numeric before
#' applying the \code{\link{rcorr}} function of \code{Hmisc}.The
#' \code{\link{corrplot}}
#' function, which draws the plot, marks with a cross all the correlations
#'  that do not reach
#' the significance threshold defined in the \code{sig.level} argument.This
#'  plot allows the user to
#' identify those sv that present significant correlations with either
#' technical and
#' biological known variables.
#' Notably, none of the sv should present signifcant correlation with
#' "class" variable.
#'
#' @return
#' A correlation plot between sv and known variables.
#'
#' @author Mattia Chiesa, Luca Piacentini
#'
#' @examples
#' # use example data:
#' data(df)
#' data(sv)
#' # Draw correlation plot:
#' DaMiR.corrplot(sv=sv, df=df, type = "pearson", sig.level=0.01)
#'
#' @seealso
#' \code{\link{DaMiR.SV}}
#'
#' @export
#'
#'
DaMiR.corrplot <- function(sv,
                           df,
                           type=c("pearson","spearman"),
                           sig.level=0.0001){

  # check missing arguments
  if (missing(sv))
    stop("'sv' argument must be provided")
  if (missing(df))
    stop("'df' argument must be provided")
  if (missing(type)){
    type <- type[1]
  }

  # check the type of argument
  if(!(is.matrix(sv)))
    stop("'sv' must be a matrix")
  if(!(is(df, "DataFrame") | is.data.frame(df)))
    stop("'df' must be a data.frame")
  if(!(is.numeric(sig.level)))
    stop("'sig.level' must be numeric")

  # check the presence of NA or Inf
  df<-as.data.frame(df)
  if (any(is.na(sv)))
    stop("NA values are not allowed in the 'sv' argument")
  if (any(is.na(df)))
    stop("NA values are not allowed in the 'df' matrix")
  if (any(is.infinite(sv)))
    stop("Inf values are not allowed in the 'sv' matrix")

  # specific checks
  if(!("class" %in% colnames(df)))
    stop("'class' info is lacking!
         Include the variable 'class'
         in the 'df' data frame and label it 'class'!")
  if(dim(sv)[1] != dim(df)[1])
    stop("nrow(df) must be equal to nrow(sv)")
  if (sig.level >1 | sig.level < 0)
    stop("'sig.level' must be between 0 and 1")

  sva_corr <- as.data.frame(cbind(sv,df))
  for (i in seq_len(ncol(sva_corr))){
    if (is.numeric(sva_corr[,i])==FALSE){
      droplevels(sva_corr[,i])
      sva_corr[,i] <- as.numeric(sva_corr[,i])
    }
  }

  if(type == "pearson"){
    corrplot::corrplot(rcorr(as.matrix(sva_corr), type='pearson')$r,
                       type = "upper",
                       is.corr = TRUE,
                       addCoef.col = "black",
                       number.cex=0.5,
                       p.mat = rcorr(as.matrix(sva_corr),
                                     type='pearson')$P,
                       sig.level = sig.level)
  }else if(type == "spearman"){
    corrplot::corrplot(rcorr(as.matrix(sva_corr), type='spearman')$r,
                       type = "upper",
                       is.corr = TRUE,
                       addCoef.col = "black",
                       number.cex=0.5,
                       p.mat = rcorr(as.matrix(sva_corr),
                                     type='spearman')$P,
                       sig.level = sig.level)
  } else{
    stop("Please set 'spearman or 'pearson' as correlation type.")
  }
}


#' @title Quality assessment and visualization of expression data
#' @description This is a helper function to easily draw (1) clustering
#'  dendrogram and
#' heatmap of a sample-per-sample correlation matrix, (2)
#' multidimensional scaling plots (MDS)  and (3) relative log expression
#'  (RLE) boxplots
#' of expression data.
#'
#' @param data A SummarizedExperiment object or an expression data frame
#' @param df A data frame with class and known variables (or a subset
#' of them); at least one column with
#' 'class' label must be included
#' @param type A character string specifing the metric to be applied to
#'  correlation
#' analysis. Either "spearman" or "pearson" is allowed; default is
#' "spearman"
#'
#' @details
#' Please be sure that NAs are not present in \code{df}'s columns.
#' Plots will not be drawn in the presence of NAs.
#' @return
#' A dendrogram and heatmap, MDS plot(s) and a RLE boxplot
#'
#' @author Mattia Chiesa, Luca Piacentini
#'
#' @examples
#' # use example data:
#' data(data_norm)
#' data(df)
#' # Draw clustering dendrogram and heatmap, MDS, RLE boxplot:
#' DaMiR.Allplot(data=data_norm, df=df)
#'
#' @export
#'
#'
DaMiR.Allplot <- function(data,
                          df,
                          type=c("spearman","pearson")){

  # check arguments
  if (missing(data))
    stop("'data' argument must be provided")
  if (missing(df))
    stop("'df' argument must be provided")
  if (missing(type)){
    type <- type[1]
  }

  # check the type of argument
  if(!(is(data, "SummarizedExperiment") | is.data.frame(data)))
    stop("'data' must be a 'data.frame' or a 'SummarizedExperiment' object")
  if(!(is(df, "DataFrame") | is.data.frame(df)))
    stop("'df' must be a data.frame")

  if (is(data, "SummarizedExperiment")) {
    count_data<-assay(data)
  } else {
    count_data<-as.matrix(data)
  }
  df<-as.data.frame(df)

  # check the presence of NA or Inf
  if (any(is.na(count_data)))
    stop("NA values are not allowed in the 'data' matrix")
  if (any(is.na(df)))
    stop("NA values are not allowed in the 'df' matrix")
  if (any(is.infinite(count_data)))
    stop("Inf values are not allowed in the 'data' matrix")

  # Specific checks
  if(!("class" %in% colnames(df)))
    stop("'class' info is lacking!
         Include the variable 'class'
         in the 'df' data frame and label it 'class'!")
  if (all((count_data %%1) == 0))
    stop("This function works with normalized count data")
  if(dim(count_data)[2] != dim(df)[1])
    stop("ncol(assay(data)) must be equal to nrow(df)")



  colors <- colorRampPalette(rev(brewer.pal(9, "RdYlGn")) )(255)
  #############
  ## Heatmap
  if(type == "spearman"){
    mydist<-as.dist(1-cor(count_data, method="spearman"))}
  else if(type == "pearson"){
    mydist<-as.dist(1-cor(count_data, method="pearson"))
  } else{
    stop("Please set 'spearman or 'pearson' as correlation type.")
  }

  sampleDistMatrix <- as.matrix(mydist)
  seque <- seq(min(sampleDistMatrix),
               max(sampleDistMatrix),
               by=max(sampleDistMatrix)/100)

  pheatmap(sampleDistMatrix,
           clustering_distance_rows=mydist,
           clustering_distance_cols=mydist,
           col=colors,
           breaks = seque,
           annotation_col = df)

  ################
  ## MDS plot
  mdsData <- data.frame(cmdscale(sampleDistMatrix))
  mds <- cbind(mdsData, df)

  for(i in seq_len(ncol(df))) {
    cov_list <- mds[,i+2,drop=FALSE]
    colnames(cov_list)<-"Vars"

    print(ggplot(mds, aes(X1, X2, shape=mds$class, color=cov_list$Vars)) +
            geom_point(size=3) +
            geom_text(aes(label=rownames(mds)),hjust=0.5,vjust=-1) +
            ggtitle(paste("Variable: ",colnames(mds[,i+2,drop=FALSE]))))
  }
  ################
  ## RLE
  colors <- brewer.pal(3, "Set2")

  plotRLE(count_data,
          k=2,
          labels=TRUE,
          isLog=TRUE,
          outline=FALSE,
          col=colors[df$class],
          main="Relative Log Expression")
}

#' @title Plot multidimentional scaling (MDS)
#'
#' @description A MDS plot is drawn in order to visualize class clustering.
#'
#' @param data A SummarizedExperiment object or an expression data.frame
#' @param df A data frame with class; it can be directly subset from data
#' @param type A character string specifing the metric to be applied to
#' correlation
#' analysis. Either "spearman" or "pearson" is allowed; default is
#' "spearman"
#'
#' @details The MDS plot is drawn taking as input a dissimilarity matrix
#' produced by either
#' a sample-per-sample Pearson's or Spearman's correlation of normalized
#' expression data.
#' @return
#' A MDS plot, using only 'class' information
#' @author Mattia Chiesa, Luca Piacentini
#'
#' @examples
#' # use example data:
#' data(data_reduced)
#' data(df)
#' # Draw MDS:
#' DaMiR.MDSplot(data=data_reduced, df=df, type="pearson")
#'
#' @export
#'
#'
DaMiR.MDSplot <- function(data,
                          df,
                          type=c("spearman","pearson")){

  # check arguments
  if (missing(data))
    stop("'data' argument must be provided")
  if (missing(df))
    stop("'df' argument must be provided")
  if (missing(type)){
    type <- type[1]
  }

  # check the type of argument
  if(!(is(data, "SummarizedExperiment") | is.data.frame(data)))
    stop("'data' must be a 'data.frame' or a 'SummarizedExperiment' object")
  if(!(is(df, "DataFrame") | is.data.frame(df)))
    stop("'df' must be a data.frame")

  if (is(data, "SummarizedExperiment")) {
    t_count_data<-t(assay(data))
  } else {
    t_count_data<-t(as.matrix(data))
  }
  df<-as.data.frame(df)

  # check the presence of NA or Inf
  if (any(is.na(t_count_data)))
    stop("NA values are not allowed in the 'data' matrix")
  if (any(is.na(df)))
    stop("NA values are not allowed in the 'df' matrix")
  if (any(is.infinite(t_count_data)))
    stop("Inf values are not allowed in the 'data' matrix")

  # Specific checks
  if(!("class" %in% colnames(df)))
    stop("'class' info is lacking!
         Include the variable 'class'
         in the 'df' data frame and label it 'class'!")
  if (all((t_count_data %%1) == 0))
    stop("This function works with normalized count data")
  if(dim(t_count_data)[2] != dim(df)[1])
    stop("ncol(assay(data)) must be equal to nrow(df)")


  sample_list <- colnames(t_count_data)
  colors <- colorRampPalette(rev(brewer.pal(9, "RdYlGn")) )(255)

  if(type == "spearman"){
    mydist<-as.dist(1-cor(t_count_data, method="spearman"))}
  else if(type == "pearson"){
    mydist<-as.dist(1-cor(t_count_data, method="pearson"))
  } else{
    stop("Please set 'spearman or 'pearson' as correlation type.")
  }

  sampleDistMatrix <- as.matrix(mydist)
  colnames(sampleDistMatrix)<-sample_list
  rownames(sampleDistMatrix)<-sample_list

  # MDS plot
  mdsData <- data.frame(cmdscale(sampleDistMatrix))
  mds <- cbind(mdsData, df)
  ggplot(mds, aes(X1, X2, shape=df$class, color=df$class)) +
    geom_point(size=3) +
    geom_text(aes(label=rownames(mds)), hjust=0.5, vjust=-1)

}

#' @title Expression data clustering and heatmap
#'
#' @description The function helps to draw a clustering dendrogram and
#'  a heatmap of expression data.
#'
#' @param data A SummarizedExpression object, or an expression data frame
#'  where
#' rows and cols should be, respectively, observations and features
#'
#' @param df A data frame with class and (optionally) known variables;
#' at least one column with
#' 'class' label must be included
#' @param type_row The metric to be used to cluster rows. Either
#' "euclidean" or "correlation"
#'  is allowed; default is "euclidean"
#' @param type_col The metric to be used to cluster cols. Either
#' "euclidean" or "correlation"
#'  is allowed; default is "euclidean"
#'
#' @return
#' A clustering dendrogram and heatmap.
#'
#' @author Mattia Chiesa, Luca Piacentini
#'
#' @examples
#' # use example data:
#' data(data_norm)
#' data(df)
#' # use the first 100 genes:
#' data_norm_red<-data_norm[1:100,]
#' # Draw heatmap: samples (cols) per genes (rows)
#' # and use variable annotation:
#' DaMiR.Clustplot(data=data_norm_red,
#' df=df, type_row="correlation", type_col="correlation")
#'
#' @export
#'
DaMiR.Clustplot <- function(data,
                            df,
                            type_row=c("euclidean","correlation"),
                            type_col=c("euclidean","correlation")){

  # check missing arguments
  if (missing(data))
    stop("'data' argument must be provided")
  if (missing(df))
    stop("'df' argument must be provided")
  if (missing(type_row)){
    type_row <- type_row[1]
  }
  if (missing(type_col)){
    type_col <- type_col[1]
  }

  # check the type of argument
  if(!(is(data, "SummarizedExperiment") | is.data.frame(data)))
    stop("'data' must be a 'data.frame' or a 'SummarizedExperiment' object")
  if(!(is(df, "DataFrame") | is.data.frame(df)))
    stop("'df' must be a data.frame")

  if (is(data, "SummarizedExperiment")) {
    count_data<-assay(data)
  } else {
    count_data<-t(data)
  }
  df<-as.data.frame(df)

  # check the presence of NA or Inf
  if (any(is.na(count_data)))
    stop("NA values are not allowed in the 'data' matrix")
  if (any(is.na(df)))
    stop("NA values are not allowed in the 'df' matrix")
  if (any(is.infinite(count_data)))
    stop("Inf values are not allowed in the 'data' matrix")

  # Specific checks
  if(!("class" %in% colnames(df)))
    stop("'class' info is lacking!
         Include the variable 'class'
         in the 'df' data frame and label it 'class'!")
  if (all((count_data %%1) == 0))
    stop("This function works with normalized count data")
  if(dim(count_data)[2] != dim(df)[1])
    stop("ncol(assay(data)) must be equal to nrow(df)")


  if(type_row == "euclidean"){d_r <- "euclidean"}
  else if(type_row == "correlation"){d_r <- "correlation"}
  else {stop("Please set 'euclidean' or 'correlation' as distance type.")}

  if(type_col == "euclidean"){d_c <- "euclidean"}
  else if(type_col == "correlation"){d_c <- "correlation"}
  else {stop("Please set 'euclidean' or 'correlation' as distance type.")}

  # Dendrogram
  colors <- colorRampPalette(rev(brewer.pal(9, "RdYlGn")) )(255)
  pheatmap(count_data,
           clustering_distance_rows = d_r,
           clustering_distance_cols = d_c,
           scale = "row",
           col=colors,
           annotation_col = df)
}




