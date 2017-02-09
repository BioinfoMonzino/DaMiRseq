#' @title Quality assessment and visualization of expression data
#' @description This is a helper function to easily draw (1) clustering dendrogram and
#' heatmap of a sample-per-sample correlation matrix, (2)
#' multidimensional scaling plots (MDS)  and (3) relative log expression (RLE) boxplots
#' of expression data.
#'
#' @param data Matrix of normalized expression data, i.e transformed counts by vst or rlog
#' @param df A data frame with class and known covariates; at least one column with
#' 'class' label must be included
#' @param type A character string specifing the metric to be applied to correlation
#' analysis. Either "spearman" or "pearson" is allowed; default is "spearman"
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
#' data(DaMiR.Data.Ex.Norm)
#' data(DaMiR.Data.Ex.dfClCov)
#' DaMiR.Allplot(DaMiR.Data.Ex.Norm, DaMiR.Data.Ex.dfClCov)
#' DaMiR.Allplot(DaMiR.Data.Ex.Norm, DaMiR.Data.Ex.dfClCov, type = "pearson")
#'
#' @export
#'
#'
DaMiR.Allplot <- function(data, df, type=c("spearman","pearson")){

  # check arguments
  if (missing(data)) stop("'data' argument must be provided")
  if (missing(df)) stop("'df' argument must be provided")
  if (missing(type)){
    type <- type[1]
  }

  if(!(is.matrix(data))) stop("'data' must be a matrix")
  if(!(is.data.frame(df))) stop("'df' must be a data frame")

  colors <- colorRampPalette(rev(brewer.pal(9, "RdYlGn")) )(255)
  #############
  ## Heatmap
  if(type == "spearman"){
    mydist<-as.dist(1-cor(data, method="spearman"))} # use dissimilarity matrix
  else if(type == "pearson"){
    mydist<-as.dist(1-cor(data, method="pearson")) # use dissimilarity matrix
  } else{
    stop("Please set 'spearman or 'pearson' as correlation type.")
  }

  sampleDistMatrix <- as.matrix(mydist)
  seque <- seq(min(sampleDistMatrix),max(sampleDistMatrix),by=max(sampleDistMatrix)/100)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=mydist,
           clustering_distance_cols=mydist,
           col=colors,
           breaks = seque,
           annotation_col = df)

  ################
  ## MDS plot
  mdsData <- data.frame(cmdscale(sampleDistMatrix)) # by dissimilarity matrix
  mds <- cbind(mdsData, df)

  for(i in 1:ncol(df)) {
    cov_list <- mds[,i+2,drop=FALSE]
    colnames(cov_list)<-"Vars"

    print(ggplot(mds, aes(X1,X2,shape=mds$class, color=cov_list$Vars)) + geom_point(size=3) + geom_text(aes(label=rownames(mds)),hjust=0.5,vjust=-1) + ggtitle(paste("Variable: ",colnames(mds[,i+2,drop=FALSE]))))
  }
  ################
  ## RLE
  colors <- brewer.pal(3, "Set2")

  plotRLE(data, k=2, labels=TRUE, isLog=TRUE, outline=FALSE, col=colors[df$class], main="Relative Log Expression")
}

#' @title Plot multidimentional scaling (MDS)
#'
#' @description A MDS plot is drawn in order to visualize class clustering.
#'
#' @param data A data frame or matrix of normalized expression data, i.e. transformed counts by vst or rlog
#' @param df A data frame with class and known covariates; at least one column with
#' 'class' label must be included
#' @param type A character string specifing the metric to be applied to correlation
#' analysis. Either "spearman" or "pearson" is allowed; default is "spearman"
#'
#' @details The MDS plot is drawn taking as input a dissimilarity matrix produced by either
#' a sample-per-sample Pearson's or Spearman's correlation of normalized expression data.
#' @return
#' A MDS plot, using only 'class' information
#' @author Mattia Chiesa, Luca Piacentini
#'
#' @examples
#' data(DaMiR.Data.Ex.t.Norm)
#' data(DaMiR.Data.Ex.dfClCov)
#' DaMiR.MDSplot(DaMiR.Data.Ex.t.Norm, DaMiR.Data.Ex.dfClCov)
#'
#' @export
#'
#'
DaMiR.MDSplot <- function(data, df, type=c("spearman","pearson")){

  # check arguments
  if (missing(data)) stop("'data' argument must be provided")
  if (missing(df)) stop("'df' argument must be provided")

  if(!(is.data.frame(df))) stop("'df' must be a data frame")

  t_data<-as.matrix(t(data))
  sample_list <- colnames(t_data)
  colors <- colorRampPalette(rev(brewer.pal(9, "RdYlGn")) )(255)
  if (missing(type)){
    type <- type[1]
  }
  if(type == "spearman"){
    mydist<-as.dist(1-cor(t_data, method="spearman"))} # use dissimilarity matrix
  else if(type == "pearson"){
    mydist<-as.dist(1-cor(t_data, method="pearson")) # use dissimilarity matrix
  } else{
    stop("Please set 'spearman or 'pearson' as correlation type.")
  }

  sampleDistMatrix <- as.matrix(mydist)
  colnames(sampleDistMatrix)<-sample_list
  rownames(sampleDistMatrix)<-sample_list

  # MDS plot
  mdsData <- data.frame(cmdscale(sampleDistMatrix)) # by dissimilarity matrix
  mds <- cbind(mdsData, df)
  ggplot(mds, aes(X1,X2,shape=df$class, color=df$class)) + geom_point(size=3) + geom_text(aes(label=rownames(mds)),hjust=0.5,vjust=-1)

}

#' @title Expression data clustering and heatmap
#'
#' @description The function helps to draw a clustering dendrogram and a heatmap of expression data.
#'
#' @param data A data frame or matrix of normalized expression data, i.e transformed counts by vst or rlog
#' @param df A data frame with class and known covariates; at least one column with
#' 'class' label must be included
#' @param type_row The metric to be used to cluster rows. Either "euclidean" or "correlation"
#'  is allowed; default is "euclidean"
#' @param type_col The metric to be used to cluster cols. Either "euclidean" or "correlation"
#'  is allowed; default is "euclidean"
#'
#' @return
#' A clustering dendrogram and heatmap.
#'
#' @author Mattia Chiesa, Luca Piacentini
#'
#' @examples
#' data(DaMiR.Data.Ex.t.Norm)
#' data(DaMiR.Data.Ex.dfClCov)
#' DaMiR.Clustplot(DaMiR.Data.Ex.t.Norm, DaMiR.Data.Ex.dfClCov)
#'
#' @export
#'
DaMiR.Clustplot <- function(data, df, type_row=c("euclidean","correlation"),type_col=c("euclidean","correlation")){
  if (missing(data)) stop("'data' argument must be provided")
  if (missing(df)) stop("'df' argument must be provided")
  if (missing(type_row)){
    type_row <- type_row[1]
  }
  if (missing(type_col)){
    type_col <- type_col[1]
  }

  if(!(is.data.frame(df))) stop("'df' must be a data frame")

  if(type_row == "euclidean"){d_r <- "euclidean"}
  else if(type_row == "correlation"){d_r <- "correlation"}
  else{stop("Please set 'euclidean' or 'correlation' as distance type.")}

  if(type_col == "euclidean"){d_c <- "euclidean"}
  else if(type_col == "correlation"){d_c <- "correlation"}
  else{stop("Please set 'euclidean' or 'correlation' as distance type.")}

  # Dendrogram
  colors <- colorRampPalette(rev(brewer.pal(9, "RdYlGn")) )(255)
  pheatmap(t(data),
           clustering_distance_rows = d_r,
           clustering_distance_cols = d_c,
           scale = "row",
           col=colors,
           annotation_col = df)
}




