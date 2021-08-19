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
#' #DaMiR.corrplot(sv=sv, df=df, type = "pearson", sig.level=0.01)
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
                           sig.level=0.01){

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
  if (length(type) > 1)
    stop("length(type) must be equal to 1")
  if (!(all(type %in% c("pearson", "spearman"))))
    stop("'type' must be 'pearson' or 'spearman'")

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
                       sig.level = sig.level,
                       col = colorRampPalette(rev(c("#67001F",
                                                    "#B2182B",
                                                    "#D6604D",
                                                    "#F4A582",
                                                    "#FDDBC7",
                                                    "#FFFFFF",
                                                    "#D1E5F0",
                                                    "#92C5DE",
                                                    "#4393C3",
                                                    "#2166AC",
                                                    "#053061")))(200))
  }else if(type == "spearman"){
    corrplot::corrplot(rcorr(as.matrix(sva_corr), type='spearman')$r,
                       type = "upper",
                       is.corr = TRUE,
                       addCoef.col = "black",
                       number.cex=0.5,
                       p.mat = rcorr(as.matrix(sva_corr),
                                     type='spearman')$P,
                       sig.level = sig.level,
                       col = colorRampPalette(rev(c("#67001F",
                                                    "#B2182B",
                                                    "#D6604D",
                                                    "#F4A582",
                                                    "#FDDBC7",
                                                    "#FFFFFF",
                                                    "#D1E5F0",
                                                    "#92C5DE",
                                                    "#4393C3",
                                                    "#2166AC",
                                                    "#053061")))(200))
  } else{
    stop("Please set 'spearman or 'pearson' as correlation type.")
  }
}


#' @title Quality assessment and visualization of expression data
#' @description This is a helper function to easily draw (1) clustering
#'  dendrogram and heatmap of a sample-per-sample correlation matrix, (2)
#' multidimensional scaling plots (MDS), (3) relative log expression
#'  (RLE) boxplots of expression data, (4) a sample-by-sample expression
#'  value distribution, and (5) a class average expression value
#'  distribution
#'
#' @param data A SummarizedExperiment object or a matrix or a data.frame
#'  where rows and cols should be, respectively, observations and features
#' @param df A data frame with class and known variables (or a subset
#' of them); at least one column with
#' 'class' label must be included
#' @param type A character string specifing the metric to be applied to
#'  correlation
#' analysis. Either "spearman" or "pearson" is allowed; default is
#' "spearman"
#' @param what A character string specifing the plots to be shown
#'  'all', 'all_w_PCA', 'MDS','PCA','heatmap','RLEbox', 'distr',
#'  'avg_distr' are allowed; default is "all"
#'
#' @details
#' Please be sure that NAs are not present in \code{df}'s columns.
#' Plots will not be drawn in the presence of NAs.
#' @return
#' A dendrogram and heatmap, MDS plot(s), a RLE boxplot, a
#' sample-by-sample expression value distribution, and a class
#' average expression value distribution
#'
#' @author Mattia Chiesa, Luca Piacentini
#'
#' @examples
#' # use example data:
#' data(data_norm)
#' data(df)
#' # Draw clustering dendrogram and heatmap, MDS, RLE boxplot:
#' DaMiR.Allplot(data=data_norm, df=df[,5,drop=FALSE])
#'
#' @export
#'
#'
DaMiR.Allplot <- function(data,
                          df,
                          type=c("spearman","pearson"),
                          what=c("all",
                                 "all_w_PCA",
                                 "MDS",
                                 "PCA",
                                 "heatmap",
                                 "RLEbox",
                                 "distr",
                                 "avg_distr")){

  # check arguments
  if (missing(data))
    stop("'data' argument must be provided")
  if (missing(df))
    stop("'df' argument must be provided")
  if (missing(type)){
    type <- type[1]
  }
  if (missing(what)){
    what <- what[1]
  }

  # check the type of argument
  if(!(
    is(data, "SummarizedExperiment") | is.data.frame(data) | is.matrix(data))
  )
    stop("'data' must be a 'matrix', a 'data.frame'
         or a 'SummarizedExperiment' object")
  if(!(is(df, "DataFrame") | is.data.frame(df)))
    stop("'df' must be a data.frame")

  if (is(data, "SummarizedExperiment")) {
    count_data<-assay(data)
  } else {
    count_data<-as.matrix(data)
  }
  df<-as.data.frame(df)

  if (!(all(what %in% c("all",
                        "all_w_PCA",
                        "MDS",
                        "PCA",
                        "heatmap",
                        "RLEbox",
                        "distr",
                        "avg_distr")))
  )
    stop("'what' must be one of
         all, 'all_w_PCA','MDS', 'PCA',
         'heatmap', 'RLEbox', 'distr', 'avg_distr'")



  # check the presence of NA or Inf
  if (any(is.na(count_data)))
    stop("NA values are not allowed in the 'data' matrix")
  # if (any(is.na(df)))
  #   stop("NA values are not allowed in the 'df' matrix")
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
  if (length(type) > 1)
    stop("length(type) must be equal to 1")
  if (!(all(type %in% c("pearson", "spearman"))))
    stop("'type' must be 'pearson' or 'spearman'")



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

  if (what == "all" | what == "all_w_PCA" | what == "heatmap"){
    if (ncol(df)< 30){
      pheatmap(sampleDistMatrix,
               clustering_distance_rows=mydist,
               clustering_distance_cols=mydist,
               #col=colors,
               #breaks = seque,
               breaks = seq(0,1,by=0.01),
               annotation_col = df,
               main="Heatmap by Dissimilarity")

    }else{
      warning("Too many 'df' variables provided for plotting heatmab by
    pheatmap. Please, split 'df' in more subsets")
    }
  }
  #options(warn = -1)
  ################
  ## MDS plot
  mdsData <- data.frame(cmdscale(sampleDistMatrix))

  if (what == "all" | what == "MDS"){
    for(i in seq_len(ncol(df))) {
      mds <- cbind(mdsData, df)
      cov_list <- mds[,i+2,drop=FALSE]
      colnames(cov_list)<-"Vars"
      # cov_list$Vars <- droplevels(cov_list$Vars)
      #checkNA
      NA_idx <- which(is.na(cov_list$Vars))

      if(length(NA_idx) != 0){
        cov_list <- cov_list[-NA_idx,,drop=FALSE]
        mds <- mds[-NA_idx,]
        #cov_list$Vars <- droplevels(cov_list$Vars)
        cat("MDS:",length(NA_idx),
            "samples not drawn (because 'NA') for variable:",
            colnames(df)[i],"\n")
      }

      # print(ggplot(mds, aes(X1, X2, shape=mds$class, color=cov_list$Vars)) +
      print(ggplot(mds, aes(X1, X2, shape=class, color=cov_list$Vars)) +
              geom_point(size=3) +
              geom_text(aes(label=rownames(mds)),hjust=0.5,vjust=-1) +
              ggtitle(paste("Variable: ",colnames(mds[,i+2,drop=FALSE]))))
    }
  }

  #### PCA
  if (what == "all_w_PCA" | what == "PCA"){

    PC_var <-prcomp(t(count_data),center = T,scale. = T)
    df_PC <- as.data.frame(PC_var$x)
    df_PC <- df_PC[,1:3]
    colnames(df_PC) <- c("PC1","PC2","PC3")


    for(i in seq_len(ncol(df))) {
      df_PC_ok <- cbind(df_PC, df)
      cov_list <- df_PC_ok[,i+3,drop=FALSE]
      colnames(cov_list)<-"Vars"
      # cov_list$Vars <- droplevels(cov_list$Vars)
      #checkNA
      NA_idx <- which(is.na(cov_list$Vars))

      if(length(NA_idx) != 0){
        cov_list <- cov_list[-NA_idx,,drop=FALSE]
        df_PC_ok <- df_PC_ok[-NA_idx,]
        #cov_list$Vars <- droplevels(cov_list$Vars)
        cat("MDS:",length(NA_idx),
            "samples not drawn (because 'NA') for variable:",
            colnames(df)[i],"\n")
      }

      print(ggplot(df_PC_ok, aes(PC1, PC2, shape=class, color=cov_list$Vars)) +
              geom_point(size=3) +
              geom_text(aes(label=rownames(df_PC_ok)),hjust=0.5,vjust=-1) +
              ggtitle(paste("Variable: ",colnames(df_PC_ok[,i+3,drop=FALSE]))))
    }


  }

  ################
  ## RLE
  colors <- brewer.pal(8, "Set2")

  if (what == "all" | what == "all_w_PCA" | what == "RLEbox"){
    plotRLE(count_data,
            k=2,
            labels=TRUE,
            isLog=TRUE,
            outline=FALSE,
            col=colors[as.factor(df$class)],
            main="Relative Log Expression",
            xaxt="n",las=1)
    axis(1,
         las=2,
         at=seq_len(dim(count_data)[2]),
         labels=colnames(count_data))
  }
  ####################
  ## Sample by sample distribution
  if (what == "all" | what == "all_w_PCA" | what == "distr"){
    acc_dotplot <- melt(as.data.frame(count_data),
                        measure.vars = colnames(count_data))
    print(ggplot(acc_dotplot, aes(value,
                                  fill = variable,
                                  colours= variable)) +
            geom_density(alpha=0.3) +
            facet_wrap(~variable)+
            theme(legend.position = "none") +
            ggtitle("Sample by Sample expression value distribution")
    )
  }

  ####################
  ## Class average expression distribution
  if (what == "all" | what == "all_w_PCA" | what == "avg_distr"){
    dataset_2 <- t(count_data)
    livelli <- levels(as.factor(df$class))

    for (i in seq_len(length(livelli))){
      dataset_tmp <- dataset_2[which(df$class %in% livelli[i] ),,drop=FALSE]
      colmeans_tmp <- colMeans(dataset_tmp)
      if (i == 1){
        dataset_3 <- as.data.frame(colmeans_tmp)
      }else{
        dataset_3 <- cbind(dataset_3,colmeans_tmp)
      }
    }

    colnames(dataset_3) <-livelli
    acc_dotplot_2 <- melt(dataset_3,
                          measure.vars = colnames(dataset_3))
    print(ggplot(acc_dotplot_2, aes(value,
                                    fill = variable,
                                    colours= variable)) +
            geom_density(alpha=0.3) +
            ggtitle("Class average expression value distribution")
    )
  }




}


#' @title Plot multidimentional scaling (MDS)
#'
#' @description A MDS plot is drawn in order to visualize class clustering.
#'
#' @param data A SummarizedExperiment object or a matrix or a data.frame
#'  where
#' rows and cols should be, respectively, observations and features
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
  if(!(
    is(data, "SummarizedExperiment") | is.data.frame(data) | is.matrix(data))
  )
    stop("'data' must be a 'matrix', a 'data.frame'
         or a 'SummarizedExperiment' object")
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
  if (length(type) > 1)
    stop("length(type) must be equal to 1")
  if (!(all(type %in% c("pearson", "spearman"))))
    stop("'type' must be 'pearson' or 'spearman'")


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
#  ggplot(mds, aes(X1, X2, shape=df$class, color=df$class)) +
   ggplot(mds, aes(X1, X2, shape=class, color=class)) +
    geom_point(size=3) +
    geom_text(aes(label=rownames(mds)), hjust=0.5, vjust=-1)

}

#' @title Expression data clustering and heatmap
#'
#' @description The function helps to draw a clustering dendrogram and
#'  a heatmap of expression data.
#'
#' @param data A SummarizedExperiment object or a matrix or a data.frame
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
  if(!(
    is(data, "SummarizedExperiment") | is.data.frame(data) | is.matrix(data))
  )
    stop("'data' must be a 'matrix', a 'data.frame'
         or a 'SummarizedExperiment' object")
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




