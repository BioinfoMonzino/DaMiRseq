#' @title Feature selection for classification
#'
#' @description This function identifies the class-correlated principal
#'  components (PCs)
#' which are then used to implement a backward variable elimination
#' procedure for the removal of non informative features.
#'
#'
#' @param data A transposed data frame or a matrix of normalized expression
#'  data.
#' Rows and Cols should be,
#' respectively, observations and features
#' @param df A data frame with known variables; at least one column
#' with
#' 'class' label must be included
#' @param th.corr Minimum threshold of correlation between class and
#' PCs; default is 0.6
#' @param type Type of correlation metric; default is "spearman"
#' @param th.VIP Threshold for \code{bve_pls} function, to remove
#' non-important variables; default is 3
#' @param nPlsIter Number of times that \link{bve_pls} has to run.
#' Each iteration produces a set of selected features, usually similar
#' to each other but not exacly the same! When nPlsIter is > 1, the
#' intersection between each set of selected features is performed;
#' so that, only the most robust features are selected. Default is 1
#'
#' @details The function aims to reduce the number of features to obtain
#' the most informative
#' variables for classification purpose. First, PCs obtained by principal
#'  component analysis (PCA)
#' are correlated with "class". The correlation is defined by the user in
#'  \code{th.corr}
#' argument. The higher is the correlation, the lower is the number of PCs
#'  returned.
#' Users should pay attention to appropriately set the \code{th.corr}
#' argument because
#' it will also affect the total number of selected features that ultimately
#'  depend on the number of PCs.
#' The \code{\link{bve_pls}} function of \code{plsVarSel} package is, then,
#'  applied.
#' This function exploits a backward variable elimination procedure coupled
#' to a partial least squares approach to remove those variable which are
#' less informative with
#' respect to class. The returned vector of variables is further reduced by
#'  the following
#' \code{\link{DaMiR.FReduct}} function in order to obtain a subset of
#' non correlated
#' putative predictors.
#'
#'
#' @return A list containing:
#' \itemize{
#'   \item An expression matrix with only informative features.
#'   \item A data frame with class and optional variables information.
#' }
#'
#' @references Tahir Mehmood, Kristian Hovde Liland, Lars Snipen and
#' Solve Saebo (2011).
#' A review of variable selection methods in Partial Least Squares
#' Regression. Chemometrics and
#' Intelligent Laboratory Systems 118, pp. 62-69.
#'
#' @author Mattia Chiesa, Luca Piacentini
#'
#' @seealso
#'  \itemize{
#'   \item \code{\link{bve_pls}}
#'   \item \code{\link{DaMiR.FReduct}}
#' }
#'
#'
#' @examples
#' # use example data:
#' data(data_norm)
#' data(df)
#' # extract expression data from SummarizedExperiment object
#' # and transpose the matrix:
#' t_data<-t(assay(data_norm))
#' t_data <- t_data[,seq_len(100)]
#' # select class-related features
#' data_reduced <- DaMiR.FSelect(t_data, df,
#' th.corr = 0.7, type = "spearman", th.VIP = 1)
#'
#' @export
#'
DaMiR.FSelect <- function(data,
                          df,
                          th.corr=0.6,
                          type=c("spearman", "pearson"),
                          th.VIP=3,
                          nPlsIter=1){
  # check missing arguments
  if (missing(data))
    stop("'data' argument must be provided")
  if (missing(df))
    stop("'df' argument must be provided")
  if (missing(type)){
    type <- type[1]
  }

  # check the type of argument
  if(!(is.matrix(data) | is.data.frame(data)))
    stop("'data' must be a matrix or a data.frame")
  if(!(is(df, "DataFrame") | is.data.frame(df)))
    stop("'df' must be a data.frame")
  if(!(is.numeric(th.corr)))
    stop("'th.corr' must be numeric")
  if(!(is.numeric(th.VIP)))
    stop("'th.VIP' must be numeric")
  if(!(is.numeric(nPlsIter)))
    stop("'nPlsIter' must be numeric")


  data<-as.data.frame(data)
  df<-as.data.frame(df)

  # check the presence of NA or Inf
  if (any(is.na(data)))
    stop("NA values are not allowed in the 'data' matrix")
  if (any(is.na(df)))
    stop("NA values are not allowed in the 'df' matrix")
  if (any(is.infinite(as.matrix(data))))
    stop("Inf values are not allowed in the 'data' matrix")

  # specific checks
  if (th.corr >1 | th.corr < 0)
    stop("'th.corr' must be between 0 and 1")
  if (th.VIP <= 0)
    stop("'th.VIP' must be positive")
  if (all((as.matrix(data) %%1) == 0))
    warning("It seems that you are using raw counts!
            This function works with normalized data")
  if(dim(as.matrix(data))[1] != dim(df)[1])
    stop("ncol(as.matrix(data)) must be equal to nrow(sv)")
  if(!("class" %in% colnames(df)))
    stop("'class' info is lacking!
         Include the variable 'class'
         in the 'df' data frame and label it 'class'!")
  if (nPlsIter <= 0)
    stop("'nPlsIter' must be positive")
  if ((nPlsIter %%1) != 0)
    stop("nPlsIter must be integer")
  if (length(type) > 1)
    stop("length(type) must be equal to 1")
  if (!(all(type %in% c("pearson", "spearman"))))
    stop("'type' must be 'pearson' or 'spearman'")

  correlation <-0
  features<-dim(data)[2]

  # identify PCs
  pca_plot<-FactoMineR::PCA(data,
                            graph = FALSE,
                            ncp = round(nrow(data)/2))

  # class vs PCs correlation
  if(type == "spearman"){
    class.vs.PCs<-function(x){
      correlation <- abs(cor(x, (as.numeric(df$class)-1),
                             method = "spearman"))
    }
  } else if(type == "pearson"){
    class.vs.PCs<-function(x){
      correlation <- abs(cor(x,(as.numeric(df$class)-1)))
    }
  } else{
    stop("Please set 'spearman or 'pearson' as correlation type.")
  }
  correlation<-apply(pca_plot$ind$coord, 2, class.vs.PCs)

  # identify informative PCs
  if ( length(which(correlation>=th.corr)) == 0 ){
    m_stop <- paste("There are not PCs with a correlation grater than",
                    th.corr,
                    "with class. Please decrease 'th.corr'")
    stop(m_stop)
  }else{
  num_PC <- max(which(correlation>=th.corr))
  }

  # filter out uninformative features
  if (nPlsIter == 1){
    featureslist<-bve_pls(df$class,
                               as.matrix(data),
                               ncomp=num_PC,
                               VIP.threshold=th.VIP)
    selectedCandidates <- names(data)[unlist(featureslist)]
  }else{
    for (i in seq_len(nPlsIter)){
      if (i== 1){
        featureslist<-bve_pls(df$class,
                              as.matrix(data),
                              ncomp=num_PC,
                              VIP.threshold=th.VIP)
        selectedCandidates <- names(data)[unlist(featureslist)]
      }else{
       featureslist<-bve_pls(df$class,
                      as.matrix(data),
                      ncomp=num_PC,
                      VIP.threshold=th.VIP)
       selectedCandidates_iter <- names(data)[unlist(featureslist)]
       selectedCandidates <- intersect(selectedCandidates,
                                       selectedCandidates_iter)
      }
    }
  }
  data_reduced<-data[, which(names(data) %in% selectedCandidates)]
  if (length(colnames(data_reduced))==0) {
    cat("All the genes have been discarded!!!
        th.VIP argument is too high. Please choose a lower level of th.VIP")
  } else {

  cat(features-dim(data_reduced)[2],
      "Genes have been discarded for classification",
      dim(data_reduced)[2],
      "Genes remained.","\n")

  return(list(data = data_reduced, variables = df))

  }
}

#' @title Remove highly correlated features, based on feature-per-feature
#'  correlation.
#'
#' @description This function allows the user to remove highly correlated
#'  features.
#'
#' @param data A transposed data frame or matrix of normalized expression
#' data. Rows and Cols should be,
#' respectively, observations and features
#' @param th.corr Feature-per-feature correlation threshold; default is
#'  0.85
#' @param type Type of correlation metric to be applied; default is
#' "spearman"
#'
#' @details This function produces an absolute correlation matrix that
#' it is then used to reduce pair-wise correlations.
#' When two features present a correlation higher than that defined by
#' the user in \code{th.corr} argument,
#' the function, first, calculates the mean absolute correlation of each
#'  feature and, then, removes the feature
#' with the largest mean absolute correlation.
#'
#'
#' @return An expression matrix without highly correlated features.
#'
#' @author Mattia Chiesa, Luca Piacentini
#'
#' @seealso
#' \code{\link{rcorr}}, \code{\link{findCorrelation}}
#'
#' @examples
#' # use example data:
#' data(data_reduced)
#' # reduce the number of features:
#' data_Reduced <- DaMiR.FReduct(data_reduced,
#' th.corr = 0.75, type = "pearson")
#'
#' @export
#'
DaMiR.FReduct <- function(data,
                          th.corr=0.85,
                          type=c("spearman","pearson")){

  # check arguments
  if (missing(data)) stop("'data' argument must be provided")
  if (missing(type)){
    type <- type[1]
  }

  # check the type of argument
  if(!(is.numeric(th.corr)))
    stop("'th.corr' must be numeric")
  if(!(is.data.frame(data)))
    stop("'data' must be a data.frame")

  # check the presence of NA or Inf
  if (any(is.na(data)))
    stop("NA values are not allowed in the 'data' matrix")
  if (any(is.infinite(as.matrix(data))))
    stop("Inf values are not allowed in the 'data' matrix")

  # specific checks
  if (th.corr >1 | th.corr < 0)
    stop("'th.corr must be between 0 and 1")
  if (all((as.matrix(data) %%1) == 0))
    warning("It seems that you are using raw counts!
            This function works with normalized data")

  features<-dim(data)[2]

  # remove redundancy
  if(type == "spearman"){
  cormatrix <- abs(rcorr(as.matrix(data), type='spearman')$r)
  } else if (type == "pearson"){
    cormatrix <- abs(rcorr(as.matrix(data), type='pearson')$r)
  } else {
    stop("Please set 'spearman or 'pearson' as correlation type.")
  }

  index_geneHighCorr<- findCorrelation(cormatrix, cutoff = th.corr)
  data_reduced<-data[, -index_geneHighCorr, drop=FALSE]

  cat(features-dim(data_reduced)[2],
      "Highly correlated features have been discarded for classification.",
      "\n",
      dim(data_reduced)[2],
      "Features remained.",
      "\n")
  return(data_reduced)
}

#' @title Order features by importance, using RReliefF filter
#'
#' @description This function implements a procedure in order to rank
#'  features by their
#' importance evaluated by RReliefF score.
#'
#' @param data A transposed data frame of expression data, i.e.
#' transformed
#' counts by vst or rlog. A log2 transformed expression matrix is also
#'  accepted.
#' Rows and Cols should be, respectively, observations and features
#' @param df A data frame with class and known variables; at least one
#' column with
#' 'class' label must be included
#' @param fSample Fraction of sample to be used for the implementation
#' of RReliefF algorithm; default is 1
#'
#' @return A data frame with two culmuns, where features are sorted by
#' importance scores:
#' \itemize{
#'   \item RReliefF score - Calculated by \code{\link{relief}} function,
#'    implemented in
#'   \code{FSelector} package;
#'   \item scaled.RReliefF  score - Z-score value, computed for each
#'   RReliefF score.
#' }
#' A plot with the first 50 features ordered by their importance.
#'
#' @details This function is very time-consuming when the number of
#' features is
#' high. We observed there is a quadratic relationship between execution
#'  time
#' and the number of features. Thus, we have also implemented a formula
#'  which allows
#' the users to estimate the time to perform this step, given the number
#'  of features.
#' The formula is:
#' \deqn{T = 0.0011 * N^2 - 0.1822 * N + 27.092}
#' where T = Time and N = Number of genes.
#' We strongly suggest to filter out non informative features before
#' performing
#' this step.
#'
#' @author Mattia Chiesa, Luca Piacentini
#'
#' @seealso
#' \code{\link{relief}}, \code{\link{DaMiR.FSelect}},
#'  \code{\link{DaMiR.FReduct}}
#'
#' @references
#' Marko Robnik-Sikonja, Igor Kononenko: An adaptation of Relief for
#' attribute
#' estimation in regression. In: Fourteenth International Conference
#' on Machine
#' Learning, 296-304, 1997
#'
#' @examples
#' # use example data:
#' data(data_reduced)
#' data(df)
#' # rank features by importance:
#' df.importance <- DaMiR.FSort(data_reduced[,1:10],
#'  df, fSample = 0.75)
#'
#' @export
#'
#'
DaMiR.FSort <- function(data, df, fSample=1){

  # check missing arguments
  if (missing(data))
    stop("'data' argument must be provided")
  if (missing(df))
    stop("'df' argument must be provided")

  # check the type of argument
  if(!(is.data.frame(data)))
    stop("'data' must be a data.frame")
  if(!(is(df, "DataFrame") | is.data.frame(df)))
    stop("'df' must be a data.frame")
  if(!(is.numeric(fSample)))
    stop("'fSample' must be a numeric")

  df<-as.data.frame(df)

  # check the presence of NA or Inf
  if (any(is.na(data)))
    stop("NA values are not allowed in the 'data' matrix")
  if (any(is.na(df)))
    stop("NA values are not allowed in the 'df' matrix")
  if (any(is.infinite(as.matrix(data))))
    stop("Inf values are not allowed in the 'data' matrix")

  # specific checks
  if (fSample >1 | fSample <= 0)
    stop("'fSample must be between 0 and 1")
  if (all((as.matrix(data) %%1) == 0))
    warning("It seems that you are using raw counts!
            This function works with normalized data")
  if(dim(as.matrix(data))[1] != dim(df)[1])
    stop("nrow(as.matrix(data)) must be equal to nrow(df)")
  if(!("class" %in% colnames(df)))
    stop("'class' info is lacking!
         Include the variable 'class'
         in the 'df' data frame and label it 'class'!")



  # estimated time for reliefF
  x <- dim(data)[2]
  y <-0.0011*x^2 - 0.1822*x + 27.092

  cat("Please wait. This operation will take about", round(y),
      "seconds (i.e. about",round(y/60),"minutes).")

  classes <- levels(df$class)

  # find the min number of sample per class
  min_sub <-0
  for (i in seq_len(length(classes))){
    min_sub[i] <- length(which(levels(df$class)[i]==df$class))
  }
    n_sample_relief <- round(min(min_sub)*fSample)
  if (n_sample_relief == 0 )
    stop("No subject has been selected for RRelieF.
         Please increase 'fSample'.")

  # find the indexes for sampling
  sample_index_list <- ""
  for (i in seq_len(length(classes))){
    sample_index <- sample(which(levels(df$class)[i]==df$class),replace = FALSE)
    sample_index <- sample_index[seq_len(n_sample_relief)]
    sample_index_list <- c(sample_index_list,sample_index)
  }
  sample_index_list <- as.numeric(sample_index_list[-1])
  
  # bug fixed. Thanks to Dr. Priti Prasad!
  if (length(sample_index_list) < 3 )
    stop("Few subjects have been selected for RRelieF.")

  # subset data for relieF
  dataset.relief <- data[sample_index_list,,drop=FALSE]
  classes.relief <- df$class[sample_index_list]
  dataset.relief$classes <- classes.relief

  # reliefF
  formula.relief <- as.formula(paste("classes",
                                     paste(names(data),
                                           collapse = "+"),
                                     sep = " ~ "))
  if (length(sample_index_list) < 10 ){
    rel_res <- relief(formula = formula.relief, data = dataset.relief,
	                  sample.size = length(sample_index_list),
					  neighbours.count = floor(length(sample_index_list)/2))
  }else{
    rel_res <- relief(formula = formula.relief, data = dataset.relief)
  }
  

  # Importance plot
  imp_attrib <- rel_res[as.numeric(order(rel_res$attr_importance,
                                         decreasing = TRUE)),, drop=FALSE]
  imp_attrib$attr_importance_scaled <-scale(as.matrix(
    imp_attrib$attr_importance))
  colnames(imp_attrib) <- c("RReliefF","scaled.RReliefF")

  if (dim(imp_attrib)[1] < 50){
    top_50_imp <- imp_attrib[seq_len(dim(imp_attrib)[1]),1, drop=FALSE]
    top_50_imp <- top_50_imp[rev(rownames(top_50_imp)),, drop=FALSE]
    colnames(top_50_imp) <- "Top features"
    dotchart(as.matrix(top_50_imp), xlab = "RReliefF importance",
             main = "Attributes importance by RReliefF")
  } else{
    top_50_imp <- imp_attrib[seq_len(50),1, drop=FALSE]
    top_50_imp <- top_50_imp[rev(rownames(top_50_imp)),, drop=FALSE]
    colnames(top_50_imp) <- "Top50 features"
    dotchart(as.matrix(top_50_imp), xlab = "RReliefF importance",
           main = "Attributes importance by RReliefF")
  }
  return(imp_attrib)
}

#' @title Select best predictors to build Classification Model
#'
#' @description This function allows the user to select a subset
#' of predictors; the number
#' of predictors can be defined by user or selected automatically.
#'
#' @param data A transposed data frame of expression data.
#' Rows and Cols should be, respectively, observations and features
#' @param ranking A data frame with importance score for each feature,
#' generated by
#'  \code{\link{DaMiR.FSort}}
#' @param autoselect A flag to specify how to select predictors:
#' \itemize{
#'     \item{"no" (default) - Manually: users can specify the number
#'     of best predictors,
#'     setting \code{n.pred} argument}
#'     \item{"yes" - Automatically: users have to specify the importance
#'      threshold
#'     defined by the \code{th.zscore} argument; features will be
#'     accordingly selected}
#' }
#' @param n.pred If \code{autoselect="no"} then the user have to specify
#'  the number of
#' predictors; default is 10
#' @param th.zscore Threshold of scaled importance score (Z-score);
#' default value is 2
#'
#' @return A list containing:
#' \itemize{
#'   \item A data frame of normalized expression data of the most important
#'   selected predictors.
#'   \item A vector with predictors name.
#' }
#'
#' @author Mattia Chiesa, Luca Piacentini
#'
#' @examples
#' # use example data:
#' data(data_reduced)
#' data(data_relief)
#' # select the first 8 predictors rankad by imporatance:
#' selected_features <- DaMiR.FBest(data_reduced, data_relief, n.pred = 8)
#' # select predictors by importance but automatically:
#' selected_features <- DaMiR.FBest(data_reduced, data_relief,
#' autoselect = "yes", th.zscore = 1.5)
#'
#' @seealso
#' \code{\link{DaMiR.FSort}}
#'
#' @export
#'
#'
DaMiR.FBest <- function(data,
                        ranking,
                        autoselect=c("no","yes"),
                        n.pred=10,
                        th.zscore=2){

  # check missing arguments
  if (missing(data))
    stop("'data' argument must be provided")
  if (missing(ranking))
    stop("'ranking' argument must be provided")
  if (missing(autoselect)){
    autoselect <- autoselect[1]
  }

  # check the type of argument
  if(!(is.data.frame(data)))
    stop("'data' must be a data frame")
  if(!(is.data.frame(ranking)))
    stop("'ranking' must be a data frame")
  if(!(is.numeric(n.pred)))
    stop("'n.pred' must be numeric")
  if(!(is.numeric(th.zscore)))
    stop("'th.zscore' must be numeric")


  # check the presence of NA or Inf
  if (any(is.na(data)))
    stop("NA values are not allowed in the 'data' matrix")
  if (any(is.na(ranking)))
    stop("NA values are not allowed in the 'ranking' matrix")
  if (any(is.infinite(as.matrix(data))))
    stop("Inf values are not allowed in the 'data' matrix")

  # specific checks
  if (all((as.matrix(data) %%1) == 0))
    warning("It seems that you are using raw counts!
            This function works with normalized data")
  if(dim(as.matrix(data))[2] < dim(ranking)[1])
    stop("ncol(as.matrix(data)) have not to be lower than nrow(ranking)")
  if(dim(ranking)[2] != 2)
    stop("'ranking' must have 2 columns")



  if (autoselect == "no"){
    if(missing(th.zscore)){
    if(n.pred < 1)
      stop("'n.pred' must be greater than 0")
    if(n.pred > dim(ranking)[1])
      stop("'n.pred' must be lower than dim(ranking)[1]" )
    predictors <- rownames(ranking[seq_len(n.pred),,drop=FALSE])}
    else{
      stop("'th.zscore' must be set only with 'autoselect = yes'")
    }

  } else if (autoselect == "yes"){
    if(missing(n.pred)){
    predictors <- rownames(
    ranking[ranking[,2]>th.zscore,, drop=FALSE])
    } else {
      stop("'n.pred' must be set only with 'autoselect = no'")
    }
  } else {
    stop("Please set 'yes or 'no' in 'autoselect' option.")
    }

  data_for_class <- data[,colnames(data) %in% predictors]
  cat(length(predictors),
      "Predictors have been selected for classification",
      "\n")
  #cat("Predictors:","\n", predictors)

  return(list(data = data_for_class, predictors = predictors))

}



