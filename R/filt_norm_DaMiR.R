#' @title Filter non Expressed and 'Hypervariant' features and Data Normalization
#'
#' @description Features will be firstly filtered based on their expression value and/or by
#'  their variability across samples; features will be then normalized.
#'
#' @param dds A \code{DESeq} object
#' @param minCounts Minimum reads counts; default is 10
#' @param fSample Fraction of samples with \code{minCounts} counts; default is 0.5
#' @param hyper Flag to enable gene filtering by Coefficient of
#' Variation (CV); default is "yes"
#' @param th.cv Threshold of minimum CV to consider a feature 'Hypervariant'
#'  accross samples; default is 3
#' @param type Type of normalization to be applied: \code{varianceStabilizingTransformation}
#' (\code{vst}) or \code{rlog} are allowed; default is "\code{vst}"
#'
#' @details
#' Before normalization step, this function allows the user to filter features by:
#' \itemize{
#'   \item Expression - Features will be filtered out whether their reads count
#'   do not reach a \code{minCounts} in at least \code{fSample} of samples;
#'   \item CV - The CV of each feature is individually calculated for each sample class.
#'   Featurers with both class CV greater than \code{th.cv} will be discarded.
#'   Computing a class restricted CV may prevent the removal of hypervariant features that
#'   may be specifically associated with a certain class. This could be important, for example, for
#'   immune genes whose expression under definite conditions may unveil peculiar class-gene
#'   association.
#' }
#' Finally, expressed features will be normalized by \code{varianceStabilizingTransformation}
#'  (default) or \code{rlog}, both implemented in \code{DESeq2} package. We suggest to
#'  use \code{varianceStabilizingTransformation} to speed up the normalization process
#'  because \code{rlog} is very time-consuming despite the two methods produce quite
#'  similar results.
#'
#' @return A normalized expression matrix.
#'
#' @references Michael I Love, Wolfgang Huber and Simon Anders (2014): Moderated estimation of
#'  fold change and dispersion for RNA-Seq data with DESeq2. Genome Biology
#'
#' @author Mattia Chiesa, Luca Piacentini
#'
#' @seealso
#' \code{\link{varianceStabilizingTransformation}, \link{rlog} }
#'
#' @examples
#' rawdata.path <- system.file(package = "DaMiRseq","extdata")
#' filecounts <- list.files(rawdata.path,full.names = TRUE)[1]
#' filecovariates <- list.files(rawdata.path,full.names = TRUE)[2]
#' list.res <- DaMiR.read(filecounts,filecovariates)
#' dds <- list.res$dds
#'
#' data_norm <- DaMiR.normalization(dds)
#' data_norm <- DaMiR.normalization(dds, fSample = 0.6)
#' data_norm <- DaMiR.normalization(dds, th.cv = 2.5, type = "rlog")
#'
#' @export
#'
#'
DaMiR.normalization <- function(dds, minCounts=10, fSample=0.5, hyper=c("yes","no"), th.cv=3, type=c("vst","rlog")){
  # check arguments
  if (missing(dds)) stop("'dds' argument must be provided")
  if (missing(hyper)){
    hyper <- hyper[1]
  }
  if (missing(type)){
    type <- type[1]
  }

  if(!(is(dds, "DESeqDataSet"))) stop("'dds' must be a 'DESeqDataSet' object")
  if(!(is.numeric(minCounts))) stop("'minCounts' must be numeric")
  if(!(is.numeric(fSample))) stop("'fSample' must be numeric")
  if(!(is.numeric(th.cv))) stop("'th.cv' must be numeric")

  # start execution
  init_lenght<-dim(dds)[1]

  #filtering by Expression Level
  minSamplesExpr<-round((dim(dds)[2])*fSample)
  dds <- dds[rowSums(counts(dds) > minCounts) > minSamplesExpr,]
  dds_exp_count <- counts(dds)
  cat(init_lenght-dim(dds)[1],"Features have been filtered out by espression.",dim(dds)[1], "Features remained.","\n")

  #filtering by CV
  if (hyper == "yes"){
    init_lenght_cv<-dim(dds)[1]

    label<-dds@colData$class
    patt_cl1<-paste(c("^",levels(dds@colData$class)[1],"$"),sep="",collapse="")
    patt_cl2<-paste(c("^",levels(dds@colData$class)[2],"$"),sep="",collapse="")

    cv_value1 <- 0
    cv_value2 <- 0
    funct_cv <- function(matr_4_cv){
      cv_val<- sd(matr_4_cv)/abs(mean(matr_4_cv))
    }
    cv_value1 <- apply(dds_exp_count[,grep(patt_cl1,dds@colData$class)],1,funct_cv)
    cv_value2 <- apply(dds_exp_count[,grep(patt_cl2,dds@colData$class)],1,funct_cv)

    dds <- dds[cv_value1 < th.cv & cv_value2 < th.cv ,]
    dds_exp_count <- counts(dds)
    cat(init_lenght_cv-dim(dds)[1],"'Hypervariant' Features have been filtered out.",dim(dds)[1], "Features remained.","\n")

  }

  # Normalizion by VST
  if(type == "vst"){
    cat("Performing Normalization by 'vst'","\n")
    data_norm <-varianceStabilizingTransformation(dds_exp_count)
  } else if (type == "rlog"){
    cat("Performing Normalization by 'rlog'. For large dataset it could be very time-consuming.","\n")
    data_norm <-rlog(dds_exp_count)}
  else{
    stop("Please set 'vst or 'rlog' as normalization type.")
  }

  return(data_norm)
}

#' @title Filter Samples by Mean Correlation Distance Metric
#'
#' @description This function implements a sample-per-sample correlation. Samples with a mean
#' correlation lower than a user's defined threshold will be filtered out.
#'
#' @param data Matrix of normalized expression data, i.e transformed counts by vst or rlog.
#' A log2 transformed expression matrix is also accepted
#' @param df A data frame with known covariates; at least one column with 'class'
#' label must be included
#' @param th.corr Threshold of mean correlation; default is 0.9
#' @param type Type of correlation metric; default is "spearman"
#'
#' @details This step introduces a sample quality checkpoint. Global gene expression should,
#' in fact, exhibit a high correlation among biological replicates; conversely, low correlated
#' samples may be suspected to bear some technical artifact (e.g. poor RNA or library
#' preparation quality), despite they may have passed sequencing quality checks. If not assessed,
#' these samples may, thus, negatively affect all the downstream analysis. This function looks at
#' the mean absolute correlation of each sample and removes those samples with a mean correlation
#' lower than the value set in \code{th.corr} argument. This threshold may be specific for
#' different experimental setting but should be as high as possible. For sequencing data we
#' suggest to set \code{th.corr} greater than 0.85.
#'
#' @return A list containing:
#' \itemize{
#'   \item A normalized and filtered expression matrix.
#'   \item A data frame with class and optional covariates information.
#' }
#'
#' @author Mattia Chiesa, Luca Piacentini
#'
#' @examples
#' data(DaMiR.Data.Ex.Norm)
#' data(DaMiR.Data.Ex.dfClCov)
#' list.res <- DaMiR.sampleFilt(DaMiR.Data.Ex.Norm, DaMiR.Data.Ex.dfClCov)
#' list.res <- DaMiR.sampleFilt(DaMiR.Data.Ex.Norm, DaMiR.Data.Ex.dfClCov, th.corr = 0.8)
#' list.res <- DaMiR.sampleFilt(DaMiR.Data.Ex.Norm, DaMiR.Data.Ex.dfClCov, type = "pearson")
#'
#' @export
#'
#'
DaMiR.sampleFilt <- function(data, df, th.corr=0.9, type=c("spearman","pearson")){

  # check arguments
  if (missing(data)) stop("'data' argument must be provided")
  if (missing(df)) stop("'df' argument must be provided")
  if (missing(type)){
    type <- type[1]
  }

  if(!(is.matrix(data))) stop("'data' must be a matrix")
  if(!(is.data.frame(df))) stop("'df' must be a data frame")
  if(!(is.numeric(th.corr))) stop("'th.corr' must be numeric")

  init_length_corr<-dim(data)[2]
  # Sample-by-Sample correlation
  if(type == "spearman"){
    cormatrix <- abs(rcorr(as.matrix(data), type='spearman')$r)}
  else if(type == "pearson"){
    cormatrix <- abs(rcorr(as.matrix(data), type='pearson')$r)
  } else{
    stop("Please set 'spearman or 'pearson' as correlation type.")
  }

  mean_corr <- rowMeans(cormatrix)

  #check
  if(length(which(mean_corr>=th.corr)) == 0){
    stop("You are removing all samples. Please decrease 'th.corr' value.")
  }else{

  # update matrices
  data <- data[,mean_corr>=th.corr]
  df <-df[mean_corr>=th.corr,]
  cat(init_length_corr-dim(data)[2],"Samples have been excluded by averaged Sample-per-Sample correlation.","\n",dim(data)[2], "Samples remained.","\n")

  return(list(data = data, classCov = df))
  }
}


