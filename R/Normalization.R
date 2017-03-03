#' @title Filter non Expressed and 'Hypervariant' features and Data
#' Normalization
#'
#' @description Features will be firstly filtered based on their
#' expression value and/or by
#'  their variability across samples; features will be then normalized.
#'
#' @param data A \code{SummarizedExperiment} object
#' @param minCounts Minimum reads counts; default is 10
#' @param fSample Fraction of samples with \code{minCounts} counts; default
#'  is 0.5
#' @param hyper Flag to enable gene filtering by Coefficient of
#' Variation (CV); default is "yes"
#' @param th.cv Threshold of minimum CV to consider a feature 'Hypervariant'
#'  accross samples; default is 3
#' @param type Type of normalization to be applied:
#' \code{varianceStabilizingTransformation}
#' (\code{vst}) or \code{rlog} are allowed; default is "\code{vst}"
#'
#' @details
#' Before normalization step, this function allows the user to filter
#' features by:
#' \itemize{
#'   \item Expression - Features will be filtered out whether their reads
#'    count
#'   do not reach a \code{minCounts} in at least \code{fSample} of samples;
#'   \item CV - The CV of each feature is individually calculated for each
#'    sample class.
#'   Featurers with both class CV greater than \code{th.cv} will be
#'   discarded.
#'   Computing a class restricted CV may prevent the removal of hypervariant
#'   features that
#'   may be specifically associated with a certain class. This could be
#'   important, for example, for
#'   immune genes whose expression under definite conditions may unveil
#'   peculiar class-gene
#'   association.
#' }
#' Finally, expressed features will be normalized by
#' \code{varianceStabilizingTransformation}
#'  (default) or \code{rlog}, both implemented in \code{DESeq2} package.
#'  We suggest to
#'  use \code{varianceStabilizingTransformation} to speed up the
#'  normalization process
#'  because \code{rlog} is very time-consuming despite the two methods
#'   produce quite
#'  similar results.
#'
#' @return A \code{SummarizedExperiment} object which contains a normalized
#' expression matrix (log2 scale) and the data frame with 'class' and
#' (optionally) variables.
#'
#' @references Michael I Love, Wolfgang Huber and Simon Anders (2014):
#' Moderated estimation of
#'  fold change and dispersion for RNA-Seq data with DESeq2. Genome Biology
#'
#' @author Mattia Chiesa, Luca Piacentini
#'
#' @seealso
#' \code{\link{varianceStabilizingTransformation}, \link{rlog} }
#'
#' @examples
#' # use example data:
#' data(SE)
#' # perform normalization on a subset of data:
#' SE_sub<-SE[1:10000, c(1:3, 21:23)]
#' data_norm <- DaMiR.normalization(SE_sub, minCounts=10, fSample=0.8,
#' hyper="yes", th.cv = 2.5)
#'
#' @export
#'
#'
DaMiR.normalization <- function(data,
                                minCounts=10,
                                fSample=0.5,
                                hyper=c("yes","no"),
                                th.cv=3,
                                type=c("vst","rlog")){
  # check missing arguments
  if (missing(data))
    stop("'data' argument must be provided")
  if (missing(hyper)){
    hyper <- hyper[1]
  }
  if (missing(type)){
    type <- type[1]
  }

  # check the type of argument
  if(!(is(data, "SummarizedExperiment")))
    stop("'data' must be a 'SummarizedExperiment' object")
  if(!(is.numeric(minCounts)))
    stop("'minCounts' must be numeric")
  if(!(is.numeric(fSample)))
    stop("'fSample' must be numeric")
  if(!(is.numeric(th.cv)))
    stop("'th.cv' must be numeric")

  # check the presence of NA or Inf
  if (any(is.na(assay(data))))
    stop("NA values are not allowed in the 'data' matrix")
  if (any(is.infinite(assay(data))))
    stop("Inf values are not allowed in the 'data' matrix")

  # specific checks
  if (any(assay(data) < 0))
    stop( "'data' contains negative values" )
  if (all(assay(data) == 0))
    stop("All genes have 0 counts.")
  if (any((assay(data) %%1) != 0))
    stop("Some values are not raw counts. Check the dataset")
  if (all(!(colnames(assay(data)) %in% rownames(colData(data)))))
    stop("colnames(assay(data)) must be equal to rownames(colData(data))")
  if(!("class" %in% colnames(colData(data))))
    stop("'class' info is lacking!
         Include the variable 'class'
         in colData(data) and label it 'class'!")

  # start execution
  init_lenght<-dim(data)[1]

  #filtering by Expression Level
  minSamplesExpr<-round((dim(data)[2])*fSample)
  data <- data[rowSums(assay(data) > minCounts) > minSamplesExpr,]
  exprs_counts <- assay(data)
  cat(init_lenght-dim(data)[1],
      "Features have been filtered out by espression.",
      dim(data)[1], "Features remained.","\n")

  #filtering by CV
  if (hyper == "yes"){
    init_lenght_cv<-dim(data)[1]

    label<-data@colData$class
    patt_cl1<-paste(c("^", levels(data@colData$class)[1],
                      "$"), sep="", collapse="")
    patt_cl2<-paste(c("^", levels(data@colData$class)[2],
                      "$"), sep="", collapse="")

    cv_value1 <- 0
    cv_value2 <- 0
    funct_cv <- function(matr_4_cv){
      cv_val<- sd(matr_4_cv)/abs(mean(matr_4_cv))
    }
    cv_value1 <- apply(exprs_counts[, grep(patt_cl1,
                                           data@colData$class)], 1, funct_cv)
    cv_value2 <- apply(exprs_counts[, grep(patt_cl2,
                                           data@colData$class)], 1, funct_cv)

    exprs_counts <- exprs_counts[cv_value1 < th.cv & cv_value2 < th.cv ,]
    # exprs_counts <- assay(data)
    cat(init_lenght_cv-dim(exprs_counts)[1],
        "'Hypervariant' Features have been filtered out.",
        dim(exprs_counts)[1], "Features remained.","\n")

  } else if(hyper == "no"){}
  else{stop("'hyper' must be set with 'yes' or 'no' ")}

  # Normalization
  if(type == "vst"){
    cat("Performing Normalization by 'vst'","\n")
    data_norm <-varianceStabilizingTransformation(exprs_counts)
  } else if (type == "rlog"){
    cat("Performing Normalization by 'rlog'.
        For large dataset it could be very time-consuming.","\n")
    data_norm <-rlog(exprs_counts)}
  else{
    stop("Please set 'vst or 'rlog' as normalization type.")
  }

  data_norm<-SummarizedExperiment(assays=as.matrix(data_norm),
                                  colData=as.data.frame(colData(data)))

  return(data_norm)
}

#' @title Filter Samples by Mean Correlation Distance Metric
#'
#' @description This function implements a sample-per-sample correlation.
#'  Samples with a mean
#' correlation lower than a user's defined threshold will be filtered out.
#'
#' @param data A SummarizedExpression object
#' @param th.corr Threshold of mean correlation; default is 0.9
#' @param type Type of correlation metric; default is "spearman"
#'
#' @details This step introduces a sample quality checkpoint. Global gene
#'  expression should,
#' in fact, exhibit a high correlation among biological replicates;
#' conversely, low correlated
#' samples may be suspected to bear some technical artifact (e.g. poor RNA
#'  or library
#' preparation quality), despite they may have passed sequencing quality
#'  checks. If not assessed,
#' these samples may, thus, negatively affect all the downstream analysis.
#'  This function looks at
#' the mean absolute correlation of each sample and removes those samples
#' with a mean correlation
#' lower than the value set in \code{th.corr} argument. This threshold may
#'  be specific for
#' different experimental setting but should be as high as possible.
#' For sequencing data we
#' suggest to set \code{th.corr} greater than 0.85.
#'
#' @return A \code{SummarizedExperiment} object which contains a normalized
#'  and filtered
#' expression matrix (log2 scale) and a filtered data frame with 'class'
#' and (optionally) variables.
#'
#' @author Mattia Chiesa, Luca Piacentini
#'
#' @examples
#' # use example data:
#' data(data_norm)
#' # filter out samples with Pearson's correlation <0.92:
#' data_filt<- DaMiR.sampleFilt(data_norm, th.corr=0.92, type ="pearson")
#'
#' @export
#'
#'
DaMiR.sampleFilt <- function(data,
                             th.corr=0.9,
                             type=c("spearman","pearson")){

  # check missing arguments
  if (missing(data)) stop("'data' argument must be provided")
  if (missing(type)){
    type <- type[1]
  }

  # check the type of argument
  if(!(is(data, "SummarizedExperiment")))
    stop("'data' must be a 'SummarizedExperiment' object")
  if(!(is.numeric(th.corr)))
    stop("'th.corr' must be numeric")

  # check the presence of NA or Inf
  if (any(is.na(assay(data))))
    stop("NA values are not allowed in the 'data' matrix")
  if (any(is.infinite(assay(data))))
    stop("Inf values are not allowed in the 'data' matrix")

  # specific checks
  if (all(assay(data) == 0))
    stop("All genes have 0 values")
  if (all((assay(data) %%1) == 0))
    warning("It seems that you are using raw counts!
            This function works with normalized data")
  if (all(!(colnames(assay(data)) %in% rownames(colData(data)))))
    stop("colnames(assay(data)) must be equal to rownames(colData(data))")
  if(!("class" %in% colnames(colData(data))))
    stop("'class' info is lacking!
         Include the variable 'class'
         in colData(data) and label it 'class'!")
  if (th.corr >1 | th.corr < 0)
    stop("'th.corr' must be between 0 and 1")

  count_data<-assay(data)


  number_of_samples<-dim(count_data)[2]
  # Sample-by-Sample correlation
  if(type == "spearman"){
    cormatrix <- abs(rcorr(as.matrix(count_data), type='spearman')$r)}
  else if(type == "pearson"){
    cormatrix <- abs(rcorr(as.matrix(count_data), type='pearson')$r)
  } else {
    stop("Please set 'spearman or 'pearson' as correlation type.")
  }

  mean_corr <- rowMeans(cormatrix)

  #check
  if(length(which(mean_corr>=th.corr)) == 0){
    stop("You are removing all samples. Please decrease 'th.corr' value.")
  } else {

  # update SummarizedExperiment object
  data_filt<-data[, mean_corr>=th.corr]

  cat(number_of_samples-dim(data_filt)[2],
      "Samples have been excluded by averaged Sample-per-Sample correlation.",
      "\n",dim(data_filt)[2], "Samples remained.","\n")
  if(number_of_samples>dim(data_filt)[2]){
  cat ("Filtered out samples :",
       "\n",
       setdiff(rownames(colData(data)), rownames(colData(data_filt))))
  }
  return(data_filt)

  }
}


