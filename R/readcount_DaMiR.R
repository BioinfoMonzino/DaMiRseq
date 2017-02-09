#' @title Import RNA-Seq count data and covariates
#'
#' @description This function allows the user to simultaneously import counts, class (mandatory) and
#' covariate (optional) data, and creates a \code{DESeqDataSet} object.
#'
#' @param counts_file The name (with path, if needed) of a tab-delimited file which contains
#'  RNA-Seq count data. Each row is a feature (i.e. gene, transcript, exon etc.) and
#'  each column is a sample
#' @param classcov_file The name (with path, if needed) of a tab-delimited file which
#'  contains experiment information. Each row is a sample and each column is a variable.
#'  This file must contain at least one column which represent 'class' information for data
#'  adjustment and classification; the class column must be labeled as 'class'
#'
#' @return A list containing:
#' \itemize{
#'   \item A \code{DESeqDataSet} object as created by DESeqDataSetFromMatrix function of DESeq2.
#'   \item A data frame with class and covariates information.
#' }
#'
#' @references Michael I Love, Wolfgang Huber and Simon Anders (2014): Moderated estimation of
#'  fold change and dispersion for RNA-Seq data with DESeq2. Genome Biology
#'
#' @author Mattia Chiesa, Luca Piacentini
#'
#' @examples
#' rawdata.path <- system.file(package = "DaMiRseq","extdata")
#' filecounts <- list.files(rawdata.path,full.names = TRUE)[1]
#' filecovariates <- list.files(rawdata.path,full.names = TRUE)[2]
#' list.res <- DaMiR.read(filecounts,filecovariates)
#'
#' @seealso
#' \code{\link{DESeqDataSetFromMatrix}}
#'
#' @export
#'
#'
DaMiR.read <- function(counts_file, classcov_file){
  # check arguments
  if (missing(counts_file)) stop("'counts_file' argument must be provided")
  if (missing(classcov_file)) stop("'classcov_file' argument must be provided")

  datasetInput <- read.delim(counts_file)
  covariates <- read.delim(classcov_file)

  label <- as.factor(covariates$class)
  if(length(label)==0) {stop("Please set 'class' as label for class vector")}

  dds <- DESeqDataSetFromMatrix(countData = datasetInput, colData = covariates, design = ~ class)

  #class info
  patt_cl1<-paste(c("^",levels(dds@colData$class)[1],"$"),sep="",collapse="")
  patt_cl2<-paste(c("^",levels(dds@colData$class)[2],"$"),sep="",collapse="")
  l_cl1<-length(grep(patt_cl1,dds@colData$class))
  l_cl2<-length(grep(patt_cl2,dds@colData$class))

  cat("Your dataset has:","\n")
  cat(dim(dds)[1],"Features;","\n")
  cat(dim(dds)[2],"Samples:",l_cl1,levels(dds@colData$class)[1],"class","and",l_cl2,levels(dds@colData$class)[2],"class",";\n")
  cat(dim(covariates)[2],"Covariates, 'class' included.","\n")

  return(list(dds = dds, classCov = covariates))
}
