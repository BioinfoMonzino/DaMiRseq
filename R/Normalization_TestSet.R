#' @title Normalization of Independent Test Set
#'
#' @description This function aims to normalize properly an actual
#' independent test set by taking information from the Learning set
#' that will be used to transform the new sample(s).
#'
#' @param Learning_set A SummarizedExperiment object or a data frame/matrix
#' of raw count data. The learning set is supposed to be a raw counts
#' dataset of the expressed features (not all features).
#' Rows and Cols should be features and samples, respectively.
#' @param Ind_Test_set A SummarizedExperiment object or a data frame/matrix
#' of raw count data. The independent test set is supposed to be a raw counts
#' dataset with the same features of 'Learning_set'.
#' Rows and Cols should be features and samples, respectively.
#' @param normtype Type of normalization to be applied:
#' \code{varianceStabilizingTransformation}
#' (\code{vst}), \code{rlog} or \code{logcpm} are allowed;
#'  default is "\code{vst}".
#' @param method Type of method to estimate the dispersion, applied to the
#' independent test set to normalize data. Only 'precise' and 'quick' are
#' allowed. In the first case, the dispersion is estimated by the Learning
#' set and applied to the independent test set. In the second case, is
#' estimated from the independent test set. Default is "precise".
#' See details in \link{dispersionFunction}
#'
#' @details
#' The Learning_set is supposed to be a raw counts dataset of the
#' expressed features. Moreover, the independent test set is supposed
#' to be a raw counts dataset with the same features of 'Learning_set'.
#' The independent test set is normalized, taking into account the dispersion
#' parameter, estimated by the Learning set ('precise' method) or by the
#' independent test set itself ('quick' method).
#'
#' @return A matrix containing a normalized expression matrix (log2 scale)
#'
#' @references Michael I Love, Wolfgang Huber and Simon Anders (2014):
#' Moderated estimation of
#'  fold change and dispersion for RNA-Seq data with DESeq2. Genome Biology
#'
#' @author Mattia Chiesa, Luca Piacentini
#'
#' @seealso
#' \code{\link{varianceStabilizingTransformation}, \link{rlog} \link{cpm}}
#'
#' @examples
#' # use example data:
#' data(SE)
#'
#' @export
#'
#'
DaMiR.iTSnorm <- function(Learning_set,
                                Ind_Test_set,
                                normtype=c("vst","rlog","logcpm"),
                                method=c("precise","quick")){

  # check missing arguments
  if (missing(Learning_set))
    stop("'Learning_set' argument must be provided")
  if (missing(Ind_Test_set))
    stop("'Ind_Test_set' argument must be provided")
  if (missing(normtype)){
    normtype <- normtype[1]
  }
  if (missing(method)){
    method <- method[1]
  }

  # check the type of argument
  if (!(
    is(Learning_set, "SummarizedExperiment") |
    is.data.frame(Learning_set) |
    is.matrix(Learning_set))
  )
    stop("'Learning_set' must be a 'data.frame', a 'matrix'
       or a 'SummarizedExperiment' object")

  if (!(
    is(Ind_Test_set, "SummarizedExperiment") |
    is.data.frame(Ind_Test_set) |
    is.matrix(Ind_Test_set))
  )
    stop("'Ind_Test_set' must be a 'data.frame', a 'matrix'
     or a 'SummarizedExperiment' object")

  if (length(normtype) > 1)
    stop("length(normtype) must be equal to 1")
  if (!(all(normtype %in% c("vst", "rlog", "logcpm"))))
    stop("'normtype' must be 'vst', 'rlog' or 'logcpm' ")

  if (length(method) > 1)
    stop("length(method) must be equal to 1")
  if (!(all(method %in% c("precise", "quick"))))
    stop("'method' must be 'precise' or 'quick' ")

  # Specific check
  if (is(Learning_set, "SummarizedExperiment")){
    Lset <- assay(Learning_set)
    Ldf <- as.data.frame(colData(Learning_set))
  }else{
    Lset <- as.matrix(Learning_set)
    Ldf <- as.data.frame(colnames(Learning_set))
    colnames(Ldf) <- "Samplenames"
    rownames(Ldf) <- Ldf$Samplenames
  }

  if (is(Ind_Test_set, "SummarizedExperiment")){
    ITset <- assay(Ind_Test_set)
    ITdf <- as.data.frame(colData(Ind_Test_set))
  }else{
    ITset <- as.matrix(Ind_Test_set)
    ITdf <- as.data.frame(colnames(Ind_Test_set))
    colnames(ITdf) <- "Samplenames"
    rownames(ITdf) <- ITdf$Samplenames
  }
  # check the presence of NA or Inf
  if (any(is.na(Lset)))
    stop("NA values are not allowed in 'Learning_set'")
  if (any(is.infinite(as.matrix(Lset))))
    stop("Inf values are not allowed in 'Learning_set'")
  if (any(is.na(ITset)))
    stop("NA values are not allowed in 'Ind_Test_Set'")
  if (any(is.infinite(as.matrix(ITset))))
    stop("Inf values are not allowed in 'Ind_Test_Set'")

  # Lset and ITset must have the same rownames
  if(any(rownames(ITset) != rownames(Lset)))
    stop("Learning_set and Ind_Test_set must have the same rownames")

  # Sort matrices by rownames (Genes must be in the same position)
  ITset <- ITset[match(rownames(Lset), rownames(ITset)), ]

  # Lset and ITset must have the same rownames
  if(any(nrow(ITset) != nrow(Lset)))
    stop("Learning_set and Ind_Test_set must have the same features")

  # data must be raw counts
  if (any((Lset %%1) != 0))
    stop("Check 'Learning_set': some values are not raw counts.")
  if (any((ITset %%1) != 0))
    stop("Check 'Ind_test_set': some values are not raw counts.")

  # n. sample must be > 1 if normtype == VST | rlog
  if (normtype == "vst" & ncol(ITset) == 1)
    stop("With 'vst' and 'rlog' at least 2 samples must be provided")
  if (normtype == "rlog" & ncol(ITset) == 1)
    stop("With 'vst' and 'rlog' at least 2 samples must be provided")

###################### Body
if (normtype == "vst" | normtype == "rlog"){
  cat("You selected the",normtype, "normalization and the",method,"method.",
      "\n")
}else{
  cat("You selected the logcpm normalization.
  The dispersion parameter will not be estimated.", "\n" )
}


options( warn = -1 )

## Estimate the learning set dispersion
Lset_dds <- DESeqDataSetFromMatrix(countData = Lset,
                                           colData = Ldf,
                                           design = ~1)
Lset_dds <- DESeq(Lset_dds, quiet = TRUE)

## Estimate the independent test set dispersion
ITset_dds <- DESeqDataSetFromMatrix(countData = ITset,
                                             colData = ITdf,
                                             design = ~1)
ITset_dds <- DESeq(ITset_dds, quiet = TRUE)

####### Normalization

if(normtype == "vst" & method == "precise"){

  dispersionFunction(ITset_dds) <- dispersionFunction(Lset_dds)
  norm_ITset <- varianceStabilizingTransformation(ITset_dds,
                                                  blind = FALSE)
  norm_ITset <- assay(norm_ITset)

}else if(normtype == "vst" & method == "quick"){

  norm_ITset <- varianceStabilizingTransformation(ITset_dds,
                                                  blind = FALSE)
  norm_ITset <- assay(norm_ITset)

}else if(normtype == "rlog" & method == "precise"){

  dispersionFunction(ITset_dds) <- dispersionFunction(Lset_dds)
  norm_ITset <- rlog(ITset_dds, blind = FALSE)
  norm_ITset <- assay(norm_ITset)

}else if(normtype == "rlog" & method == "quick"){

  norm_ITset <- rlog(ITset_dds, blind = FALSE)
  norm_ITset <- assay(norm_ITset)

}else if(normtype == "logcpm"){

  norm_ITset <- cpm(assay(ITset_dds),log = TRUE, prior.count = 1)

}

################ Plots
## RLE
colors <- brewer.pal(12, "Set3")

plotRLE(norm_ITset,
        k=1,
        labels=TRUE,
        isLog=TRUE,
        outline=FALSE,
        col=colors[10],
        main="Relative Log Expression",
        xaxt="n",las=1)
axis(1,
     las=2,
     at=seq_len(dim(norm_ITset)[2]),
     labels=colnames(norm_ITset))

####################
## Sample by sample distribution
acc_dotplot <- melt(as.data.frame(norm_ITset),
                    measure.vars = colnames(norm_ITset))
print(ggplot(acc_dotplot, aes(value,
                              fill = variable,
                              colours= variable)) +
        geom_density(alpha=0.3) +
        facet_wrap(~variable)+
        theme(legend.position = "none") +
        ggtitle("Sample by Sample expression value distribution")
)

return(norm_ITset)

}

#' @title Batch correction of normalized Independent Test Set
#'
#' @description This function aims to perform a batch correction on
#' a normalized independent test set, exploiting the \link{ComBat}
#' function of the sva package.
#'
#' @param adj_Learning_set A SummarizedExperiment object or a
#' data frame/matrix of adjusted and normalized data, obtained by
#' the \link{DaMiR.SVadjust} function.
#' @param norm_Ind_Test_set A data frame or a matrix of normalized data.
#'  The independent test set is supposed to be already normlaized by
#'  the \link{DaMiR.iTSnorm} function
#' @param iTS_batch (Optional). A factor or a data.frame, containing
#' information regarding experimental batches of the independent test set.
#' Users can ignore this argument, if the independent test set is deemed
#' a single experimental batch.
#'
#' @details
#' The function applied a batch correction procedure to the independent test
#' set, normalized by \link{DaMiR.iTSnorm}.
#'
#' @return A matrix containing a normalized and adjusted expression matrix
#'  (log2 scale).
#'
#' @references Jeffrey T. Leek, W. Evan Johnson, Hilary S. Parker, Elana J.
#' Fertig, Andrew E. Jaffe and John D. Storey (2016).
#' sva: Surrogate Variable Analysis. R package version 3.22.0.
#'
#' @author Mattia Chiesa, Luca Piacentini
#'
#' @seealso
#' \code{\link{ComBat}}
#'
#' @examples
#' # use example data:
#' data(SE)
#'
#' @export
#'
#'
DaMiR.iTSadjust <- function(adj_Learning_set,
                            norm_Ind_Test_set,
                            iTS_batch){

  # check missing arguments
  if (missing(adj_Learning_set))
    stop("'adj_Learning_set' argument must be provided")
  if (missing(norm_Ind_Test_set))
    stop("'norm_Ind_Test_set' argument must be provided")

  # check the type of argument
  if (!(
    is(adj_Learning_set, "SummarizedExperiment") |
    is.data.frame(adj_Learning_set) |
    is.matrix(adj_Learning_set))
  )
    stop("'adj_Learning_set' must be a 'data.frame', a 'matrix'
       or a 'SummarizedExperiment' object")

  if (!(
    is.data.frame(norm_Ind_Test_set) |
    is.matrix(norm_Ind_Test_set))
  )
    stop("'norm_Ind_Test_set' must be a 'data.frame' or a 'matrix'")


  # Specific check
  if (is(adj_Learning_set, "SummarizedExperiment")){
    Lset <- assay(adj_Learning_set)
    Ldf <- as.data.frame(colData(adj_Learning_set))
  }else{
    Lset <- as.matrix(adj_Learning_set)
    Ldf <- as.data.frame(colnames(adj_Learning_set))
    colnames(Ldf) <- "Samplenames"
    rownames(Ldf) <- Ldf$Samplenames
  }

  ITset <- norm_Ind_Test_set
  # check the presence of NA or Inf
  if (any(is.na(Lset)))
    stop("NA values are not allowed in 'Learning_set'")
  if (any(is.infinite(as.matrix(Lset))))
    stop("Inf values are not allowed in 'Learning_set'")
  if (any(is.na(ITset)))
    stop("NA values are not allowed in 'Ind_Test_Set'")
  if (any(is.infinite(as.matrix(ITset))))
    stop("Inf values are not allowed in 'Ind_Test_Set'")

  # Lset and ITset must have the same rownames
  if(any(rownames(ITset) != rownames(Lset)))
    stop("Learning_set and Ind_Test_set must have the same rownames")

  # Sort matrices by rownames (Genes must be in the same position)
  ITset <- ITset[match(rownames(Lset), rownames(ITset)), ]

  # Lset and ITset must have the same rownames
  if(any(nrow(ITset) != nrow(Lset)))
    stop("Learning_set and Ind_Test_set must have the same features")

  # batch information
  if(missing(iTS_batch)){
    batch <- as.data.frame(c(rep("b_a_t_c_h_LS_19041986", dim(Lset)[2]),
                             rep("batch_2", dim(ITset)[2])))
  }else{
    iTS_batch <- as.data.frame(iTS_batch)
    colnames(iTS_batch) <- "batch"

    if(dim(iTS_batch)[1] != dim(ITset)[2])
      stop("dim(iTS_batch)[1] must be equal to ncol(norm_Ind_Test_set)")

    LS_batch<- as.data.frame(c(rep("b_a_t_c_h_LS_19041986", dim(Lset)[2])))
    colnames(LS_batch) <- "batch"

    batch <- rbind(LS_batch,iTS_batch)
  }



  ############################ Body

  colnames(batch) <- "batch"
  data_tot <- rbind(t(Lset), t(ITset)) # data adjust

  # use ComBat method to adjust for known batch:
  modcombat <- model.matrix(~1, data=batch)
suppressMessages(combat_edata <- ComBat(dat=t(data_tot),
                         batch=batch$batch,
                         mod=modcombat,
                         par.prior=TRUE,
                         prior.plots=FALSE))

  adj_norm_Ind_test_set <- combat_edata[,-which(batch$batch %in% "b_a_t_c_h_LS_19041986")]


  ################ Plots
  ## RLE
  colors <- brewer.pal(12, "Set3")

  plotRLE(adj_norm_Ind_test_set,
          k=1,
          labels=TRUE,
          isLog=TRUE,
          outline=FALSE,
          col=colors[10],
          main="Relative Log Expression",
          xaxt="n",las=1)
  axis(1,
       las=2,
       at=seq_len(dim(adj_norm_Ind_test_set)[2]),
       labels=colnames(adj_norm_Ind_test_set))

  ####################
  ## Sample by sample distribution
  acc_dotplot <- melt(as.data.frame(adj_norm_Ind_test_set),
                      measure.vars = colnames(adj_norm_Ind_test_set))
  print(ggplot(acc_dotplot, aes(value,
                                fill = variable,
                                colours= variable)) +
          geom_density(alpha=0.3) +
          facet_wrap(~variable)+
          theme(legend.position = "none") +
          ggtitle("Sample by Sample expression value distribution")
  )



  return(adj_norm_Ind_test_set)

}


