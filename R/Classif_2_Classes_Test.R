#' @title Test Binary Classifiers
#'
#' @description This function tests the models learned by the
#' \link{DaMiR.EnsL_Train} function, on a test
#' set
#'
#' @param data A SummarizedExperiment object or a data frame/matrix
#' of normalized expression data. Rows and Cols should be
#' observations and features, respectively.
#' @param classes A class vector with \code{nrow(data)} elements.
#'  Each element represents the class label for each observation.
#'  Two different class labels are allowed. Note. this argument should
#'  not be set when 'data' is a SummarizedExperiment object
#' @param EnsL_model A list with the models trained by
#' \link{DaMiR.EnsL_Train} function.
#'
#' @return A dataframe containing the predictions on the testset
#'
#' @details
#' This function implements the test step of
#' \link{DaMiR.EnsembleLearning2cl} function
#'
#' @author Mattia Chiesa, Luca Piacentini
#'
#' @examples
#' # use example data:
#' data(selected_features)
#' data(df)
#' set.seed(1)
#' # only for the example:
#' # speed up the process setting a low 'iter' argument value;
#' # for real data set use default 'iter' value (i.e. 100) or higher:
#' #  Tr_res <- DaMiR.EnsL_Train(
#' #  selected_features,classes=df$class, fSample.tr.w=0.6, iter=3,
#' # cl_type=c("RF","LR"))
#' # DaMiR.EnsembleLearning2cl_Test(selected_features,
#' #classes=df$class,Tr_res)
#'
#' @export
#'
#'
DaMiR.EnsL_Test <- function(data,
                            classes,
                            EnsL_model){


  # check missing arguments
  if (missing(data))
    stop("'data' argument must be provided")
  # if (missing(classes))
  #   stop("'classes' argument must be provided")
  if (missing(EnsL_model))
    stop("'EnsL_model' argument must be provided")

  # check the type of argument
  if (!(
    is(data, "SummarizedExperiment") | is.data.frame(data) | is.matrix(data))
  )
    stop("'data' must be a 'data.frame', a 'matrix'
         or a 'SummarizedExperiment' object")

  if (is(data, "SummarizedExperiment")){

    if(!("class" %in% colnames(colData(data))))
      stop("'class' info is lacking! Include the variable 'class'
         in colData(data) and label it 'class'!")

    classes <- colData(data)$class
    data <- t(assay(data))

  }else{
    if(missing(classes))
      stop("'classes' argument must be provided when
           data is a 'data.frame', or a 'matrix'")
  }

  data <- as.data.frame(data)

   if(!(is.factor(classes)))
     stop("'classes' must be a factor")

  if(!(is.list(EnsL_model)))
    stop("'EnsL_model' must be a list")

  # specific checks

  if(length(classes) != dim(data)[1])
    stop("length(classes) must be equal to dim(data)[1]")

  # check the presence of NA or Inf
  if (any(is.na(data)))
    stop("NA values are not allowed in the 'data' matrix")
  if (any(is.infinite(as.matrix(data))))
    stop("Inf values are not allowed in the 'data' matrix")

  options( warn = -1 )

  ## body

  if (!(all(names(EnsL_model) %in% c("RF",
                                   "SVM",
                                   "LDA",
                                   "LR",
                                   "NB",
                                   "NN",
                                   "PLS",
                                   "Ensemble",
                                   "classes",
                                   "positiveClass"))))
    stop("'names(EnsL_model)' must be
         'RF','SVM','LDA','LR','NB','NN','PLS','Ensemble',
         'classes' or 'positiveClass'")


  # check nomi variabili modelli = nomi colonne testset
   var_model <- attr(EnsL_model[[1]][[1]][["terms"]],"term.labels")
   var_testset <- colnames(data)
   if(!(all(var_model %in% var_testset)))
     stop(" 'data' and 'EnsL_model' must have the same features")

  # start test
   tPred <- matrix(nrow = dim(data)[1],
                   ncol = (length(EnsL_model)-2))
   colnames(tPred) <- names(EnsL_model)[1:ncol(tPred)]

    for (ii in seq_len(dim(data)[1])){

      if (any(colnames(tPred) %in% "RF")){
        if(unlist(predict(EnsL_model[["RF"]],data[ii,])) == levels(classes)[1])
        {tPred[ii,which(colnames(tPred) %in% "RF")] <- 1 }
        else{tPred[ii,which(colnames(tPred) %in% "RF")] <- 0 }
      }
      if (any(colnames(tPred) %in% "SVM")){
        if(unlist(predict(EnsL_model[["SVM"]],data[ii,])) == levels(classes)[1])
        {tPred[ii,which(colnames(tPred) %in% "SVM")] <- 1 }
        else{tPred[ii,which(colnames(tPred) %in% "SVM")] <- 0 }
      }
      if (any(colnames(tPred) %in% "NB")){
        if(unlist(predict(EnsL_model[["NB"]],data[ii,])) == levels(classes)[1])
        {tPred[ii,which(colnames(tPred) %in% "NB")] <- 1 }
        else{tPred[ii,which(colnames(tPred) %in% "NB")] <- 0 }
      }
      if (any(colnames(tPred) %in% "LDA")){
        #if(predict(EnsL_model[["LDA"]],data[ii,])$class == levels(classes)[1])
        if(unlist(predict(EnsL_model[["LDA"]],data[ii,]))[1] == 1)
        {tPred[ii,which(colnames(tPred) %in% "LDA")] <- 1 }
        else{tPred[ii,which(colnames(tPred) %in% "LDA")] <- 0 }
      }
      if (any(colnames(tPred) %in% "LR")){
        if(unlist(predict(EnsL_model[["LR"]],data[ii,])) == levels(classes)[1])
        {tPred[ii,which(colnames(tPred) %in% "LR")] <- 1 }
        else{tPred[ii,which(colnames(tPred) %in% "LR")] <- 0 }
      }
      if (any(colnames(tPred) %in% "NN")){
        if(unlist(predict(EnsL_model[["NN"]],data[ii,])) == levels(classes)[1])
        {tPred[ii,which(colnames(tPred) %in% "NN")] <- 1 }
        else{tPred[ii,which(colnames(tPred) %in% "NN")] <- 0 }
      }
      if (any(colnames(tPred) %in% "PLS")){
        if(unlist(predict(EnsL_model[["PLS"]],data[ii,])) == levels(classes)[1])
        {tPred[ii,which(colnames(tPred) %in% "PLS")] <- 1 }
        else{tPred[ii,which(colnames(tPred) %in% "PLS")] <- 0 }
      }
      # if (any(colnames(tPred) %in% "kNN")){
      #
      #   if (kknn(formula_DM,
      #            trainingSet_DM,
      #            testSet_DM[ii,],
      #            k=3)$CL[,3] == levels(classes)[1])
      #   {tPred[ii,which(colnames(tPred) %in% "kNN")] <- 1 }
      #   else{tPred[ii,which(colnames(tPred) %in% "kNN")] <- 0 }
      # }
      #

      # # ensemble Learning
      tPred[ii,
            which(
              colnames(tPred) %in% "Ensemble")] <- round(sum(
                unlist(EnsL_model[["Ensemble"]])*tPred[ii,
                                                     seq_len(dim(
                                                       tPred)[2]-1)]
                )
                )
    }
    #Ensemble must be the first column
    class_labels <- colnames(tPred)
    idx_ens <- which(class_labels %in% "Ensemble")
    tPred <- tPred[, c(idx_ens,1:(idx_ens-1)), drop=FALSE]

    # calculate scores

    acc.Class <- matrix(nrow=1, ncol = length(colMeans(tPred)))
    MCC.Class <- acc.Class
    TP.Class <- acc.Class
    TN.Class <- acc.Class
    FP.Class <- acc.Class
    FN.Class <- acc.Class
    Sensit.Class <- acc.Class
    Specif.Class <- acc.Class
    PPV.Class <- acc.Class
    NPV.Class <- acc.Class

    # TP,FN, FP,TN
    TP_Class <- colSums(tPred[
      which(classes == levels(classes)[1]),,drop=FALSE] == 1)
    TN_Class <- colSums(tPred[
      which(classes == levels(classes)[2]),,drop=FALSE] == 0)
    FP_Class <- colSums(tPred[
      which(classes == levels(classes)[1]),,drop=FALSE] == 0)
    FN_Class <- colSums(tPred[
      which(classes == levels(classes)[2]),,drop=FALSE] == 1)

    TP.Class <- TP_Class
    TN.Class <- TN_Class
    FP.Class <- FP_Class
    FN.Class <- FN_Class

    acc.Class <- (TP_Class + TN_Class)/(
      TP_Class + TN_Class +  FP_Class + FN_Class)

    MCC.Class <- (TP_Class * TN_Class - FP_Class * FN_Class) /
      sqrt((TP_Class + FP_Class) * (TP_Class + FN_Class) *
             (TN_Class + FP_Class) * (TN_Class + FN_Class))


    Sensit.Class <- TP_Class / (TP_Class + FN_Class)
    Specif.Class <- TN_Class / (TN_Class + FP_Class)

    PPV.Class <- TP_Class / (TP_Class + FP_Class)
    NPV.Class <- TN_Class / (TN_Class + FN_Class)

    acc.Class[which(is.nan(acc.Class))] <- 0
    MCC.Class[which(is.nan(MCC.Class))] <- 0
    Sensit.Class[which(is.nan(Sensit.Class))] <- 0
    Specif.Class[which(is.nan(Specif.Class))] <- 0
    PPV.Class[which(is.nan(PPV.Class))] <- 0
    NPV.Class[which(is.nan(NPV.Class))] <- 0
  # cat("Accuracy [%]:",
  #     "\n",
  #     colnames(tPred),
  #     "\n",
  #     "Mean:",round(acc.Class,2),"\n")
  # cat("MCC score:",
  #     "\n",
  #     colnames(tPred),
  #     "\n",
  #     "Mean:",round(MCC.Class,2),"\n")
  # cat("Sensitivity:",
  #     "\n",
  #     colnames(tPred),
  #     "\n",
  #     "Mean:",round(Sensit.Class,2),"\n")
  # cat("Specificity:",
  #     "\n",
  #     colnames(tPred),
  #     "\n",
  #     "Mean:",round(Specif.Class,2),"\n")
  #  cat("PPV:",
  #      "\n",
  #      colnames(tPred),
  #      "\n",
  #      "Mean:",round(PPV.Class,2),"\n")
  #  cat("NPV:",
  #      "\n",
  #      colnames(tPred),
  #      "\n",
  #      "Mean:",round(NPV.Class,2),"\n")

  tPred[tPred==1] <- levels(classes)[1]
  tPred[tPred==0] <- levels(classes)[2]

  output_list <- list()
  output_list$Prediction <- tPred
  output_list$accuracy <- acc.Class
  output_list$MCC <- MCC.Class
  output_list$sensitivity <- Sensit.Class
  output_list$Specificty <- Specif.Class
  output_list$PPV <- PPV.Class
  output_list$NPV <- NPV.Class

  output_list$positiveClass <- levels(classes)[1]

  output_list$predictors <- colnames(data)

  return(output_list = output_list)
}
