#' @title Test Binary Classifiers
#'
#' @description This function tests the models learned by the
#' \link{DaMiR.EnsembleLearning2cl_Training} function, on a test set
#'
#' @param data A data frame of normalized expression data.
#' Rows and Cols should be, respectively, observations and features
#' @param classes A class vector with \code{nrow(data)} elements.
#'  Each element represents the class label for each observation.
#'  Two different class labels are allowed
#' @param models_List A list with the models trained by
#' \link{DaMiR.EnsembleLearning2cl_Training} function.
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
#' #  Tr_res <- DaMiR.EnsembleLearning2cl_Training(
#' #  selected_features,classes=df$class, fSample.tr.w=0.6, iter=3,
#' # cl_type=c("RF","LR"))
#' # DaMiR.EnsembleLearning2cl_Test(selected_features, classes=df$class,
#' # Tr_res)
#'
#' @export
#'
#'
DaMiR.EnsembleLearning2cl_Test <- function(data,
                                               classes,
                                               models_List){
  # check missing arguments
  if (missing(data))
    stop("'data' argument must be provided")
  if (missing(classes))
    stop("'classes' argument must be provided")
  if (missing(models_List))
    stop("'models_List' argument must be provided")

  # check the type of argument
   if(!(is.data.frame(data)))
     stop("'data' must be a data frame")

   if(!(is.factor(classes)))
     stop("'classes' must be a factor")

  if(!(is.list(models_List)))
    stop("'models_List' must be a list")

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

  if (!(all(names(models_List) %in% c("RF",
                                   "SVM",
                                   "LDA",
                                   "LR",
                                   "NB",
                                   "NN",
                                   "PLS",
                                   "Ensemble",
                                   "classes"))))
    stop("'names(models_List)' must be
         'RF','SVM','LDA','LR','NB','NN','PLS','Ensemble' or 'classes'")


  # check nomi variabili modelli = nomi colonne testset
   var_model <- attr(models_List[[1]][[1]][["terms"]],"term.labels")
   var_testset <- colnames(data)
   if(!(all(var_model %in% var_testset)))
     stop(" 'data' and 'models_List' must have the same features")

  # start test
   tPred <- matrix(nrow = dim(data)[1],
                   ncol = (length(models_List)-1))
   colnames(tPred) <- names(models_List)[-length(names(models_List))]

    for (ii in seq_len(dim(data)[1])){

      if (any(colnames(tPred) %in% "RF")){
        if(unlist(predict(models_List[["RF"]],data[ii,])) == levels(classes)[1])
        {tPred[ii,which(colnames(tPred) %in% "RF")] <- 1 }
        else{tPred[ii,which(colnames(tPred) %in% "RF")] <- 0 }
      }
      if (any(colnames(tPred) %in% "SVM")){
        if(unlist(predict(models_List[["SVM"]],data[ii,])) == levels(classes)[1])
        {tPred[ii,which(colnames(tPred) %in% "SVM")] <- 1 }
        else{tPred[ii,which(colnames(tPred) %in% "SVM")] <- 0 }
      }
      if (any(colnames(tPred) %in% "NB")){
        if(unlist(predict(models_List[["NB"]],data[ii,])) == levels(classes)[1])
        {tPred[ii,which(colnames(tPred) %in% "NB")] <- 1 }
        else{tPred[ii,which(colnames(tPred) %in% "NB")] <- 0 }
      }
      if (any(colnames(tPred) %in% "LDA")){
        #if(predict(models_List[["LDA"]],data[ii,])$class == levels(classes)[1])
        if(unlist(predict(models_List[["LDA"]],data[ii,]))[1] == 1)
        {tPred[ii,which(colnames(tPred) %in% "LDA")] <- 1 }
        else{tPred[ii,which(colnames(tPred) %in% "LDA")] <- 0 }
      }
      if (any(colnames(tPred) %in% "LR")){
        if(unlist(predict(models_List[["LR"]],data[ii,])) == levels(classes)[1])
        {tPred[ii,which(colnames(tPred) %in% "LR")] <- 1 }
        else{tPred[ii,which(colnames(tPred) %in% "LR")] <- 0 }
      }
      if (any(colnames(tPred) %in% "NN")){
        if(unlist(predict(models_List[["NN"]],data[ii,])) == levels(classes)[1])
        {tPred[ii,which(colnames(tPred) %in% "NN")] <- 1 }
        else{tPred[ii,which(colnames(tPred) %in% "NN")] <- 0 }
      }
      if (any(colnames(tPred) %in% "PLS")){
        if(unlist(predict(models_List[["PLS"]],data[ii,])) == levels(classes)[1])
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
                unlist(models_List[["Ensemble"]])*tPred[ii,
                                                     seq_len(dim(
                                                       tPred)[2]-1)]))
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
  cat("Accuracy [%]:",
      "\n",
      colnames(tPred),
      "\n",
      "Mean:",round(acc.Class,2),"\n")
  cat("MCC score:",
      "\n",
      colnames(tPred),
      "\n",
      "Mean:",round(MCC.Class,2),"\n")
  cat("Sensitivity:",
      "\n",
      colnames(tPred),
      "\n",
      "Mean:",round(Sensit.Class,2),"\n")
  cat("Specificity:",
      "\n",
      colnames(tPred),
      "\n",
      "Mean:",round(Specif.Class,2),"\n")
   cat("PPV:",
       "\n",
       colnames(tPred),
       "\n",
       "Mean:",round(PPV.Class,2),"\n")
   cat("NPV:",
       "\n",
       colnames(tPred),
       "\n",
       "Mean:",round(NPV.Class,2),"\n")

  tPred[tPred==1] <- levels(classes)[1]
  tPred[tPred==0] <- levels(classes)[2]

  return(Predictions = tPred)
}
