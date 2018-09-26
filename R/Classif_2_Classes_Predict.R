#' @title Predict new samples class
#'
#' @description The models learned by the
#' \link{DaMiR.EnsembleLearning2cl_Training} functionn are applied to
#'  a dataset, in order to predict the samples class
#'
#' @param data A data frame of normalized expression data.
#' Rows and Cols should be, respectively, observations and features
#' @param models_List A list with the models trained by
#' \link{DaMiR.EnsembleLearning2cl_Training} function.
#'
#' @return A dataframe containing the predictions
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
#' # DaMiR.EnsembleLearning2cl_Predict(selected_features, Tr_res)
#'
#' @export
#'
#'
DaMiR.EnsembleLearning2cl_Predict <- function(data,
                                              models_List){
  # check missing arguments
  if (missing(data))
    stop("'data' argument must be provided")
  if (missing(models_List))
    stop("'models_List' argument must be provided")

  # check the type of argument
  if(!(is.data.frame(data)))
    stop("'data' must be a data frame")

  if(!(is.list(models_List)))
    stop("'models_List' must be a list")

  # specific checks

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
      tPred[ii,which(colnames(tPred) %in% "RF")] <- unlist(predict(models_List[["RF"]],data[ii,]))
    }
    if (any(colnames(tPred) %in% "SVM")){
      tPred[ii,which(colnames(tPred) %in% "SVM")] <- unlist(predict(models_List[["SVM"]],data[ii,]))
    }
    if (any(colnames(tPred) %in% "NB")){
      tPred[ii,which(colnames(tPred) %in% "NB")] <- unlist(predict(models_List[["NB"]],data[ii,]))
    }
    if (any(colnames(tPred) %in% "LDA")){
      tPred[ii,which(colnames(tPred) %in% "LDA")] <- unlist(predict(models_List[["LDA"]],data[ii,]))[1]
    }
    if (any(colnames(tPred) %in% "LR")){
      tPred[ii,which(colnames(tPred) %in% "LR")] <- unlist(predict(models_List[["LR"]],data[ii,]))
    }
    if (any(colnames(tPred) %in% "NN")){
      tPred[ii,which(colnames(tPred) %in% "NN")] <- unlist(predict(models_List[["NN"]],data[ii,]))
    }
    if (any(colnames(tPred) %in% "PLS")){
      tPred[ii,which(colnames(tPred) %in% "PLS")] <- unlist(predict(models_List[["PLS"]],data[ii,]))
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
  tPred <- tPred[,c(idx_ens,1:(idx_ens-1)),drop=FALSE]

  tPred[tPred!=1] <- levels(models_List$classes)[2]
  tPred[tPred==1] <- levels(models_List$classes)[1]


  return(Predictions = tPred)
}
