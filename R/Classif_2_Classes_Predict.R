#' @title Predict new samples class
#'
#' @description The best model learned by the
#' \link{DaMiR.EnsL_Train} functionn is tested on
#'  a new dataset, in order to predict the samples class
#'
#' @param data A SummarizedExperiment object or a data frame/matrix
#' of normalized expression data. Rows and Cols should be
#' observations and features, respectively.
#' @param bestModel The best model, selected between those trained by
#' the \link{DaMiR.EnsL_Train} function.
#'
#' @return A matrix containing the predictions
#'
#' @details
#' This function implements the prediction step on new data,
#' given a model learned by \link{DaMiR.EnsL_Train}
#'
#' @author Mattia Chiesa, Luca Piacentini
#'
#' @examples
#' # use example data:
#' data(selected_features)
#' data(df)
#'
#' @export
#'
#'
DaMiR.EnsL_Predict <- function(data,
                               bestModel){
  # check missing arguments
  if (missing(data))
    stop("'data' argument must be provided")
  if (missing(bestModel))
    stop("'bestModel' argument must be provided")


  # check the type of argument
  if (!(
    is(data, "SummarizedExperiment") | is.data.frame(data) | is.matrix(data))
  )
    stop("'data' must be a 'data.frame', a 'matrix'
         or a 'SummarizedExperiment' object")

  if (is(data, "SummarizedExperiment")){
    data <- t(assay(data))

  }

  data <- as.data.frame(data)


  if(!(is.list(bestModel)))
    stop("'bestModel' must be a list")

  # specific checks

  # check the presence of NA or Inf
  if (any(is.na(data)))
    stop("NA values are not allowed in the 'data' matrix")
  if (any(is.infinite(as.matrix(data))))
    stop("Inf values are not allowed in the 'data' matrix")

  options( warn = -1 )

  ## body

  if (!(all(names(bestModel) %in% c("RF",
                                      "SVM",
                                      "LDA",
                                      "LR",
                                      "NB",
                                      "NN",
                                      "PLS",
                                      "Ensemble",
                                      "classes",
                                      "positiveClass"))))
    stop("'names(bestModel)' must be
         'RF','SVM','LDA','LR','NB','NN','PLS','Ensemble',
         'classes' or 'positiveClass'")


  # check nomi variabili modelli = nomi colonne testset
  var_model <- attr(bestModel[[1]][[1]][["terms"]],"term.labels")
  var_testset <- colnames(data)
  if(!(all(var_model %in% var_testset)))
    stop(" 'data' and 'bestModel' must have the same features.
         Please, remove irrelevant features from 'data'")

  # start test
  tPred <- matrix(nrow = dim(data)[1],
                  ncol = (length(bestModel)-2))
  colnames(tPred) <- names(bestModel)[1:ncol(tPred)]

  for (ii in seq_len(dim(data)[1])){

    if (any(colnames(tPred) %in% "RF")){
      tPred[ii,which(colnames(tPred) %in% "RF")] <- unlist(predict(bestModel[["RF"]],data[ii,]))
    }
    if (any(colnames(tPred) %in% "SVM")){
      tPred[ii,which(colnames(tPred) %in% "SVM")] <- unlist(predict(bestModel[["SVM"]],data[ii,]))
    }
    if (any(colnames(tPred) %in% "NB")){
      tPred[ii,which(colnames(tPred) %in% "NB")] <- unlist(predict(bestModel[["NB"]],data[ii,]))
    }
    if (any(colnames(tPred) %in% "LDA")){
      tPred[ii,which(colnames(tPred) %in% "LDA")] <- unlist(predict(bestModel[["LDA"]],data[ii,]))[1]
    }
    if (any(colnames(tPred) %in% "LR")){
      tPred[ii,which(colnames(tPred) %in% "LR")] <- unlist(predict(bestModel[["LR"]],data[ii,]))
    }
    if (any(colnames(tPred) %in% "NN")){
      tPred[ii,which(colnames(tPred) %in% "NN")] <- unlist(predict(bestModel[["NN"]],data[ii,]))
    }
    if (any(colnames(tPred) %in% "PLS")){
      tPred[ii,which(colnames(tPred) %in% "PLS")] <- unlist(predict(bestModel[["PLS"]],data[ii,]))
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
              unlist(bestModel[["Ensemble"]])*tPred[ii,
                                                      seq_len(dim(
                                                        tPred)[2]-1)]))
  }
  #Ensemble must be the first column
  class_labels <- colnames(tPred)
  idx_ens <- which(class_labels %in% "Ensemble")
  tPred <- tPred[,c(idx_ens,1:(idx_ens-1)),drop=FALSE]

  tPred[tPred!=1] <- levels(bestModel$classes)[2]
  tPred[tPred==1] <- levels(bestModel$classes)[1]


  return(Predictions = tPred)
}
