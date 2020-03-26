#' @title Train a Binary Classifier using 'Staking' Learning strategy.
#'
#' @description This function learn a meta learner by a 'Stacking'
#' strategy.
#' Users can provide heterogeneous features
#' (other than genomic features)
#' which will be taken into account during
#' classification model building. A 'two-classes' classification task
#'  isaddressed.
#'
#' @param data A SummarizedExperiment object or a data frame/matrix
#' of normalized expression data. Rows and Cols should be
#' observations and features, respectively.
#' @param classes A class vector with \code{nrow(data)} elements.
#'  Each element represents the class label for each observation.
#'  Two different class labels are allowed. Note. this argument should
#'  not be set when 'data' is a SummarizedExperiment object
#' @param variables An optional data frame containing other variables
#' (but without 'class' column). Each column represents a different
#' covariate to be considered in the model
#' @param fSample.tr.w Fraction of samples of training set to be used
#' during weight estimation; default is 0.7
#' @param cl_type List of weak classifiers that will compose the
#' meta-learners. "RF", "SVM", "LDA", "LR", "NB", "NN", "PLS"
#' are allowed. Default is c("RF", "LR", "LDA", "NB", "SVM")
#'
#' @return A list containing:
#' \itemize{
#'   \item The models of each classifier used to build the Ensemble
#'   meta-learner with the median or the best accuracy
#'   (over the iteration) for the Ensemble classifier;
#'   \item the weights associated to each weak classifier;
#' }
#'
#' @details
#' This function implements the training step of
#' \link{DaMiR.EnsembleLearning2cl} function
#'
#' @author Mattia Chiesa, Luca Piacentini
#'
#' @examples
#' # use example data:
#' data(selected_features)
#' data(df)
#' set.seed(1)
#' # For the example:
#' # speed up the process setting a low 'iter' argument value;
#' # for real data set use default 'iter' value (i.e. 100) or higher:
#' #  Classification_res <- DaMiR.EnsL_Train(
#' #  selected_features,classes=df$class, fSample.tr.w=0.6, iter=3,
#' # cl_type=c("RF","LR"))
#'
#' @export
#'
#'
DaMiR.EnsL_Train <- function(data,
                             classes,
                             variables,
                             fSample.tr.w=0.7,
                             cl_type=c("RF",
                                       "SVM",
                                       "LDA",
                                       "LR",
                                       "NB",
                                       "NN",
                                       "PLS")){

  # check missing arguments
  if (missing(data))
    stop("'data' argument must be provided")
  # if (missing(classes))
  #   stop("'classes' argument must be provided")
  if (missing(cl_type)){
    cl_type <- c("RF", "LR", "LDA", "NB", "SVM")
  }

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

 # check the type of argument
  if(!(is.numeric(fSample.tr.w)))
    stop("'fSample.tr.w' must be numeric")

  if(!(is.factor(classes)))
    stop("'classes' must be a factor")


  # specific checks
  if (fSample.tr.w >0.9 | fSample.tr.w < 0.5)
    stop("'th.corr' must be between 0.5 and 1")

  if((dim(data)[1]-round(dim(data)[1]*fSample.tr.w)) == 0)
    stop("A Test Set is not available to weight estimation
         Decrease 'fSample.tr.w' or increase the number of
         observation.")
  if(length(classes) != dim(data)[1])
    stop("length(classes) must be equal to dim(data)[1]")

  # check balanced dataset (Thanks to Dr. Pawel Karpinski)
  class_lev_ckeck <- levels(classes)
  count_cl1 <- length(which(classes %in% class_lev_ckeck[1]))
  count_cl2 <- length(which(classes %in% class_lev_ckeck[2]))
  # if(count_cl1 != count_cl2)
  #   stop("Trainingset must be balanced (same n. of samples in
  #        both classes)")

  if (missing(variables)){
    data <- data
  } else {
    variables<-as.data.frame(variables)
    if(!(is.data.frame(variables)))
      stop("'variables' must be a data frame") ###

    for (ic in seq_len(dim(variables)[2])){
      if(isTRUE(is.factor(variables[,ic]))){
        variables[,ic]<-as.numeric(variables[,ic])
      }
    }
    data <- cbind(data, variables)
  }

  # check the presence of NA or Inf
  if (any(is.na(data)))
    stop("NA values are not allowed in the 'data' matrix")
  if (any(is.infinite(as.matrix(data))))
    stop("Inf values are not allowed in the 'data' matrix")

  options( warn = -1 )

  ## body

  if (!(all(cl_type %in% c("RF", "SVM","LDA","LR","NB","NN","PLS")))
      )
    stop("'cl_type' must be
         'RF', 'SVM', 'LDA', 'LR', 'NB', 'NN', 'PLS'")

  # training and test
  # cat("You select:",cl_type, "weak classifiers for creating
  #     the Ensemble meta-learner.","\n")



  acc.Class<- matrix()

  wei_all <- matrix(nrow = 1, ncol = length(cl_type))
  model_rf_all <- list()
  model_svm_all <- list()
  model_nb_all <- list()
  model_lr_all <- list()
  model_nn_all <- list()
  model_lda_all <- list()
  model_pls_all <- list()


  # start training


    # Sample Training for weighting
    sample_index_train2 <- ""
    class_level2 <- levels(classes)
    tr_sample2 <- round(
      dim(data)[1]*fSample.tr.w/length(class_level2))
    for (i in seq_len(length(class_level2))){
      sample_index2 <- sample(which(class_level2[i]==classes),
                              replace = FALSE)[seq_len(tr_sample2)]
      sample_index_train2 <- c(sample_index_train2,sample_index2)
    }
    sample_index_train2 <- as.numeric(sample_index_train2[-1])

    # create TrainingSet for weighting
    trainingSet2 <- data[sample_index_train2,, drop=FALSE]
    trainingSetClasses2 <- classes[sample_index_train2,
                                              drop=FALSE]
    trainingSetClasses2 <- droplevels(trainingSetClasses2)

    # create TestSet for weighting
    sample_index_test2 <- setdiff(seq_len(dim(data)[1]),
                                  sample_index_train2)
    testSet2 <- data[sample_index_test2,, drop=FALSE]
    testSetClasses2 <- classes[sample_index_test2,
                                          drop=FALSE]
    testSetClasses2 <- droplevels(testSetClasses2)


    # create the formula and datasets for models
    trainingSet_DM <- cbind(trainingSet2,trainingSetClasses2)
    varNames <- colnames(trainingSet2)
    colnames(trainingSet_DM) <- c(varNames, "classes")
    varNames1 <- paste(varNames, collapse = "+")
    formula_DM <- as.formula(paste("classes", varNames1, sep = " ~ ")
                             )

    testSet_DM <- cbind(testSet2,testSetClasses2)
    varNames <- colnames(testSet2)
    colnames(testSet_DM) <- c(varNames, "classes")


    ####################################################
    ###############  Weak classifiers ##################
    ####################################################
    model_class <- list()
    acc_model <- c()
    kk <- 1

    # Random Forest
    if (any(cl_type %in% "RF")){
      model_rf <- randomForest(formula = formula_DM,
                               data = trainingSet_DM,
                               ntree=1000,
                               importance =TRUE)
      acc_model[kk] <- caret::confusionMatrix(
        table(predict(model_rf, testSet_DM),testSet_DM$classes),
        reference = testSet_DM$classes)$overall['Accuracy']
      names(acc_model)[kk] <- "RF"

      model_class[[kk]] <- model_rf
      kk <- kk +1

    }

    # Support Vector Machine
    if (any(cl_type %in% "SVM")){

      tune.for.model <- tune.svm(formula_DM,
                                 data = trainingSet_DM,
                                 gamma = 10^(-6:-1),
                                 cost = 10^(1:4),
                                 tunecontrol = tune.control(
                                   # cross = min(min_sub,10) ))
                                   cross = min(dim(trainingSet_DM)[1],
                                               10) ))
      model_svm <- svm(formula = formula_DM,
                       data = trainingSet_DM,
                       gamma = tune.for.model$best.parameters$gamma,
                       cost = tune.for.model$best.parameters$cost,
                       hidden=3)
      acc_model[kk] <- caret::confusionMatrix(
        predict(model_svm, testSet_DM),
        reference = testSet_DM$classes)$overall['Accuracy']
      names(acc_model)[kk] <- "SVM"

      model_class[[kk]] <- model_svm
      kk <- kk +1

    }

    # Naive Bayes
    if (any(cl_type %in% "NB")){

      model_nb <- naiveBayes(formula = formula_DM,
                             data = trainingSet_DM)
      acc_model[kk] <- caret::confusionMatrix(
        predict(model_nb, testSet_DM),
        reference = testSet_DM$classes)$overall['Accuracy']
      names(acc_model)[kk] <- "NB"

      model_class[[kk]] <- model_nb
      kk <- kk +1
    }

    # Linear Discriminant Analysis
    if (any(cl_type %in% "LDA")){
      model_lda <- lda(formula = formula_DM,
                       data = trainingSet_DM)
      acc_model[kk] <- caret::confusionMatrix(
        predict(model_lda, testSet_DM)$class,
        reference = testSet_DM$classes)$overall['Accuracy']
      names(acc_model)[kk] <- "LDA"

      model_class[[kk]] <- model_lda
      kk <- kk +1
    }

    # bayesian Logistic Regression
    if (any(cl_type %in% "LR")){
      model_lr <- caret::train(formula_DM,
                               method="bayesglm",
                               data = trainingSet_DM,
                               family=binomial(link="logit"))
      acc_model[kk] <- caret::confusionMatrix(
        predict(model_lr, testSet_DM),
        reference = testSet_DM$classes)$overall['Accuracy']
      names(acc_model)[kk] <- "LR"

      model_class[[kk]] <- model_lr
      kk <- kk +1
    }

    # Neural Networks
    if (any(cl_type %in% "NN")){
      model_nn <- caret::train(formula_DM,
                               method="mlpWeightDecay",
                               data = trainingSet_DM)

      acc_model[kk] <- caret::confusionMatrix(
        predict(model_nn, testSet_DM),
        reference = testSet_DM$classes)$overall['Accuracy']
      names(acc_model)[kk] <- "NN"

      model_class[[kk]] <- model_nn
      kk <- kk +1

    }

    # PLS
    if (any(cl_type %in% "PLS")){
      model_pls <- caret::train(formula_DM,
                                method="pls",
                                data = trainingSet_DM)
      acc_model[kk] <- caret::confusionMatrix(
        predict(model_pls, testSet_DM),
        reference = testSet_DM$classes)$overall['Accuracy']
      names(acc_model)[kk] <- "PLS"

      model_class[[kk]] <- model_pls

      #kk <- kk +1
    }

    # # k-Nearest Neighbours
    # if (any(cl_type %in% "kNN")){
    #
    #   acc_model[kk] <- caret::confusionMatrix(
    #     as.factor(kknn(formula_DM, trainingSet_DM, testSet_DM,
    #k=3)$CL[,3]),
    #     reference = testSet_DM$classes)$overall['Accuracy']
    #   names(acc_model)[kk] <- "kNN"
    #
    # }

    # Weighting
    wei <- acc_model/sum(acc_model)
    wei_all <- wei

    ###############################################
    ## Test Samples
    testSetClasses <- testSet_DM$classes
    tPred <- matrix(nrow = dim(testSet_DM)[1],
                    ncol = (length(wei) + 1))
    colnames(tPred) <- c("Ensemble",names(wei))

    for (ii in seq_len(dim(testSet_DM)[1])){

      if (any(colnames(tPred) %in% "RF")){
        if(predict(model_rf,testSet_DM[ii,]) == levels(testSetClasses)[1])
        {tPred[ii,which(colnames(tPred) %in% "RF")] <- 1 }
        else{tPred[ii,which(colnames(tPred) %in% "RF")] <- 0 }
      }
      if (any(colnames(tPred) %in% "SVM")){
        if(predict(model_svm,testSet_DM[ii,]) == levels(testSetClasses)[1])
        {tPred[ii,which(colnames(tPred) %in% "SVM")] <- 1 }
        else{tPred[ii,which(colnames(tPred) %in% "SVM")] <- 0 }
      }
      if (any(colnames(tPred) %in% "NB")){
        if(predict(model_nb,testSet_DM[ii,]) == levels(testSetClasses)[1])
        {tPred[ii,which(colnames(tPred) %in% "NB")] <- 1 }
        else{tPred[ii,which(colnames(tPred) %in% "NB")] <- 0 }
      }
      if (any(colnames(tPred) %in% "LDA")){
        if(predict(model_lda,testSet_DM[ii,])$class == levels(testSetClasses)[1])
        {tPred[ii,which(colnames(tPred) %in% "LDA")] <- 1 }
        else{tPred[ii,which(colnames(tPred) %in% "LDA")] <- 0 }
      }
      if (any(colnames(tPred) %in% "LR")){
        if(predict(model_lr,testSet_DM[ii,]) == levels(testSetClasses)[1])
        {tPred[ii,which(colnames(tPred) %in% "LR")] <- 1 }
        else{tPred[ii,which(colnames(tPred) %in% "LR")] <- 0 }
      }
      if (any(colnames(tPred) %in% "NN")){
        if(predict(model_nn,testSet_DM[ii,]) == levels(testSetClasses)[1])
        {tPred[ii,which(colnames(tPred) %in% "NN")] <- 1 }
        else{tPred[ii,which(colnames(tPred) %in% "NN")] <- 0 }
      }
      if (any(colnames(tPred) %in% "PLS")){
        if(predict(model_pls,testSet_DM[ii,]) == levels(testSetClasses)[1])
        {tPred[ii,which(colnames(tPred) %in% "PLS")] <- 1 }
        else{tPred[ii,which(colnames(tPred) %in% "PLS")] <- 0 }
      }
      if (any(colnames(tPred) %in% "kNN")){

        if (kknn(formula_DM,
                 trainingSet_DM,
                 testSet_DM[ii,],
                 k=3)$CL[,3] == levels(testSetClasses)[1])
        {tPred[ii,which(colnames(tPred) %in% "kNN")] <- 1 }
        else{tPred[ii,which(colnames(tPred) %in% "kNN")] <- 0 }
      }

      # ensemble Learning
      tPred[ii,
            which(
              colnames(tPred) %in% "Ensemble")] <- round(sum(
                wei*tPred[ii,-1]))
    }

   # calculate scores
      acc.Class <- matrix(nrow=1, ncol = ncol(tPred))
      MCC.Class <- acc.Class
      TP.Class <- acc.Class
      TN.Class <- acc.Class
      FP.Class <- acc.Class
      FN.Class <- acc.Class
      Sensit.Class <- acc.Class
      Specif.Class <- acc.Class
      PPV.Class <- acc.Class
      NPV.Class <- acc.Class

    # per MCC
    TP_Class <- colSums(tPred[
      which(testSetClasses == levels(testSetClasses)[1]),
      ,drop=FALSE] == 1)
    TN_Class <- colSums(tPred[
      which(testSetClasses == levels(testSetClasses)[2]),
      ,drop=FALSE] == 0)
    FP_Class <- colSums(tPred[
      which(testSetClasses == levels(testSetClasses)[1]),
      ,drop=FALSE] == 0)
    FN_Class <- colSums(tPred[
      which(testSetClasses == levels(testSetClasses)[2]),
      ,drop=FALSE] == 1)

    TP.Class <- TP_Class
    TN.Class <- TN_Class
    FP.Class <- FP_Class
    FN.Class <- FN_Class

    acc.Class <- 100 * (TP_Class + TN_Class)/(
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

  acc.Class <- as.data.frame(t(acc.Class))
  MCC.Class <- as.data.frame(t(MCC.Class))
  Sensit.Class <- as.data.frame(t(Sensit.Class))
  Specif.Class <- as.data.frame(t(Specif.Class))
  PPV.Class <- as.data.frame(t(PPV.Class))
  NPV.Class <- as.data.frame(t(NPV.Class))


  wei_all <- as.data.frame(t(wei_all))

  ############################ create output_data
  #output:list
  #wei: vector
  output_data_list <- list()
  label_list <- c()
  k_list <- 1
  if (any(colnames(tPred) %in% "RF")){

    #idx_model <- which(colnames(model_id_ok) %in% "RF")
    output_data_list[[k_list]] <- model_rf
    label_list <- c(label_list,"RF")

    k_list <- k_list + 1}

  if (any(colnames(tPred) %in% "SVM")){
    #idx_model <- which(colnames(model_id_ok) %in% "SVM")
    output_data_list[[k_list]] <- model_svm
    label_list <- c(label_list,"SVM")

    k_list <- k_list + 1}

  if (any(colnames(tPred) %in% "NB")){
    #idx_model <- which(colnames(model_id_ok) %in% "NB")
    output_data_list[[k_list]] <- model_nb
    label_list <- c(label_list,"NB")

    k_list <- k_list + 1}

  if (any(colnames(tPred) %in% "LDA")){
    #idx_model <- which(colnames(model_id_ok) %in% "LDA")
    output_data_list[[k_list]] <- model_lda
    label_list <- c(label_list,"LDA")

    k_list <- k_list + 1}

  if (any(colnames(tPred) %in% "LR")){
    #idx_model <- which(colnames(model_id_ok) %in% "LR")
    output_data_list[[k_list]] <- model_lr
    label_list <- c(label_list,"LR")

    k_list <- k_list + 1}

  if (any(colnames(tPred) %in% "NN")){
    #idx_model <- which(colnames(model_id_ok) %in% "NN")
    output_data_list[[k_list]] <- model_nn
    label_list <- c(label_list,"NN")

    k_list <- k_list + 1}

  if (any(colnames(tPred) %in% "PLS")){
    #idx_model <- which(colnames(model_id_ok) %in% "PLS")
    output_data_list[[k_list]] <- model_pls
    label_list <- c(label_list,"PLS")

    k_list <- k_list + 1}


  # weights
  idx_wei <- which(colnames(tPred) %in% "Ensemble")
  output_data_list[[k_list]] <- wei_all
  label_list <- c(label_list,"Ensemble")

  names(output_data_list) <-label_list
  output_data_list$classes <- as.factor(class_level2)
  output_data_list$positiveClass <- levels(classes)[1]

  return(model_list = output_data_list)
}
