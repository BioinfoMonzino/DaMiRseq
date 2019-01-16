#' @title Train a Binary Classifier using 'Staking' Learning strategy.
#'
#' @description This function learn a meta learner by a 'Stacking'
#' strategy.
#' Users can provide heterogeneous features (other than genomic features)
#' which will be taken into account during
#' classification model building. A 'two-classes' classification task is
#'  addressed.
#'
#' @param data A transposed data frame of normalized expression data.
#' Rows and Cols should be, respectively, observations and features
#' @param classes A class vector with \code{nrow(data)} elements.
#'  Each element represents the class label for each observation. Two
#'  different class labels are allowed
#' @param variables An optional data frame containing other variables
#' (but without 'class' column). Each column represents a different
#' covariate to be considered in the model
#' @param fSample.tr.w Fraction of samples of training set to be used
#' during weight estimation; default is 0.7
#' @param iter Number of iterations to assess classification accuracy;
#' default is 100
#' @param cl_type List of weak classifiers that will compose the
#' meta-learners. "RF", "SVM", "LDA", "LR", "NB", "NN", "PLS"
#' are allowed. Default is c("RF", "LR", "LDA", "NB", "SVM")
#' @param type_model Select the model with the median or best accuracy
#' over the iteration. "median" and "best" are allowed.
#' Default: median
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
#' #  Classification_res <- DaMiR.EnsembleLearning2cl_Training(
#' #  selected_features,classes=df$class, fSample.tr.w=0.6, iter=3,
#' # cl_type=c("RF","LR"))
#'
#' @export
#'
#'
DaMiR.EnsembleLearning2cl_Training <- function(data,
                                      classes,
                                      variables,
                                      fSample.tr.w=0.7,
                                      iter=100,
                                      cl_type=c("RF",
                                                "SVM",
                                                "LDA",
                                                "LR",
                                                "NB",
                                                "NN",
                                                "PLS"),
                                      type_model=c("median",
                                                   "best")){
  # check missing arguments
  if (missing(data))
    stop("'data' argument must be provided")
  if (missing(classes))
    stop("'classes' argument must be provided")
  if (missing(cl_type)){
    cl_type <- c("RF", "LR", "LDA", "NB", "SVM")
  }
  if (missing(type_model)){
    type_model <- "median"
  }



  # check the type of argument
  if(!(is.data.frame(data)))
    stop("'data' must be a data frame")
  if(!(is.numeric(fSample.tr.w)))
    stop("'fSample.tr.w' must be numeric")
  if(!(is.numeric(iter)))
    stop("'iter' must be numeric")
  if(!(is.factor(classes)))
    stop("'classes' must be a factor")

  # specific checks
  if (fSample.tr.w >0.9 | fSample.tr.w < 0.5)
    stop("'th.corr' must be between 0.5 and 1")
  if (iter < 1)
    stop("'iter' must be greater than 1")
  if((dim(data)[1]-round(dim(data)[1]*fSample.tr.w)) == 0)
    stop("A Test Set is not available to weight classifiers.
         Decrease 'fSample.tr.w' or increase the number of observation.")
  if(length(classes) != dim(data)[1])
    stop("length(classes) must be equal to dim(data)[1]")
	
  # check balanced dataset (Thanks to Dr. Pawel Karpinski)
  class_lev_ckeck <- levels(classes)
  count_cl1 <- length(which(classes %in% class_lev_ckeck[1]))
  count_cl2 <- length(which(classes %in% class_lev_ckeck[2]))
  if(count_cl1 != count_cl2)
    stop("Trainingset must be balanced (same n. of samples in both classes)")

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

  if (!(all(cl_type %in% c("RF", "SVM","LDA","LR","NB","NN","PLS"))))
    stop("'cl_type' must be
         'RF', 'SVM', 'LDA', 'LR', 'NB', 'NN', 'PLS'")

  if (!(all(type_model %in% c("median","best"))))
    stop("'type_model' must be 'median','best'")

  # training and test
  cat("You select:",cl_type, "weak classifiers for creating
      the Ensemble meta-learner.","\n")
  cat("Ensemble classification is running. ",
      iter," iterations were chosen:","\n")


  acc.Class<- matrix()

  wei_all <- matrix(nrow = iter,ncol = length(cl_type))
  model_rf_all <- list()
  model_svm_all <- list()
  model_nb_all <- list()
  model_lr_all <- list()
  model_nn_all <- list()
  model_lda_all <- list()
  model_pls_all <- list()


  # start training
  for (jj in seq_len(iter)){

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
    formula_DM <- as.formula(paste("classes", varNames1, sep = " ~ "))

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
      model_rf_all[[jj]] <- model_rf
    }


    # Support Vector Machine
    if (any(cl_type %in% "SVM")){

      tune.for.model <- tune.svm(formula_DM,
                                 data = trainingSet_DM,
                                 gamma = 10^(-6:-1),
                                 cost = 10^(1:4),
                                 tunecontrol = tune.control(
                                   # cross = min(min_sub,10) ))
                                   cross = min(dim(trainingSet_DM)[1],10) ))
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
      model_svm_all[[jj]] <- model_svm

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
      model_nb_all[[jj]] <- model_nb
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
      model_lda_all[[jj]] <- model_lda
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
      model_lr_all[[jj]] <- model_lr
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
      model_nn_all[[jj]] <- model_nn

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
      model_pls_all[[jj]] <- model_pls

      #kk <- kk +1
    }

    # # k-Nearest Neighbours
    # if (any(cl_type %in% "kNN")){
    #
    #   acc_model[kk] <- caret::confusionMatrix(
    #     as.factor(kknn(formula_DM, trainingSet_DM, testSet_DM, k=3)$CL[,3]),
    #     reference = testSet_DM$classes)$overall['Accuracy']
    #   names(acc_model)[kk] <- "kNN"
    #
    # }

    # Weighting
    wei <- acc_model/sum(acc_model)
    wei_all[jj,] <- wei

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
    if(jj==1){
      acc.Class <- matrix(nrow=iter, ncol = length(colMeans(tPred)))
      MCC.Class <- acc.Class
      TP.Class <- acc.Class
      TN.Class <- acc.Class
      FP.Class <- acc.Class
      FN.Class <- acc.Class
      Sensit.Class <- acc.Class
      Specif.Class <- acc.Class
      PPV.Class <- acc.Class
      NPV.Class <- acc.Class

    }

    # per MCC
    TP_Class <- colSums(tPred[
      which(testSetClasses == levels(testSetClasses)[1]),,drop=FALSE] == 1)
    TN_Class <- colSums(tPred[
      which(testSetClasses == levels(testSetClasses)[2]),,drop=FALSE] == 0)
    FP_Class <- colSums(tPred[
      which(testSetClasses == levels(testSetClasses)[1]),,drop=FALSE] == 0)
    FN_Class <- colSums(tPred[
      which(testSetClasses == levels(testSetClasses)[2]),,drop=FALSE] == 1)

    TP.Class[jj,] <- TP_Class
    TN.Class[jj,] <- TN_Class
    FP.Class[jj,] <- FP_Class
    FN.Class[jj,] <- FN_Class

    acc.Class[jj,] <- 100 * (TP_Class + TN_Class)/(
      TP_Class + TN_Class +  FP_Class + FN_Class)

    MCC.Class[jj,] <- (TP_Class * TN_Class - FP_Class * FN_Class) /
      sqrt((TP_Class + FP_Class) * (TP_Class + FN_Class) *
             (TN_Class + FP_Class) * (TN_Class + FN_Class))

    Sensit.Class[jj,] <- TP_Class / (TP_Class + FN_Class)
    Specif.Class[jj,] <- TN_Class / (TN_Class + FP_Class)
    PPV.Class[jj,] <- TP_Class / (TP_Class + FP_Class)
    NPV.Class[jj,] <- TN_Class / (TN_Class + FN_Class)

  }

  acc.Class[which(is.nan(acc.Class))] <- 0
  MCC.Class[which(is.nan(MCC.Class))] <- 0
  Sensit.Class[which(is.nan(Sensit.Class))] <- 0
  Specif.Class[which(is.nan(Specif.Class))] <- 0
  PPV.Class[which(is.nan(PPV.Class))] <- 0
  NPV.Class[which(is.nan(NPV.Class))] <- 0

  colnames(acc.Class) <- colnames(tPred)
  colnames(MCC.Class) <- colnames(tPred)
  colnames(Sensit.Class) <- colnames(tPred)
  colnames(Specif.Class) <- colnames(tPred)
  colnames(PPV.Class) <- colnames(tPred)
  colnames(NPV.Class) <- colnames(tPred)

  rownames(acc.Class) <- paste0("M",seq_len(iter))
  rownames(MCC.Class) <- paste0("M",seq_len(iter))
  rownames(Sensit.Class) <- paste0("M",seq_len(iter))
  rownames(Specif.Class) <- paste0("M",seq_len(iter))
  rownames(PPV.Class) <- paste0("M",seq_len(iter))
  rownames(NPV.Class) <- paste0("M",seq_len(iter))

  ##########################find good model

  if (type_model == "median"){
    if((iter %% 2) == 0){
      #even
      id_models <- iter/2
      }else{
      #odd
      id_models <- iter/2 - 0.5
    }
  }else{

    id_models <- iter
  }
  model_id_ok <-c()
  for (iii in seq_len(dim(acc.Class)[2])){
    model_id_ok[iii] <- names(sort(acc.Class[,iii])[id_models])
  }
  colnames(wei_all) <- colnames(tPred[,-1])
  model_id_ok <- as.data.frame(t(as.numeric(sub("M","",model_id_ok))))
  colnames(model_id_ok) <- c("Ensemble",colnames(tPred[,-c(1),drop=FALSE]))

  ############################ create output_data
  #output:list
  #wei: vector
  output_data_list <- list()
  label_list <- c()
  k_list <- 1
  if (any(colnames(model_id_ok) %in% "RF")){

    #idx_model <- which(colnames(model_id_ok) %in% "RF")
    idx_model <- which(colnames(model_id_ok) %in% "Ensemble")
    output_data_list[[k_list]] <- model_rf_all[model_id_ok[,idx_model]]
    label_list <- c(label_list,"RF")

    k_list <- k_list + 1}

  if (any(colnames(model_id_ok) %in% "SVM")){
    #idx_model <- which(colnames(model_id_ok) %in% "SVM")
    idx_model <- which(colnames(model_id_ok) %in% "Ensemble")
    output_data_list[[k_list]] <- model_svm_all[model_id_ok[,idx_model]]
    label_list <- c(label_list,"SVM")

    k_list <- k_list + 1}

  if (any(colnames(model_id_ok) %in% "NB")){
    #idx_model <- which(colnames(model_id_ok) %in% "NB")
    idx_model <- which(colnames(model_id_ok) %in% "Ensemble")
    output_data_list[[k_list]] <- model_nb_all[model_id_ok[,idx_model]]
    label_list <- c(label_list,"NB")

    k_list <- k_list + 1}

  if (any(colnames(model_id_ok) %in% "LDA")){
    #idx_model <- which(colnames(model_id_ok) %in% "LDA")
    idx_model <- which(colnames(model_id_ok) %in% "Ensemble")
    output_data_list[[k_list]] <- model_lda_all[model_id_ok[,idx_model]]
    label_list <- c(label_list,"LDA")

    k_list <- k_list + 1}

  if (any(colnames(model_id_ok) %in% "LR")){
    #idx_model <- which(colnames(model_id_ok) %in% "LR")
    idx_model <- which(colnames(model_id_ok) %in% "Ensemble")
    output_data_list[[k_list]] <- model_lr_all[model_id_ok[,idx_model]]
    label_list <- c(label_list,"LR")

    k_list <- k_list + 1}

  if (any(colnames(model_id_ok) %in% "NN")){
    #idx_model <- which(colnames(model_id_ok) %in% "NN")
    idx_model <- which(colnames(model_id_ok) %in% "Ensemble")
    output_data_list[[k_list]] <- model_nn_all[model_id_ok[,idx_model]]
    label_list <- c(label_list,"NN")

    k_list <- k_list + 1}

  if (any(colnames(model_id_ok) %in% "PLS")){
    #idx_model <- which(colnames(model_id_ok) %in% "PLS")
    idx_model <- which(colnames(model_id_ok) %in% "Ensemble")
    output_data_list[[k_list]] <- model_pls_all[model_id_ok[,idx_model]]
    label_list <- c(label_list,"PLS")

    k_list <- k_list + 1}


  # weights
  idx_wei <- which(colnames(model_id_ok) %in% "Ensemble")
  output_data_list[[k_list]] <- wei_all[model_id_ok[,idx_wei],]
  label_list <- c(label_list,"Ensemble")

  names(output_data_list) <-label_list
  output_data_list$classes <- as.factor(class_level2)

  ## plots
  acc.Class<-round(acc.Class,2)
  MCC.Class<-round(MCC.Class,2)
  Sensit.Class<-round(Sensit.Class,2)
  Specif.Class<-round(Specif.Class,2)
  PPV.Class<-round(PPV.Class,2)
  NPV.Class<-round(NPV.Class,2)


  ## accuracy
  acc_dotplot <- melt(as.data.frame(acc.Class),
                      measure.vars = colnames(acc.Class))

  colnames(acc_dotplot) <- c("Classifiers","Accuracy")
  print(ggplot(acc_dotplot, aes(x=Classifiers,y=Accuracy)) +
          geom_violin(aes(fill=factor(Classifiers)),na.rm = TRUE)+
          # geom_dotplot(binaxis='y',
          #              stackdir='center',
          #              stackratio=1.5,
          #              dotsize=0.2,
          #              binwidth = 0.5) +
          stat_summary(fun.data=mean_sdl,
                       fun.args = list(mult=1),
                       geom="pointrange",
                       color="white") +
          coord_cartesian(ylim=c(min(acc.Class)-5,100))
  )
  ##
  ##
  mcc_dotplot <- melt(as.data.frame(MCC.Class),
                      measure.vars = colnames(MCC.Class))

  colnames(mcc_dotplot) <- c("Classifiers","MCC")
  print(ggplot(mcc_dotplot, aes(x=Classifiers,y=MCC)) +
          geom_violin(aes(fill=factor(Classifiers)),na.rm = TRUE)+
          # geom_dotplot(binaxis='y',
          #              stackdir='center',
          #              stackratio=1.5,
          #              dotsize=0.002,
          #              binwidth = 0.5) +
          stat_summary(fun.data=mean_sdl,
                       fun.args = list(mult=1),
                       geom="pointrange",
                       color="white") +
          coord_cartesian(ylim=c(min(MCC.Class)-0.05,1))
  )
  ##
  ##
  sen_dotplot <- melt(as.data.frame(Sensit.Class),
                      measure.vars = colnames(Sensit.Class))

  colnames(sen_dotplot) <- c("Classifiers","Sensitivity")
  print(ggplot(sen_dotplot, aes(x=Classifiers,y=Sensitivity)) +
          #ylim(0,1) +
          geom_violin(aes(fill=factor(Classifiers)),na.rm = TRUE)+
          # geom_dotplot(binaxis='y',
          #              stackdir='center',
          #              stackratio=1.5,
          #              dotsize=0.002,
          #              binwidth = 0.5) +
          stat_summary(fun.data=mean_sdl,
                       fun.args = list(mult=1),
                       geom="pointrange",
                       color="white") +
          coord_cartesian(ylim=c(min(Sensit.Class)-0.05,1))
  )
  ##
  ##
  ##
  spe_dotplot <- melt(as.data.frame(Specif.Class),
                      measure.vars = colnames(Specif.Class))

  colnames(spe_dotplot) <- c("Classifiers","Specificity")
  print(ggplot(spe_dotplot, aes(x=Classifiers,y=Specificity)) +
          #ylim(0,1) +
          geom_violin(aes(fill=factor(Classifiers)),na.rm = TRUE)+
          # geom_dotplot(binaxis='y',
          #              stackdir='center',
          #              stackratio=1.5,
          #              dotsize=0.002,
          #              binwidth = 0.5) +
          stat_summary(fun.data=mean_sdl,
                       fun.args = list(mult=1),
                       geom="pointrange",
                       color="white") +
          coord_cartesian(ylim=c(min(Specif.Class)-0.05,1))
  )
  ##
  ##
  ppv_dotplot <- melt(as.data.frame(PPV.Class),
                      measure.vars = colnames(PPV.Class))

  colnames(ppv_dotplot) <- c("Classifiers","PPV")
  print(ggplot(ppv_dotplot, aes(x=Classifiers,y=PPV)) +
          #ylim(0,1) +
          geom_violin(aes(fill=factor(Classifiers)),na.rm = TRUE)+
          # geom_dotplot(binaxis='y',
          #              stackdir='center',
          #              stackratio=1.5,
          #              dotsize=0.002,
          #              binwidth = 0.5) +
          stat_summary(fun.data=mean_sdl,
                       fun.args = list(mult=1),
                       geom="pointrange",
                       color="white") +
          coord_cartesian(ylim=c(min(PPV.Class)-0.05,1))
  )
  ##
  ##
  npv_dotplot <- melt(as.data.frame(NPV.Class),
                      measure.vars = colnames(NPV.Class))

  colnames(npv_dotplot) <- c("Classifiers","NPV")
  print(ggplot(npv_dotplot, aes(x=Classifiers,y=NPV)) +
          #ylim(0,1) +
          geom_violin(aes(fill=factor(Classifiers)),na.rm = TRUE)+
          # geom_dotplot(binaxis='y',
          #              stackdir='center',
          #              stackratio=1.5,
          #              dotsize=0.002,
          #              binwidth = 0.5) +
          stat_summary(fun.data=mean_sdl,
                       fun.args = list(mult=1),
                       geom="pointrange",
                       color="white") +
          coord_cartesian(ylim=c(min(NPV.Class)-0.05,1))
  )

  cat("Accuracy [%]:",
      "\n",
      colnames(acc.Class),
      "\n",
      "Mean:",round(colMeans(acc.Class),2),"\n","St.Dev.",
      round(colSds(acc.Class),digits = 2),"\n")
  cat("MCC score:",
      "\n",
      colnames(MCC.Class),
      "\n",
      "Mean:",round(colMeans(MCC.Class),2),"\n","St.Dev.",
      round(colSds(MCC.Class),digits = 2),"\n")
  cat("Sensitivity:",
      "\n",
      colnames(Sensit.Class),
      "\n",
      "Mean:",round(colMeans(Sensit.Class),2),"\n","St.Dev.",
      round(colSds(Sensit.Class),digits = 2),"\n")
  cat("Specificity:",
      "\n",
      colnames(Specif.Class),
      "\n",
      "Mean:",round(colMeans(Specif.Class),2),"\n","St.Dev.",
      round(colSds(Specif.Class),digits = 2),"\n")
   cat("PPV:",
       "\n",
       colnames(PPV.Class),
       "\n",
       "Mean:",round(colMeans(PPV.Class),2),"\n","St.Dev.",
       round(colSds(PPV.Class),digits = 2),"\n")
   cat("NPV:",
       "\n",
       colnames(NPV.Class),
       "\n",
       "Mean:",round(colMeans(NPV.Class),2),"\n","St.Dev.",
       round(colSds(NPV.Class),digits = 2),"\n")

  return(model_list = output_data_list)
}
