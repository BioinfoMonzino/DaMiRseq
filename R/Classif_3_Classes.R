#' @title Build a Multi-Class Classifier using 'Staking' Learning strategy.
#'
#' @description This function implements a 'Stacking' ensemble learning
#' strategy.
#' Users can provide heterogeneous features (other than genomic features)
#' which will be taken into account during
#' classification model building. A 'multi-classes' classification task is
#'  addressed.
#'
#' @param data A transposed data frame of normalized expression data.
#' Rows and Cols should be, respectively, observations and features
#' @param classes A class vector with \code{nrow(data)} elements.
#'  Each element represents the class label for each observation. Only
#'  more than two different class labels are allowed
#' @param variables An optional data frame containing other variables
#' (but without 'class' column). Each column represents a different
#' covariate to be considered in the model
#' @param fSample.tr Fraction of samples to be used as training set;
#'  default is 0.7
#' @param fSample.tr.w Fraction of samples of training set to be used
#' during weight estimation; default is 0.7
#' @param iter Number of iterations to assess classification accuracy;
#' default is 100
#' @param cl_type List of weak classifiers that will compose the
#' meta-learners. Only "RF", "kNN", "SVM", "LDA", "LR", "NB", "NN", "PLS"
#' are allowed. Default is c("RF", "LR", "kNN", "LDA", "NB", "SVM")
#'
#' @return A matrix of accuracies of each classifier in each iteration.
#'
#' @details
#' To assess the robustness of a set of predictors, a specific 'Stacking'
#' strategy
#' has been implemented. First, a training set (TR1) and a test set (TS1)
#'  are generated
#' by 'bootstrap' sampling. Then, sampling again from TR1 subset, another
#' pair of training (TR2) and test set (TS2) are obtained. TR2 is used to
#' train
#' Random Forest (RF), Naive Bayes (NB), Support Vector Machines
#' (SVM), k-Nearest Neighbour (kNN), Linear Discriminant Analysis (LDA)
#' and Logistic
#' Regression (LR) classifiers, whereas TS2 is used to test their accuracy
#'  and to calculate weights.
#' The decision rule of 'Stacking' classifier is made by a linear
#' combination of the
#' product between weigths (w) and predictions (Pr) of each classifier;
#' for each sample k, the prediction
#' is computed by:
#'  \deqn{Pr_{k, Ensemble} = w_{RF} * Pr_{k, RF} + w_{NB} * Pr_{k, NB} +
#'   w_{SVM} * Pr_{k, SVM} + w_{k, kNN} * Pr_{k, kNN} +
#'  w_{k, LDA} * Pr_{k, LDA} + w_{k, LR} * Pr_{k, LR}}
#' Performance of 'Stacking' classifier is evaluated by using TS1. This
#' process is
#' repeated several times (default 100 times).
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
#' # Classification_res <- DaMiR.EnsembleLearning(selected_features,
#' # classes=df$class, fSample.tr=0.6, fSample.tr.w=0.6, iter=3,
#' # cl_type=c("RF","kNN"))
#'
#' @export
#'
#'
DaMiR.EnsembleLearningNcl <- function(data,
                                   classes,
                                   variables,
                                   fSample.tr=0.7,
                                   fSample.tr.w=0.7,
                                   iter=100,
                                   cl_type=c("RF",
                                             "kNN",
                                             "SVM",
                                             "LDA",
                                             "LR",
                                             "NB",
                                             "NN",
                                             "PLS")){
  # check missing arguments
  if (missing(data))
    stop("'data' argument must be provided")
  if (missing(classes))
    stop("'classes' argument must be provided")
  if (missing(cl_type)){
    cl_type <- c("RF", "LR", "kNN", "LDA", "NB", "SVM")
  }

  # check the type of argument
  if(!(is.data.frame(data)))
    stop("'data' must be a data frame")
  if(!(is.numeric(fSample.tr)))
    stop("'fSample.tr' must be numeric")
  if(!(is.numeric(fSample.tr.w)))
    stop("'fSample.tr.w' must be numeric")
  if(!(is.numeric(iter)))
    stop("'iter' must be numeric")
  if(!(is.factor(classes)))
    stop("'classes' must be a factor")

  # specific checks
  if (fSample.tr >0.9 | fSample.tr < 0.5)
    stop("'fSample.tr' must be between 0.5 and 1")
  if (fSample.tr.w >0.9 | fSample.tr.w < 0.5)
    stop("'th.corr' must be between 0.5 and 1")
  if (iter < 1)
    stop("'iter' must be greater than 1")
  if((dim(data)[1]-round(dim(data)[1]*fSample.tr)) == 0)
    stop("The Test Set is not available. Decrease 'fSample.tr'
         or increase the number of observation.")
  if((dim(data)[1]-round(dim(data)[1]*fSample.tr.w)) == 0)
    stop("A Test Set is not available to weight classifiers.
         Decrease 'fSample.tr.w' or increase the number of observation.")
  if(length(classes) != dim(data)[1])
    stop("length(classes) must be equal to dim(data)[1]")

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

  # find the min number of sample per class
  min_sub <-0
  class_level <- levels(classes)
  for (i in seq_len(length(class_level))){
    min_sub[i] <- length(which(class_level[i]==classes))
  }
  min_sample_each_class <- min(min_sub)
  if (min_sample_each_class < 3 )
    stop("At least 3 samples are needed!")

  if (!(all(cl_type %in% c("RF", "kNN","SVM","LDA","LR","NB","NN","PLS"))))
    stop("'cl_type' must be
         'RF', 'kNN', 'SVM', 'LDA', 'LR', 'NB', 'NN', 'PLS'")
  ###################################
  acc.Class<-matrix()
  # find an exact number of sample for each group
  # use the others for the Independent TestSet
  sample_index_list <- ""
  for (i in seq_len(length(class_level))){
    sample_index <- sample(which(class_level[i]==classes),replace = FALSE)
    sample_index <- sample_index[seq_len(min_sample_each_class)]
    sample_index_list <- c(sample_index_list,sample_index)
  }
  sample_index_list <- as.numeric(sample_index_list[-1])

  # Start to create independent TestSet
  index_indep <- setdiff(seq_len(dim(data)[1]), sample_index_list)
  if(length(index_indep) != 0){
  independentTestSet <- data[index_indep,, drop=FALSE]
  independentClasses <- classes[index_indep, drop=FALSE]
  independentClasses <- droplevels(independentClasses)
  }
  # create datset for DM
  dataset <- data[sample_index_list,, drop=FALSE]
  datasetClasses <- classes[sample_index_list, drop=FALSE]
  datasetClasses <- droplevels(datasetClasses)

  # training and test
  cat("You select:",cl_type, "weak classifiers for creating
      the Ensemble meta-learner.","\n")
  cat("Ensemble classification is running. ",
      iter," iterations were chosen:","\n")


  acc.Class<- matrix()
  # start bootstrapping
  for (jj in seq_len(iter)){

    # Sampling
    sample_index_train <- ""
    class_level <- levels(datasetClasses)
    tr_sample <- round(dim(dataset)[1]*fSample.tr/length(class_level))
    for (i in seq_len(length(class_level))){
      sample_index <- sample(which(class_level[i]==datasetClasses),
                             replace = FALSE)[seq_len(tr_sample)]
      sample_index_train <- c(sample_index_train,sample_index)
    }
    sample_index_train <- as.numeric(sample_index_train[-1])

    # create TrainingSet for DM
    trainingSet <- dataset[sample_index_train,, drop=FALSE]
    trainingSetClasses <- datasetClasses[sample_index_train,
                                         drop=FALSE]
    trainingSetClasses <- droplevels(trainingSetClasses)

    # create TestSet
    sample_index_test <- setdiff(seq_len(dim(dataset)[1]),
                                 sample_index_train)
    testSet <- dataset[sample_index_test,, drop=FALSE]
    testSetClasses <- datasetClasses[sample_index_test]
    testSetClasses <- droplevels(testSetClasses)

    # merge TestSet and Independent testSet (if exists)
    if(length(index_indep) != 0){
      testSet <- rbind(independentTestSet,testSet)
      testSetClasses <- as.factor(c(as.character(independentClasses),
                                    as.character(testSetClasses)))
      testSetClasses <- droplevels(testSetClasses)
    }

    # Sample Training for weighting
    sample_index_train2 <- ""
    class_level2 <- levels(trainingSetClasses)
    tr_sample2 <- round(
      dim(trainingSet)[1]*fSample.tr.w/length(class_level2))
    for (i in seq_len(length(class_level2))){
      sample_index2 <- sample(which(class_level2[i]==trainingSetClasses),
                             replace = FALSE)[seq_len(tr_sample2)]
      sample_index_train2 <- c(sample_index_train2,sample_index2)
    }
    sample_index_train2 <- as.numeric(sample_index_train2[-1])

    # create TrainingSet for weighting
    trainingSet2 <- trainingSet[sample_index_train2,, drop=FALSE]
    trainingSetClasses2 <- trainingSetClasses[sample_index_train2,
                                              drop=FALSE]
    trainingSetClasses2 <- droplevels(trainingSetClasses2)

    # create TestSet for weighting
    sample_index_test2 <- setdiff(seq_len(dim(trainingSet)[1]),
                                  sample_index_train2)
    testSet2 <- trainingSet[sample_index_test2,, drop=FALSE]
    testSetClasses2 <- trainingSetClasses[sample_index_test2,
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
        predict(model_rf, testSet_DM),
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
      kk <- kk +1
    }

    # k-Nearest Neighbours
    if (any(cl_type %in% "kNN")){

      acc_model[kk] <- caret::confusionMatrix(
        as.factor(kknn(formula_DM, trainingSet_DM, testSet_DM, k=3)$CL[,3]),
        reference = testSet_DM$classes)$overall['Accuracy']
      names(acc_model)[kk] <- "kNN"

    }

    # Weighting
    wei <- acc_model/sum(acc_model)

    ###############################################
    ## Test New samples

    tPred <- matrix(nrow = dim(testSet)[1],
                    ncol = (length(wei) + 1))
    colnames(tPred) <- c("Ensemble",names(wei))

    for (ii in seq_len(dim(testSet)[1])){

      if (any(colnames(tPred) %in% "RF")){
        if(predict(model_rf,testSet[ii,]) == (testSetClasses)[ii])
        {tPred[ii,which(colnames(tPred) %in% "RF")] <- 1 }
        else{tPred[ii,which(colnames(tPred) %in% "RF")] <- 0 }
      }
      if (any(colnames(tPred) %in% "SVM")){
        if(predict(model_svm,testSet[ii,]) == (testSetClasses)[ii])
        {tPred[ii,which(colnames(tPred) %in% "SVM")] <- 1 }
        else{tPred[ii,which(colnames(tPred) %in% "SVM")] <- 0 }
      }
      if (any(colnames(tPred) %in% "NB")){
        if(predict(model_nb,testSet[ii,]) == (testSetClasses)[ii])
        {tPred[ii,which(colnames(tPred) %in% "NB")] <- 1 }
        else{tPred[ii,which(colnames(tPred) %in% "NB")] <- 0 }
      }
      if (any(colnames(tPred) %in% "LDA")){
        if(predict(model_lda,testSet[ii,])$class == (testSetClasses)[ii])
        {tPred[ii,which(colnames(tPred) %in% "LDA")] <- 1 }
        else{tPred[ii,which(colnames(tPred) %in% "LDA")] <- 0 }
      }
      if (any(colnames(tPred) %in% "LR")){
        if(predict(model_lr,testSet[ii,]) == (testSetClasses)[ii])
        {tPred[ii,which(colnames(tPred) %in% "LR")] <- 1 }
        else{tPred[ii,which(colnames(tPred) %in% "LR")] <- 0 }
      }
      if (any(colnames(tPred) %in% "NN")){
        if(predict(model_nn,testSet[ii,]) == (testSetClasses)[ii])
        {tPred[ii,which(colnames(tPred) %in% "NN")] <- 1 }
        else{tPred[ii,which(colnames(tPred) %in% "NN")] <- 0 }
      }
      if (any(colnames(tPred) %in% "PLS")){
        if(predict(model_pls,testSet[ii,]) == (testSetClasses)[ii])
        {tPred[ii,which(colnames(tPred) %in% "PLS")] <- 1 }
        else{tPred[ii,which(colnames(tPred) %in% "PLS")] <- 0 }
      }
      if (any(colnames(tPred) %in% "kNN")){

        if (kknn(formula_DM,
                 trainingSet_DM,
                 testSet[ii,],
                 k=3)$CL[,3] == (testSetClasses)[ii])
        {tPred[ii,which(colnames(tPred) %in% "kNN")] <- 1 }
        else{tPred[ii,which(colnames(tPred) %in% "kNN")] <- 0 }
      }

      # ensemble Learning
      tPred[ii,
            which(
              colnames(tPred) %in% "Ensemble")] <- round(sum(
                wei*tPred[ii,-1]))
    }

    ## ok till here
    if(jj==1){
      acc.Class <- matrix(nrow=iter, ncol = length(colMeans(tPred)))
    }
      acc.Class[jj,] <- colMeans(tPred) *100


  }
  colnames(acc.Class) <- colnames(tPred)
  acc.Class<-round(acc.Class,2)


  acc_dotplot <- melt(as.data.frame(acc.Class),
                      measure.vars = colnames(acc.Class))

  colnames(acc_dotplot) <- c("Classifiers","Accuracy")
  print(ggplot(acc_dotplot, aes(x=Classifiers,y=Accuracy)) +
    #ylim(min(acc_dotplot$Accuracy)-5,100) +
    geom_violin(aes(fill=factor(Classifiers)),na.rm = TRUE)+
    geom_dotplot(binaxis='y',
                 stackdir='center',
                 stackratio=1.5,
                 dotsize=0.2,
                 binwidth = 0.5) +
    stat_summary(fun.data=mean_sdl,
                 fun.args = list(mult=1),
                 geom="pointrange",
                 color="white"))

  cat("Accuracy:",
      "\n",
      colnames(acc.Class),
      "\n",
      "Mean:",round(colMeans(acc.Class),2),"\n","St.Dev.",
      round(colSds(acc.Class),digits = 1))

   return(list(accuracy = acc.Class))
}

