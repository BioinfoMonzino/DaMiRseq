#' @title Build Classifier using 'Staking' Ensemble Learning strategy.
#'
#' @description This function implements a 'Stacking' ensemble learning strategy.
#' Users can provide covariates (other than genomic features) which will be taken into account during
#' classification model building. A 'two-classes' classification task is addressed.
#'
#' @param data A transposed data frame of normalized expression data, i.e. transformed counts
#' by vst or rlog. A log2 transformed expression matrix is also accepted.
#' Rows and Cols should be, respectively, observations and features
#' @param classes A class vector with \code{nrow(data)} elements.
#'  Each element represents the class label for each observation. Only two different class labels are
#'  allowed
#' @param covariates An optional data frame containing other covariates
#' (but without 'class' column). Each column represents a different covariate to be considered in the model
#' @param fSample.tr Fraction of samples to be used as training set;
#'  default is 0.7
#' @param fSample.tr.w Fraction of samples of training set to be used
#' during weight estimation; default is 0.7
#' @param iter Number of iterations to assess classification accuracy;
#' default is 100
#'
#' @return A list containing:
#' \itemize{
#'   \item A matrix of accuracies of each classifier in each iteration.
#'   \item A matrix of weights used for each classifier in each iteration.
#'   \item A list of all models generated in each iteration.
#'   \item A violin plot of model accuracy obtained for each iteration.
#' }
#'
#' @details
#' To assess the robustness of a set of predictors, a specific 'Stacking' strategy
#' has been implemented. First, a training set (TR1) and a test set (TS1) are generated
#' by 'bootstrap' sampling. Then, sampling again from TR1 subset, another
#' pair of training (TR2) and test set (TS2) are obtained. TR2 is used to train
#' Random Forest (RF), Naive Bayes (NB), Support Vector Machines
#' (SVM), k-Nearest Neighbour (kNN), Linear Discriminant Analysis (LDA) and Logistic
#' Regression (LR) classifiers, whereas TS2 is used to test their accuracy and to calculate weights.
#' The decision rule of 'Stacking' classifier is made by a linear combination of the
#' product between weigths (w) and predictions (Pr) of each classifier; for each sample k, the prediction
#' is computed by:
#'  \deqn{Pr_{k, Ensemble} = w_{RF} * Pr_{k, RF} + w_{NB} * Pr_{k, NB} + w_{SVM} * Pr_{k, SVM} + w_{k, kNN} * Pr_{k, kNN} +
#'  w_{k, LDA} * Pr_{k, LDA} + w_{k, LR} * Pr_{k, LR}}
#' Performance of 'Stacking' classifier is evaluated by using TS1. This process is
#' repeated several times (default 100 times).
#'
#' @author Mattia Chiesa, Luca Piacentini
#'
#' @examples
#' data(DaMiR.Data.Ex.classif)
#' data(DaMiR.Data.Ex.dfClCov)
#' classes <- DaMiR.Data.Ex.dfClCov$class
#' set.seed(1)
#' list.res <- DaMiR.EnsembleLearning(DaMiR.Data.Ex.classif, classes)
#'
#' @export
#'
#'
DaMiR.EnsembleLearning <- function(data, classes, covariates, fSample.tr=0.7, fSample.tr.w=0.7, iter=100){
  # check arguments
  if (missing(data)) stop("'data' argument must be provided")
  if (missing(classes)) stop("'classes' argument must be provided")

  if(!(is.data.frame(data))) stop("'data' must be a data frame")
  if(!(is.numeric(fSample.tr))) stop("'fSample.tr' must be numeric")
  if(!(is.numeric(fSample.tr.w))) stop("'fSample.tr.w' must be numeric")
  if(!(is.numeric(iter))) stop("'iter' must be numeric")

  # suppress warnings
  options( warn = -1 )

  # Initialize output tables
  accuracy_class <- matrix(nrow = iter, ncol = 7)
  weight_class <- matrix(nrow = iter, ncol = 6)
  model_class <- list()

  sample_cl1 <-length(classes[classes == levels(classes)[1]])
  sample_cl2 <-length(classes[classes == levels(classes)[2]])
  num_sample_pergroup_cl1 <- round(sample_cl1*fSample.tr)
  num_sample_pergroup_cl2 <- round(sample_cl2*fSample.tr)

  if (missing(covariates)){
    data <- data
  } else {
    if(!(is.data.frame(covariates))) stop("'covariates' must be a data frame") ###
    data <- cbind(data,covariates)
  }

  patt_cl1<-paste(c("^",levels(classes)[1],"$"),sep="",collapse="")
  patt_cl2<-paste(c("^",levels(classes)[2],"$"),sep="",collapse="")
  varNames <- colnames(data)
  varNames1 <- paste(varNames, collapse = "+")
  formula_DM <- as.formula(paste("classes", varNames1, sep = " ~ "))

  cat("Ensemble classification is running. ",iter," bootstrap iterations were chosen:","\n")
  # start bootstrapping
  for (jj in 1:iter){
    #cat(jj,"\n")

    #Random sampling
    sample_index_cl1_DM <- sample(length(grep(patt_cl1,classes)),replace = FALSE)
    sample_index_cl2_DM <- sample(length(grep(patt_cl2,classes)),replace = FALSE)+length(grep(patt_cl1,classes))

    # TrainingSet and TestSet
    classes_ts <- classes[c(sample_index_cl1_DM[-(1:num_sample_pergroup_cl1)],sample_index_cl2_DM[-(1:num_sample_pergroup_cl2)])]
    testSet <- data[c(sample_index_cl1_DM[-(1:num_sample_pergroup_cl1)],sample_index_cl2_DM[-(1:num_sample_pergroup_cl2)]),]
    testSet$classes <- classes_ts
    testSet_lr <- testSet
    testSet_lr$classes <- as.integer(classes_ts)-1

    classes_tr <- classes[c(sample_index_cl1_DM[1:num_sample_pergroup_cl1],sample_index_cl2_DM[1:num_sample_pergroup_cl2])]
    trainingSet <- data[c(sample_index_cl1_DM[1:num_sample_pergroup_cl1],sample_index_cl2_DM[1:num_sample_pergroup_cl2]),]
    trainingSet$classes <- classes_tr
    trainingSet_lr <- trainingSet
    trainingSet_lr$classes <- as.integer(classes_tr)-1

    # Random Sampling x Weight assignment
    sample_index_cl1_DM2 <- sample(length(grep(patt_cl1,classes_tr)),replace = FALSE)
    sample_index_cl2_DM2 <- sample(length(grep(patt_cl2,classes_tr)),replace = FALSE)+length(grep(patt_cl1,classes_tr))
    sample_cl1 <-length(classes_tr[classes_tr == levels(classes_tr)[1]])
    sample_cl2 <-length(classes_tr[classes_tr == levels(classes_tr)[2]])
    num_sample_pergroup2_cl1 <- round(sample_cl1*fSample.tr.w)
    num_sample_pergroup2_cl2 <- round(sample_cl2*fSample.tr.w)

    classes_ts_DM <- classes_tr[c(sample_index_cl1_DM2[-(1:num_sample_pergroup2_cl1)],sample_index_cl2_DM2[-(1:num_sample_pergroup2_cl2)])]
    testSet_DM <- trainingSet[c(sample_index_cl1_DM2[-(1:num_sample_pergroup2_cl1)],sample_index_cl2_DM2[-(1:num_sample_pergroup2_cl2)]),]
    testSet_lr_DM <- trainingSet_lr[c(sample_index_cl1_DM2[-(1:num_sample_pergroup2_cl1)],sample_index_cl2_DM2[-(1:num_sample_pergroup2_cl2)]),]

    classes_tr_DM <- classes_tr[c(sample_index_cl1_DM2[1:num_sample_pergroup2_cl1],sample_index_cl2_DM2[1:num_sample_pergroup2_cl2])]
    trainingSet_DM <- trainingSet[c(sample_index_cl1_DM2[1:num_sample_pergroup2_cl1],sample_index_cl2_DM2[1:num_sample_pergroup2_cl2]),]
    trainingSet_lr_DM <- trainingSet_lr[c(sample_index_cl1_DM2[1:num_sample_pergroup2_cl1],sample_index_cl2_DM2[1:num_sample_pergroup2_cl2]),]

    ###############################################
    ## Training several models

    # Random Forest
    model_rf <- randomForest(formula = formula_DM, data = trainingSet_DM, ntree=1000, importance =TRUE)

    # Support Vector Machine
    tune.for.model <- tune.svm(formula_DM,data = trainingSet_DM,gamma = 10^(-6:-1),cost = 10^(1:4))
    model_svm <- svm(formula = formula_DM,data = trainingSet_DM, gamma = tune.for.model$best.parameters$gamma, cost = tune.for.model$best.parameters$cost,hidden=3)

    # Naive Bayes
    model_nb <- naiveBayes(formula = formula_DM,data = trainingSet_DM)

    # Logistic Regression
    invisible(capture.output(model_lr <- glm(formula = formula_DM, data = trainingSet_lr_DM, family = binomial("logit"), maxit = 100)))

    # Linear Discriminant Analysis
    model_lda <- lda(formula = formula_DM, data = trainingSet_DM)

    list_model<-list(model_rf,model_svm, model_nb,model_lr,model_lda)
    ##################################################################################
    ## Calculate weight for each classifier based on accuracy in training prediction
    acc_model_rf <- confusionMatrix(table(predict(model_rf,testSet_DM),testSet_DM$classes), reference = testSet_DM$classes)$overall['Accuracy']
    acc_model_svm <- confusionMatrix(table(predict(model_svm,testSet_DM),testSet_DM$classes), reference = testSet_DM$classes)$overall['Accuracy']
    acc_model_nb <- confusionMatrix(table(predict(model_nb,testSet_DM),testSet_DM$classes), reference = testSet_DM$classes)$overall['Accuracy']
    acc_model_lr <- confusionMatrix(table(round(predict(model_lr,testSet_lr_DM,type='response')),testSet_lr_DM$classes))$overall['Accuracy']
    acc_model_lda <- confusionMatrix(table(predict(model_lda,testSet_DM)$class,testSet_DM$classes), reference = testSet_DM$classes)$overall['Accuracy']
    acc_model_3kNN <- confusionMatrix(table(kknn(formula_DM,trainingSet_DM,testSet_DM,k=3)$CL[,3],testSet_DM$classes), reference = testSet_DM$classes)$overall['Accuracy']


    names(acc_model_rf) <- 'RandomForest'
    names(acc_model_svm) <- 'SupportVectorMachine'
    names(acc_model_nb) <- 'NaiveBayes'
    names(acc_model_lr) <- 'LogisticRegression'
    names(acc_model_lda) <- 'LDA'
    names(acc_model_3kNN) <- '3kNN'
    ###############################################
    ## Unbalance weighting
    sort_acc <- sort(c(acc_model_rf,acc_model_svm,acc_model_nb,acc_model_lr,acc_model_lda,acc_model_3kNN),decreasing = TRUE)
    acc_all <- as.data.frame(sort_acc)
    #acc_all <- as.data.frame(sort_acc * rank(sort_acc))
    names(acc_all) <- 'unbalance weigh'
    acc_model_rf <- acc_all[rownames(acc_all) %in% 'RandomForest',]
    acc_model_svm <- acc_all[rownames(acc_all) %in% 'SupportVectorMachine',]
    acc_model_nb <- acc_all[rownames(acc_all) %in% 'NaiveBayes',]
    acc_model_lr <- acc_all[rownames(acc_all) %in% 'LogisticRegression',]
    acc_model_lda <- acc_all[rownames(acc_all) %in% 'LDA',]
    acc_model_3kNN <- acc_all[rownames(acc_all) %in% '3kNN',]

    #################################################
    ## create classification weight
    weight_rf <- as.numeric(acc_model_rf/(acc_model_3kNN+acc_model_rf+acc_model_svm+acc_model_nb+acc_model_lr+acc_model_lda))
    weight_svm <- as.numeric(acc_model_svm/(acc_model_3kNN+acc_model_rf+acc_model_svm+acc_model_nb+acc_model_lr+acc_model_lda))
    weight_nb <- as.numeric(acc_model_nb/(acc_model_3kNN+acc_model_rf+acc_model_svm+acc_model_nb+acc_model_lr+acc_model_lda))
    weight_lr <- as.numeric(acc_model_lr/(acc_model_3kNN+acc_model_rf+acc_model_svm+acc_model_nb+acc_model_lr+acc_model_lda))
    weight_lda <- as.numeric(acc_model_lda/(acc_model_3kNN+acc_model_rf+acc_model_svm+acc_model_nb+acc_model_lr+acc_model_lda))
    weight_3kNN <- as.numeric(acc_model_3kNN/(acc_model_3kNN+acc_model_rf+acc_model_svm+acc_model_nb+acc_model_lr+acc_model_lda))


    ###############################################
    ## Test New samples
    table_prediction_ensemble <- matrix(nrow = length(testSet[,1]), ncol = 2)
    rownames(table_prediction_ensemble) <- rownames(testSet)
    colnames(table_prediction_ensemble) <- c('EnsClass_predicted','sumOfPrediction')

    table_prediction <- matrix(nrow = length(testSet[,1]), ncol = 6)
    rownames(table_prediction) <- rownames(testSet)
    colnames(table_prediction) <- c('RF','SVM','NB','LDA','LR','3kNN')


    cl_comparison <-as.numeric(classes_ts)-1
    for (i in 1:length(testSet[,1])){
      # Prediction

      if (predict(model_rf,testSet[i,]) == levels(classes_ts)[1]){pr_mod_rf <- 0} else {pr_mod_rf <- 1}
      if (predict(model_svm,testSet[i,]) == levels(classes_ts)[1]){pr_mod_svm <- 0} else {pr_mod_svm <- 1}
      if (predict(model_nb,testSet[i,]) == levels(classes_ts)[1]){pr_mod_nb <- 0} else {pr_mod_nb <- 1}
      if (predict(model_lda,testSet[i,])$class == levels(classes_ts)[1]){pr_mod_lda <- 0} else {pr_mod_lda <- 1}
      if (kknn(formula_DM,trainingSet_DM,testSet[i,],k=3)$CL[,3] == levels(classes_ts)[1]){pr_mod_3kNN <- 0} else {pr_mod_3kNN <- 1}
      pr_mod_lr <- as.numeric(round(predict(model_lr,testSet_lr[i,],type='response')))

      # Decision
      table_prediction[i,1] <- pr_mod_rf
      table_prediction[i,2] <- pr_mod_svm
      table_prediction[i,3] <- pr_mod_nb
      table_prediction[i,4] <- pr_mod_lda
      table_prediction[i,5] <- pr_mod_lr
      table_prediction[i,6] <- pr_mod_3kNN



      table_prediction_ensemble[i,1] <- round(weight_3kNN * pr_mod_3kNN + weight_rf * pr_mod_rf + weight_svm * pr_mod_svm + weight_nb * pr_mod_nb + weight_lr * pr_mod_lr + weight_lda * pr_mod_lda)
      if (table_prediction_ensemble[i,1] == 1){
        table_prediction_ensemble[i,2] <- weight_3kNN * pr_mod_3kNN + weight_rf * pr_mod_rf + weight_svm * pr_mod_svm + weight_nb * pr_mod_nb + weight_lr * pr_mod_lr + weight_lda * pr_mod_lda
      }else{
        table_prediction_ensemble[i,2] <- 1-(weight_3kNN * pr_mod_3kNN + weight_rf * pr_mod_rf + weight_svm * pr_mod_svm + weight_nb * pr_mod_nb + weight_lr * pr_mod_lr + weight_lda * pr_mod_lda)
      }
    }
    accuracy_class[jj,1] <-(1-sum(abs(table_prediction_ensemble[,1]-cl_comparison))/length(cl_comparison))*100
    accuracy_class[jj,2] <-(1-sum(abs(table_prediction[,1]-cl_comparison))/length(cl_comparison))*100
    accuracy_class[jj,3] <-(1-sum(abs(table_prediction[,2]-cl_comparison))/length(cl_comparison))*100
    accuracy_class[jj,4] <-(1-sum(abs(table_prediction[,3]-cl_comparison))/length(cl_comparison))*100
    accuracy_class[jj,5] <-(1-sum(abs(table_prediction[,4]-cl_comparison))/length(cl_comparison))*100
    accuracy_class[jj,6] <-(1-sum(abs(table_prediction[,5]-cl_comparison))/length(cl_comparison))*100
    accuracy_class[jj,7] <-(1-sum(abs(table_prediction[,6]-cl_comparison))/length(cl_comparison))*100

    weight_class[jj,1] <- weight_rf
    weight_class[jj,2] <- weight_svm
    weight_class[jj,3] <- weight_nb
    weight_class[jj,4] <- weight_lr
    weight_class[jj,5] <- weight_lda
    weight_class[jj,6] <- weight_3kNN

    model_class[[jj]] <- list_model
    rm("list_model")
  }

  accuracy_class<-round(accuracy_class,2)
  colnames(accuracy_class) <- c('Ensemble','RF','SVM','NB','LDA','LR','3kNN')

  acc_dotplot <- melt(as.data.frame(accuracy_class), measure.vars = c('Ensemble','RF','SVM','NB','LDA','LR','3kNN'))
  colnames(acc_dotplot) <- c("Predictors","Accuracy")
  print(ggplot(acc_dotplot, aes(x=Predictors,y=Accuracy)) +
    geom_violin(aes(fill=factor(Predictors)),na.rm = TRUE)+
    geom_dotplot(binaxis='y', stackdir='center',stackratio=1.5, dotsize=0.2,binwidth = 0.5) +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color="white"))

  cat("Accuracy:", "\n","Ensemble RF SVM NB LDA LR 3kNN","\n","Mean:",round(colMeans(accuracy_class),2),"\n","St.Dev.", round(sd(accuracy_class[,1]),2), round(sd(accuracy_class[,2]),2), round(sd(accuracy_class[,3]),2), round(sd(accuracy_class[,4]),2), round(sd(accuracy_class[,5]),2), round(sd(accuracy_class[,6]),2), round(sd(accuracy_class[,7]),2))

   return(list(accuracy = accuracy_class, weight = weight_class, models = model_class))
}
