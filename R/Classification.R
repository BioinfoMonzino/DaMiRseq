#' @title Build Classifier using 'Staking' Ensemble Learning strategy.
#'
#' @description This function implements a 'Stacking' ensemble learning
#' strategy.
#' Users can provide heterogeneous features (other than genomic features)
#' which will be taken into account during
#' classification model building.
#'
#' @param data A transposed data frame of normalized expression data.
#' Rows and Cols should be, respectively, observations and features
#' @param classes A class vector with \code{nrow(data)} elements.
#'  Each element represents the class label for each observation.
#'  More than two different class labels are handled.
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
#' @return A list containing:
#' \itemize{
#'   \item A matrix of accuracies of each classifier in each iteration.
#'   \item A matrix of weights used for each classifier in each iteration.
#'   \item A list of all models generated in each iteration.
#'   \item A violin plot of model accuracy obtained for each iteration.
#' }
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
#'  \deqn{Pr_{k, Ensemble} = sum(w[RF] * Pr[k,i]), i = 1, N}
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
#'  classes=df$class, fSample.tr=0.6, fSample.tr.w=0.6, iter=3,
#'  cl_type=c("RF","kNN"))
#'
#' @export
#'
#'
DaMiR.EnsembleLearning <- function(data,
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

  ## body
  # check number of classes
  class_level <- levels(classes)

  if (length(class_level) < 2)
    stop("At least 2 classes must be provided")

  if (length(class_level) == 2){
    out2cl <- DaMiR.EnsembleLearning2cl(data=data,
                                        classes=classes,
                                        fSample.tr=fSample.tr,
                                        fSample.tr.w=fSample.tr.w,
                                        iter=iter,
                                        cl_type=cl_type)
    return(out2cl)
  }else{
    outNcl <- DaMiR.EnsembleLearningNcl(data=data,
                                        classes=classes,
                                        fSample.tr=fSample.tr,
                                        fSample.tr.w=fSample.tr.w,
                                        iter=iter,
                                        cl_type=cl_type)
    return(outNcl)
  }






}

