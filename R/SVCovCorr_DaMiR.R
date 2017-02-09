#' @title Correlation Plot
#'
#' @description This function easily draws the correlation plot of surrogate variables (sv)
#' and covariates.
#'
#' @param sv The matrix of sv identified by \code{\link{DaMiR.SV}} function
#' @param df A data frame with class and known covariates; at least one column with
#' 'class' label must be included
#' @param type Type of correlation metric to be applied; default is "pearson"
#' @param sig.level The significance level of the correlation; default is 0.0001
#'
#' @details Factorial variables are allowed. They will be tranformed as numeric before
#' applying the \code{\link{rcorr}} function of \code{Hmisc}.The \code{\link{corrplot}}
#' function, which draws the plot, marks with a cross all the correlations that do not reach
#' the significance threshold defined in the \code{sig.level} argument.This plot allows the user to
#' identify those sv that present significant correlations with either technical and
#' biological known variables.
#' Notably, none of the sv should present signifcant correlation with "class" variable.
#'
#' @return
#' A correlation plot between sv and known covariates.
#'
#' @author Mattia Chiesa, Luca Piacentini
#'
#' @examples
#' data(DaMiR.Data.Ex.dfClCov)
#' data(DaMiR.Data.Ex.SV)
#' DaMiR.corrplot(DaMiR.Data.Ex.SV, DaMiR.Data.Ex.dfClCov)
#' DaMiR.corrplot(DaMiR.Data.Ex.SV, DaMiR.Data.Ex.dfClCov, type = "spearman", sig.level=0.01)
#'
#' @seealso
#' \code{\link{DaMiR.SV}}
#'
#' @export
#'
#'
DaMiR.corrplot <- function(sv, df, type=c("pearson","spearman"), sig.level=0.0001){
  if (missing(sv)) stop("'sv' argument must be provided")
  if (missing(df)) stop("'df' argument must be provided")
  if (missing(type)){
    type <- type[1]
  }
  if(!(is.matrix(sv))) stop("'sv' must be a matrix")
  if(!(is.data.frame(df))) stop("'df' must be a data frame")
  if(!(is.numeric(sig.level))) stop("'sig.level' must be numeric")

  sva_corr <- as.data.frame(cbind(sv,df))
  for (i in 1:ncol(sva_corr)){
    if (is.numeric(sva_corr[,i])==FALSE){
      droplevels(sva_corr[,i])
      sva_corr[,i] <- as.numeric(sva_corr[,i])
    }
  }

  if(type == "pearson"){
    corrplot::corrplot(rcorr(as.matrix(sva_corr), type='pearson')$r, type = "upper", is.corr = TRUE, addCoef.col = "black", number.cex=0.5, p.mat = rcorr(as.matrix(sva_corr), type='pearson')$P, sig.level = sig.level)
    }else if(type == "spearman"){
      corrplot::corrplot(rcorr(as.matrix(sva_corr), type='spearman')$r, type = "upper", is.corr = TRUE, addCoef.col = "black", number.cex=0.5, p.mat = rcorr(as.matrix(sva_corr), type='spearman')$P, sig.level = sig.level)
    } else{
      stop("Please set 'spearman or 'pearson' as correlation type.")
    }
}
