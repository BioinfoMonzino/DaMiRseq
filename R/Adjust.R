#' @title Remove variable effects from expression data
#'
#' @description This function removes surrogate or other confounding
#' variable effects from
#' normalized expression data by the usage of
#' \code{\link{removeBatchEffect}} function
#'  of \code{limma} package.
#'
#' @param data A SummarizedExpression object
#' @param sv The matrix of surrogate variables identified by
#' \code{\link{DaMiR.SV}} function
#' @param n.sv The number of surrogate variables to be used to adjust
#'  the data
#'
#' @return A SummarizedExpression object containing a matrix of
#' log-expression
#' values with sv effects removed and the data frame of the covariates.
#'
#' @author Mattia Chiesa, Luca Piacentini
#'
#' @examples
#' # use example data:
#' data(data_norm)
#' data(sv)
#' data_adjust <- DaMiR.SVadjust(data_norm, sv = sv, n.sv = 3)
#'
#' @seealso
#' \code{\link{removeBatchEffect}},
#' \code{\link{DaMiR.SV}}
#'
#' @export
#'
DaMiR.SVadjust <- function(data, sv, n.sv){

  # check arguments
  if (missing(data)) stop("'data' argument must be provided")
  # if (missing(df)) stop("'df' argument must be provided")
  if (missing(sv)) stop("'sv' argument must be provided")
  if (missing(n.sv)){
    n.sv=dim(sv)[2]
    surr_variables<-sv[, 1:n.sv]

  } else {
    if (n.sv> dim(sv)[2])
      stop ("'n.sv' must be <= to the number of variables in sv matrix")
    surr_variables<-sv[, 1:n.sv]
  }

  if(!(is(data, "SummarizedExperiment")))
    stop("'data' must be a 'SummarizedExperiment' object")
  if(!(is.matrix(sv))) stop("'sv' must be a matrix")
  if(!(is.numeric(n.sv))) stop("'n.sv' must be numeric")

  data_adjust <- removeBatchEffect(assay(data),
                                   design = model.matrix(~data@colData$class),
                                   covariates = surr_variables)
  data_adjust <- SummarizedExperiment(assays=data_adjust,
                                    colData=as.data.frame(colData(data)))

  return(data_adjust)
}
