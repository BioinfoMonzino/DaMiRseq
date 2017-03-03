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
#' values with sv effects removed and the data frame of the variables.
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

  # check missing arguments
  if (missing(data))
    stop("'data' argument must be provided")
  if (missing(sv))
    stop("'sv' argument must be provided")
  if (missing(n.sv)){
    n.sv=dim(sv)[2]
    cat("All the sv have been used to adjust the data")
  }

  # check the type of argument
  if(!(is(data, "SummarizedExperiment")))
    stop("'data' must be a 'SummarizedExperiment' object")
  if(!(is.matrix(sv)))
    stop("'sv' must be a matrix")
  if(!(is.numeric(n.sv)))
    stop("'n.sv' must be numeric")

  # check the presence of NA or Inf
  if (any(is.na(assay(data))))
    stop("NA values are not allowed in the 'data' matrix")
  if (any(is.na(sv)))
    stop("NA values are not allowed in the 'sv' matrix")
  if (any(is.na(n.sv)))
    stop("NA values are not allowed in the 'n.sv' argument")
  if (any(is.infinite(assay(data))))
    stop("Inf values are not allowed in the 'data' matrix")
  if (any(is.infinite(sv)))
    stop("Inf values are not allowed in the 'sv' matrix")
  if (any(is.infinite(n.sv)))
    stop("Inf values are not allowed in the 'n.sv' argument")

  # specific checks
  if ((n.sv %% 1) != 0)
    stop("'n.sv' must be a number without decimals")
  if (n.sv > dim(sv)[2])
    stop ("'n.sv' must be <= to the number of variables in sv matrix")
  if (n.sv == 0)
    stop ("At least 1 sv must be provided")
  if (n.sv < 1)
    stop ("Negative values for 'n.sv' are not allowed")
  if(!("class" %in% colnames(colData(data))))
    stop("'class' info is lacking!
         Include the variable 'class'
         in colData(data) and label it 'class'!")
  if (all((assay(data) %%1) == 0))
    warning("It seems that you are using raw counts!
            This function works with normalized data")
  if(dim(assay(data))[2] != dim(sv)[1])
    stop("ncol(assay(data)) must be equal to nrow(sv)")



  surr_variables<-sv[, seq_len(n.sv)]

  data_adjust <- removeBatchEffect(assay(data),
                                   design = model.matrix(~data@colData$class),
                                   covariates = surr_variables)
  data_adjust <- SummarizedExperiment(assays=data_adjust,
                                    colData=as.data.frame(colData(data)))

  return(data_adjust)
}
