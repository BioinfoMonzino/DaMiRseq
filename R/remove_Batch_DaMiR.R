#' @title Remove variable effects from expression data
#'
#' @description This function removes surrogate or other confounding variable effects from
#' normalized expression data by the usage of \code{\link{removeBatchEffect}} function
#'  of \code{limma} package.
#'
#' @param data Matrix of normalized expression data, i.e. transformed counts by vst or rlog.
#' A log2 transformed expression matrix is also accepted
#' @param df A data frame with class and known covariates; at least one column with
#' 'class' label must be included
#' @param sv The matrix of surrogate variables identified by \code{\link{DaMiR.SV}} function
#' @param n.sv The number of surrogate variables to be used to adjust the data
#'
#' @return A list containing:
#' \itemize{
#'   \item A normalized and adjusted expression matrix.
#'   \item A data frame with class and optional covariates information.
#' }
#'
#' @author Mattia Chiesa, Luca Piacentini
#'
#' @examples
#' data(DaMiR.Data.Ex.Norm)
#' data(DaMiR.Data.Ex.dfClCov)
#' data(DaMiR.Data.Ex.SV)
#' data.adjusted <- DaMiR.SVadjust(DaMiR.Data.Ex.Norm, DaMiR.Data.Ex.dfClCov, sv = DaMiR.Data.Ex.SV)
#' data.adjusted <- DaMiR.SVadjust(DaMiR.Data.Ex.Norm, DaMiR.Data.Ex.dfClCov, sv = DaMiR.Data.Ex.SV, n.sv = 3)
#'
#' @seealso
#' \code{\link{removeBatchEffect}},
#' \code{\link{DaMiR.SV}}
#'
#' @export
#'
DaMiR.SVadjust <- function(data, df, sv, n.sv){

  # check arguments
  if (missing(data)) stop("'data' argument must be provided")
  if (missing(df)) stop("'df' argument must be provided")
  if (missing(sv)) stop("'sv' argument must be provided")
  if (missing(n.sv)){
    n.sv=dim(sv)[2]
    cov_correz<-sv[,1:n.sv]
    # cov_correz<-sv
  }else{
    if (n.sv> dim(sv)[2]) stop ("'n.sv' must be <= to the number of variables in sv matrix")
    cov_correz<-sv[,1:n.sv]
  }

  if(!(is.matrix(data))) stop("'data' must be a matrix")
  if(!(is.data.frame(df))) stop("'df' must be a data frame")
  if(!(is.matrix(sv))) stop("'sv' must be a matrix")
  if(!(is.numeric(n.sv))) stop("'n.sv' must be numeric")

  data_norm_cleaned <- removeBatchEffect(data, design = model.matrix(~df$class), covariates = cov_correz)
  colnames(data_norm_cleaned) <- colnames(data)

  return(list(data = data_norm_cleaned, classCov = df))
}
