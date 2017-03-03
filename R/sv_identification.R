#' @title Identification of Surrogate Variables
#'
#' @description This function returns a matrix of surrogate variables (sv)
#'  using the implementation
#'  by Chiesa-Piacentini or the sva method by Leek et al.
#'
#' @param data A SummarizedExpression object
#' @param method The method used to identify sv. If missing, the "fve"
#' method
#' will be selected. Otherwise the method "leek" or "be" should be choosen
#' @param th.fve This argument sets the threshold of maximum fraction of
#' variance
#' explained (fve) to be used in conjunction with "fve" method; default is
#'  0.95
#'
#' @details
#' This function helps the user to identify the appropriate number of sv:
#' it is possible to select a different strategy to be used by changing the
#'  option in \code{method}
#' argument. Three methods are available:
#'  \itemize{
#'   \item "be" - this option uses the \code{num.sv} function of \code{sva}
#'    package with default parameters;
#'   \item "leek" - The same of before but with asymptotic approach proposed
#'    by Leek;
#'   \item "fve" - This method is introduced in \code{DaMiRseq} package,
#'   and integrates part
#'    of \code{sva} function with custom code. Briefly, we computed
#'    eigenvalues of \code{data}
#'    using code already implemented in \code{sva} function and then, we
#'    calculated the squared
#'    of each eigenvalues. Thus, the ratio between each "squared eigenvalue"
#'     and the sum of them
#'    were calculated. These values represent a surrogate measure of the
#'    "Percentage of
#'    Explained Variance" (pve) obtained by principal component analysis (PCA),
#'     and their cumulative
#'    sum can be used to select sv.
#' }
#'
#' @return A matrix of sv.
#' A plot with the sv identified by "fve" method is also returned. A red dot shows
#'  the
#'  maximum number of variables to be included for a specific "fve".
#'
#' @references Jeffrey T. Leek, W. Evan Johnson, Hilary S. Parker, Elana J. Fertig,
#'  Andrew E.
#' Jaffe and John D. Storey (2016). sva: Surrogate Variable Analysis. R package
#' version 3.22.0.
#'
#' @author Mattia Chiesa, Luca Piacentini
#'
#' @examples
#' # use example data:
#' data(data_norm)
#' sv <- DaMiR.SV(data_norm, method = "fve", th.fve=0.95)
#'
#' @seealso
#' \code{\link{sva}}
#'
#' @export
#'
#'
DaMiR.SV <- function(data,
                     method=c("fve", "leek", "be"),
                     th.fve=0.95) {

  # check arguments
  if (missing(data))
    stop("'data' argument must be provided")

  if (missing(method)){
    method <- method[1]
  }

  # check the type of argument
  if(!(is(data, "SummarizedExperiment")))
    stop("'data' must be a 'SummarizedExperiment' object")
  if(!(is.numeric(th.fve)))
    stop("'th.fve' must be numeric")

  # check the presence of NA or Inf
  if (any(is.na(assay(data))))
    stop("NA values are not allowed in the 'data' matrix")
  if (any(is.infinite(assay(data))))
    stop("Inf values are not allowed in the 'data' matrix")

  # specific checks

  if (all(assay(data) == 0))
    stop("All genes have 0 counts.")
  if (all((assay(data) %%1) == 0))
    stop("It seems that you are using raw counts!
            This function works with normalized data")
  if(!("class" %in% colnames(colData(data))))
    stop("'class' info is lacking!
         Include the variable 'class'
         in colData(data) and label it 'class'!")
  if (th.fve >1 | th.fve < 0)
    stop("'th.fve must be between 0 and 1")

  data_filt<-assay(data)

  mod <- model.matrix(~data@colData$class)
  mod0 <- cbind(mod[,1])

  if (method == "fve"){
    invisible(capture.output(svaobj <- sva(data_filt, mod, mod0, round(ncol(data)/2))))
    pprob <- svaobj$pprob.gam*(1-svaobj$pprob.b)
    dats <- data_filt*pprob
    dats <- dats - rowMeans(dats)
    uu <- eigen(t(dats)%*%dats)
    #uu_val <-uu$values / sum(uu$values)
    uu_val2 <-uu$values^2 / sum(uu$values^2)
    x_val <- seq_len(ncol(dats))
    uu_val2_dec <- round(uu_val2, 3)
    pve_cumsum <- cumsum(uu_val2_dec)
    # https://support.bioconductor.org/p/88553/#88690
    expl_var_plot <- as.data.frame(cbind(x_val, uu_val2, pve_cumsum))
    n.sv <- order(which(pve_cumsum<=th.fve), decreasing = TRUE)[1]
    if(is.na(n.sv))
      stop("No SV identified. Increase the value of 'th.fve'")
    if (n.sv>round(ncol(data)/2)){
      n.sv<-round(ncol(data)/2)
    }
    print(ggplot(expl_var_plot, aes(x_val,pve_cumsum)) +
            geom_point(size=3, color="blue")  +
            geom_text(aes(label=rownames(expl_var_plot)), hjust=0.5, vjust=2, size=3) +
            geom_point(data = expl_var_plot[n.sv, ], aes(x_val, pve_cumsum), color="red", size=4) +
            xlab("SV") +
            ylab("Fraction of Variance Explained")+
            ggtitle("Fraction of Variance Explained"))
    # extract SVs
    sv_matrix <- as.matrix(svaobj$sv[, seq_len(n.sv)])
    cat("The number of SVs identified, which explain", th.fve*100, "% of Variance, is:", n.sv, "\n")
  } else if(method == "leek") {
    n.sv <- num.sv(data_filt, mod = mod, method = "leek")
    invisible(capture.output(svaobj <- sva(data_filt, mod, mod0, n.sv)))
    # extract SVs
    sv_matrix <- as.matrix(svaobj$sv)
    cat("The number of SVs identified by SVA (with method = 'leek')is:", n.sv, "\n")
  } else if(method == "be"){
    n.sv <- num.sv(data_filt , mod = mod, method = "be")
    invisible(capture.output(svaobj <- sva(data_filt, mod, mod0, n.sv)))
    # extract SVs
    sv_matrix <- as.matrix(svaobj$sv)
    cat("The number of SVs identified by SVA (with method = 'be')is:", n.sv, "\n")
  } else {
    stop("Please set 'fve', 'leek' or 'be' as SV identification method.")

  }

  return(sv_matrix)

}
