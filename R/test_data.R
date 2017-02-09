#' Example gene-expression dataset for DaMiRseq package
#'
#' A dataset with count matrix to test several DaMiRseq functions
#'
#' @format A matrix with 13998 genes (rows) and 80 samples (columns):
#'
#' @return
#' An example dataset for \code{DaMiRseq} package
#'
"DaMiR.Data.Ex.Norm"

#' Example gene-expression dataset for DaMiRseq package
#'
#' A dataset with class and covariates information
#'
#' @format A dataframe with 80 samples (rows) and 5 variables(columns):
#' \describe{
#'  \item{center}{center where sample has been collected}
#'  \item{sex}{sample's gender}
#'  \item{age}{sample's age}
#'  \item{death}{kind of sample's death, based on Hardy scale}
#'  \item{class}{sample's class}
#' }
#'
#' @return
#' An example dataset for \code{DaMiRseq} package
#'
"DaMiR.Data.Ex.dfClCov"

#' Example ranking dataset for DaMiRseq package
#'
#' A dataset with relieF and scaled reliefF scores for each gene
#'
#' @format A dataframe with 160 genes (rows) and 2 variables(columns):
#' \describe{
#'  \item{reliefF Score}{reliefF score for each gene}
#'  \item{scaled reliefF Score}{scaled reliefF score for each gene, by z-score}
#' }
#' @return
#' An example dataset for \code{DaMiRseq} package
"DaMiR.Data.Ex.relief"

#' Example Surrogate Variables dataset for DaMiRseq package
#'
#' A dataset with Surrogate Variables to test DaMiRseq functions
#'
#' @format A matrix with 80 samples (rows) and 12 Surrogate Variables (columns):
#' @return
#' An example dataset for \code{DaMiRseq} package
#'
"DaMiR.Data.Ex.SV"

#' Example gene-expression dataset for DaMiRseq package
#'
#' A dataset with normalized values to test classification in DaMiRseq package
#'
#' @format A dataframe with 80 samples (rows) and 10 variables(genes):
#' @return
#' An example dataset for \code{DaMiRseq} package
"DaMiR.Data.Ex.classif"

#' Example gene-expression dataset for DaMiRseq package
#'
#' A dataset with normalized values to test feature selection in DaMiRseq package
#'
#' @format A dataframe with 80 samples (rows) and 160 variables(genes):
#' @return
#' An example dataset for \code{DaMiRseq} package
"DaMiR.Data.Ex.t.Norm"

