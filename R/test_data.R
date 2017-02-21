#' Example gene-expression dataset for DaMiRseq package
#'
#' A dataset with count matrix to test several DaMiRseq functions.
#' To show package functionality in a reasonably execution time,
#' sample data are a subset of Genotype-Tissue Expression (GTEx)
#' RNA-Seq database (dbGap Study Accession: phs000424.v6.p1). Samples incude
#' 20 Anterior Cingulate Cortex (ACC) tissues and 20 Frontal Cortex (FC) tissues.
#' 21363 genes have been preaviously selected to have 5 read counts
#' in at least 60% of samples.
#'
#' @format A SummarizedExperiment object containing an assay of
#' 21363 randomly selected genes (rows) and 40 samples (columns)
#' and a colData with 5 variables
#'
#' @return
#' An example dataset for \code{DaMiRseq} package
#'
"SE"

#' A dataset with a normalized matrix to test several DaMiRseq functions:
#' sample data are a subset of Genotype-Tissue Expression (GTEx)
#' RNA-Seq database (dbGap Study Accession: phs000424.v6.p1)
#'
#' @format A SummarizedExperiment object containing an assay of
#' 4897 genes (rows) and 40 samples (columns) and a colData with 5 variables
#'
#' @return
#' An example dataset for \code{DaMiRseq} package
#'
"data_norm"

#' Example gene-expression dataset for DaMiRseq package
#'
#' A data frame with class and covariates information
#'
#' @format A dataframe with 40 samples (rows) and 5 variables (columns):
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
"df"

#' Example Surrogate Variables dataset for DaMiRseq package
#'
#' A dataset with surrogate variables to test DaMiRseq functions
#'
#' @format A matrix with 40 samples (rows) and 4 surrogate variables (columns):
#' @return
#' An example dataset for \code{DaMiRseq} package
#'
"sv"

#' Example gene-expression dataset for DaMiRseq package
#'
#' A dataset with a small dimension of normalized expression data
#' in DaMiRseq package
#'
#' @format A list with:
#' \describe{
#'  \item{data}{reduced expression matrix}
#'  \item{variables}{a data frame with variables}
#' }
#' @return
#' An example dataset for \code{DaMiRseq} package
"data_reduced"

#' Example gene-expression dataset for DaMiRseq package
#'
#' A dataset with a small dimension of normalized expression data
#' in DaMiRseq package
#'
#' @format A data frame with 40 samples (rows) and 87 genes (columns)
#'
#' @return
#' An example dataset for \code{DaMiRseq} package
"data_min"

#' Example ranking dataset for DaMiRseq package
#'
#' A data frame with relieF and scaled reliefF scores for each gene
#'
#' @format A dataframe with 87 genes (rows) and 2 variables (columns):
#' \describe{
#'  \item{reliefF Score}{reliefF score for each gene}
#'  \item{scaled reliefF Score}{scaled reliefF score for each gene, by z-score}
#' }
#' @return
#' An example dataset for \code{DaMiRseq} package
"data_relief"

#' Example gene-expression dataset for DaMiRseq package
#'
#' A dataset with normalized expression data to build classification models
#' in DaMiRseq package
#'
#' @format A dataframe with 40 samples (rows) and 7 variables (genes):
#' @return
#' An example dataset for \code{DaMiRseq} package
"selected_features"

#' A sample dataset with a normalized count matrix for "testthat"
#' functions.
#'
#' @format A SummarizedExperiment object containing an assay of
#' 100 genes (rows) and 11 samples (columns) and a colData with 5 variables
#'
#' @return
#' An example dataset for \code{DaMiRseq} package
#'
"SEtest_norm"
