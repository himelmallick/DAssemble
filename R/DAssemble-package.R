#' DAssemble: Ensemble Models for Differential Analysis
#'
#' DAssemble provides an ensemble framework for differential-abundance and
#' differential-expression analysis across bulk RNA-seq, single-cell RNA-seq,
#' and microbiome studies. The package combines a user-selected core method
#' with optional enhancer tests and aggregates evidence using the Cauchy
#' Combination Test.
#'
#' @details
#' The main user-facing function is [DAssemble()]. It accepts either a
#' `MultiAssayExperiment` object or a feature table paired with sample metadata.
#' Core methods can be combined with enhancer tests including Wilcoxon rank-sum,
#' Kolmogorov-Smirnov, and presence-absence logistic regression.
#'
#' @seealso [DAssemble()]
#'
#' @keywords package
"_PACKAGE"
