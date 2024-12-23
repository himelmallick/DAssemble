#' Combine P-Value Data Frames Using Different Methods
#'
#' The `DAssemble` function combines multiple p-value data frames using specified methods
#' (Stouffer's method, CCT, Fisher, etc.) and applies corrections to compute adjusted p-values.
#'
#' @param dflist A list of data frames to combine. Each data frame must have columns named `ID` and `pvalue`.
#' @param combine.method A character string specifying the p-value combination method.
#'   Must be one of `"stouffer"`, `"CCT"`, `"fisher"`, `"invchisq"`, `"binomtest"`, `"bonferroni"`,
#'   `"snippet"`, `"harmonic"`, or `"minP"`. Default is `"stouffer"`.
#' @param correction A character string specifying the method for p-value correction.
#'   Default is `"BH"` (Benjamini-Hochberg). Other valid values are `"bonferroni"`, `"holm"`, `"hochberg"`,
#'   `"hommel"`, `"BH"`, `"BY"`, or `"none"`.
#'
#' @details
#' - The input data frames are merged by the `ID` column.
#' - Rows with missing (`NA`) or abnormal p-values (outside the range 0 to 1) are removed.
#' - Different combination methods are applied based on the `combine.method` parameter.
#'
#' @return A data frame with the following columns:
#' - `ID`: IDs shared across all input data frames.
#' - Combined p-values (as a new column named `pval.combined`).
#' - Additional columns for corrected p-values based on the `correction` parameter.
#'
#' @examples
#' pvec1 <- data.frame(ID = c('gene1', 'gene2', 'gene3', 'gene4'),
#'                     pvalue = c(0.01, 0.09, 0.02, 0.06))
#'
#' pvec2 <- data.frame(ID = c('gene1', 'gene2', 'gene3', 'gene4'),
#'                     pvalue = c(0.02, 0.03, 0.04, 0.07))
#'
#' pvec3 <- data.frame(ID = c('gene1', 'gene2', 'gene3', 'gene4'),
#'                     pvalue = c(0.05, 0.01, 0.03, 0.08))
#'
#' dflist <- list(pvec1, pvec2, pvec3)
#'
#' DAssemble(dflist, combine.method = 'harmonic', correction = 'hommel')
#' @export
DAssemble <- function(dflist, combine.method = "stouffer", correction = "BH"){

  if(!combine.method %in% c("stouffer", "CCT", "fisher", "invchisq", "binomtest", "bonferroni", "snippet", "harmonic", "minP")) stop("combine.method must be one of stouffer, CCT, invchisq, binomtest, bonferroni, snippet, harmonic, minP")
  if(!correction %in% c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "none")) stop("correction must be one of bonferroni, holm, hochberg, hommel, BH, BY, none")

  if (!requireNamespace("poolr", quietly = TRUE)) {
    stop("Package 'poolr' is required. Please install it.")
  }

  if (!requireNamespace("metapod", quietly = TRUE)) {
    stop("Package 'metapod' is required. Please install it.")
  }

  if (!requireNamespace("harmonicmeanp", quietly = TRUE)) {
    stop("Package 'harmonicmeanp' is required. Please install it.")
  }

  suppressPackageStartupMessages(library("poolr"))
  suppressPackageStartupMessages(library("metapod"))
  suppressPackageStartupMessages(library("harmonicmeanp"))

  # Ensure at least two dflist are provided
  if (length(dflist) < 2) {
    stop("At least two pvec data frames must be provided.")
  }

  # Check that each dflist has 'ID' and 'pvalue' columns
  for (i in seq_along(dflist)) {
    if (!all(c('ID', 'pvalue') %in% colnames(dflist[[i]]))) {
      stop(paste("Table", i, "must have 'ID' and 'pvalue' columns"))
    }
    # Rename the 'pvalue' column to have unique names
    colnames(dflist[[i]])[colnames(dflist[[i]]) == 'pvalue'] <- paste0('pvalue', i)
  }

  # Merge all dflist on 'ID'
  p.combined <- Reduce(function(x, y) merge(x, y, by = 'ID', all = FALSE), dflist)

  # Check if there are any IDs left after merging
  if (nrow(p.combined) == 0) {
    stop("No common IDs found across all dflist.")
  }

  # Collect all the p-values into a matrix
  pval_columns <- grep('pvalue', colnames(p.combined))
  pvals_matrix <- as.matrix(p.combined[, pval_columns])

  # Check for NA values in the p-value matrix
  if (anyNA(pvals_matrix)) {
    warning("There are NA values in the p-value matrix. IDs with NA p-values will be removed.")
    # Remove rows with any NA p-values
    na_rows <- apply(pvals_matrix, 1, function(x) any(is.na(x)))
    p.combined <- p.combined[!na_rows, ]
    pvals_matrix <- pvals_matrix[!na_rows, ]
  }

  # Check for abnormal p-values (outside the range [0, 1])
  if (any(pvals_matrix < 0 | pvals_matrix > 1)) {
    stop("Abnormal p-values detected: p-values must be between 0 and 1.")
  }

  # p-value combination
  if(combine.method=='stouffer'){
    pvals_list <- split(pvals_matrix, col(pvals_matrix))
    # Combine p-values using the specified method
    p.combined$pval.combined<-metapod::combineParallelPValues(pvals_list, method = combine.method)$p.value

  } else if (combine.method=='CCT'){
    p.combined$pval.combined<-combinePvalCCT(pvals_matrix)

  } else if (combine.method=='fisher'){
    pval.combined <- NULL
    for(i in 1 : nrow(pvals_matrix)){
      pval.combined[i] <- poolr::fisher(pvals_matrix[i,])$p
    }
    p.combined$pval.combined <- pval.combined

  } else if (combine.method=='invchisq'){
    pval.combined <- NULL
    for(i in 1 : nrow(pvals_matrix)){
      pval.combined[i] <- poolr::invchisq(pvals_matrix[i,])$p
    }
    p.combined$pval.combined <- pval.combined

  } else if (combine.method=='binomtest'){
    pval.combined <- NULL
    for(i in 1 : nrow(pvals_matrix)){
      pval.combined[i] <- poolr::binomtest(pvals_matrix[i,])$p
    }
    p.combined$pval.combined <- pval.combined

  } else if (combine.method=='bonferroni'){
    pval.combined <- NULL
    for(i in 1 : nrow(pvals_matrix)){
      pval.combined[i] <- poolr::bonferroni(pvals_matrix[i,])$p
    }
    p.combined$pval.combined <- pval.combined

  } else if (combine.method=='snippet'){
    pval.combined <- NULL
    for(i in 1 : nrow(pvals_matrix)){
      pval.combined[i] <- poolr::snippet(pvals_matrix[i,])$p
    }
    p.combined$pval.combined <- pval.combined

  } else if (combine.method=='harmonic'){
    pval.combined <- NULL
    for(i in 1 : nrow(pvals_matrix)){
      pval.combined[i] <- harmonicmeanp::p.hmp(pvals_matrix[i,], L = ncol(pvals_matrix))
    }
    p.combined$pval.combined <- pval.combined

  } else{
    ## minP
    p.combined$pval.combined<-apply(pvals_matrix, 1, min, na.rm = FALSE)
  }

  p.combined<-append_qvalues(p.combined, correction)
  p.combined<-p.combined[order(p.combined[,ncol(p.combined)]),]

  return(p.combined)
}
