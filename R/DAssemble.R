#' Ensemble models for differential analysis
#'
#' `DAssemble`  is a lightweight implementation of the stacking method as applied to two or more
#' differential analysis (DA) results tables. It takes as inputs a list of DA tables with p-values
#' and returns a single table with omnibus p-values and q-values on a per-feature basis.
#' Several p-value combination methods such as Stouffer's method, CCT, and Fisher are supported.
#'
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
#' - `pvalue` columns are validated to ensure they are numeric. Non-numeric columns are coerced
#'   to numeric where possible; otherwise, an error is raised.
#' - Rows with missing (`NA`) or abnormal p-values (outside the range 0 to 1) are removed.
#' - The `combine.method` parameter determines the method used for combining p-values.
#' - The `correction` parameter applies the specified multiple testing correction method
#'   to the combined p-values.
#'
#' @return A data frame with the following columns:
#' - `ID`: IDs shared across all input data frames.
#' - Combined p-values (as a new column named `pval.combined`).
#' - Additional columns for corrected p-values based on the `correction` parameter.
#' @examples
#' # Install packages if not already installed (uncomment if needed)
#' # install.packages('tidyverse')
#' # BiocManager::install("airway")
#' # BiocManager::install("edgeR")
#' # BiocManager::install("DESeq2")
#' # BiocManager::install("limma")
#'
#' library(tidyverse)
#' library(airway)
#' library(edgeR)
#' library(DESeq2)
#' library(limma)
#' library(DAssemble)
#'
#' #############
#' # Load Data #
#' #############
#' data("airway")
#'
#' # Filter for genes with at least 10 reads in at least 10 samples
#' filter <- filterByExpr(airway)
#' filtered <- airway[filter, ]
#' counts_full <- assay(filtered)
#' sample_info <- colData(filtered)
#'
#' # Subset to 1000 random genes
#' set.seed(123)
#' genes_subset <- sample(rownames(counts_full), size = 1000)
#' counts_sub <- counts_full[genes_subset, ]
#' treatment <- factor(sample_info$dex, levels = c("untrt", "trt"))
#'
#' #############
#' # Run edgeR #
#' #############
#' dge <- DGEList(counts = counts_sub, group = treatment)
#' dge <- calcNormFactors(dge)
#' design_edgeR <- model.matrix(~ treatment)
#' dge <- estimateDisp(dge, design_edgeR)
#' fit_edgeR <- glmQLFit(dge, design_edgeR)
#' qlf_edgeR <- glmQLFTest(fit_edgeR, coef = 2)
#' res_edgeR <- topTags(qlf_edgeR, n = Inf)
#' res_edgeR_df <- as.data.frame(res_edgeR)
#'
#' ##############
#' # Run DESeq2 #
#' ##############
#' dds <- DESeqDataSetFromMatrix(
#'   countData = counts_sub,
#'   colData = as.data.frame(sample_info),
#'   design = ~ dex
#' )
#' dds <- DESeq(dds)
#' res_DESeq2 <- results(dds, contrast = c("dex", "trt", "untrt"))
#' res_DESeq2_df <- as.data.frame(res_DESeq2)
#'
#' #################
#' # Run limmavoom #
#' #################
#' design_limma <- model.matrix(~ treatment)
#' v <- voom(counts_sub, design_limma, plot = FALSE)
#' fit_limma <- lmFit(v, design_limma)
#' fit_limma <- eBayes(fit_limma)
#' res_limma <- topTable(fit_limma, coef = 2, n = Inf)
#' res_limma_df <- as.data.frame(res_limma)
#'
#' ###############################################
#' # Convert the results to desired input format #
#' ###############################################
#' res_edgeR_df <- res_edgeR_df %>%
#'   rownames_to_column("ID") %>%
#'   dplyr::rename(pvalue = PValue)
#'
#' res_DESeq2_df <- res_DESeq2_df %>%
#'   rownames_to_column("ID")
#'
#' res_limma_df <- res_limma_df %>%
#'   rownames_to_column("ID") %>%
#'   dplyr::rename(pvalue = P.Value)
#'
#' ################
#' # Run DAssemble#
#' ################
#' dflist <- list(res_edgeR_df, res_DESeq2_df, res_limma_df)
#' assemble_res <- DAssemble(dflist, combine.method = "stouffer", correction = "BH")
#' head(assemble_res)
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

# Internal Helper Function 1
#' Append Adjusted P-Values to Data Frame
#'
#' Internal helper function to append corrected p-values to the combined data frame.
#' @param pvals_matrix A data frame or matrix containing combined p-values.
#' @param correction A string specifying the p-value correction method.
#' @keywords internal

append_qvalues<-function(pvals_matrix, correction = "BH"){

  if(correction == "bonferroni"){
    pvals_matrix$qval_bonferroni <- tryCatch(p.adjust(pvals_matrix$pval.combined, method = 'bonferroni'),
                                             error = function(err){NA})
  }

  if(correction == "holm"){
    pvals_matrix$qval_holm <- tryCatch(p.adjust(pvals_matrix$pval.combined, method = 'holm'),
                                       error = function(err){NA})
  }

  if(correction == "hochberg"){
    pvals_matrix$qval_hochberg <- tryCatch(p.adjust(pvals_matrix$pval.combined, method = 'hochberg'),
                                           error = function(err){NA})
  }

  if(correction == "hommel"){
    pvals_matrix$qval_hommel <- tryCatch(p.adjust(pvals_matrix$pval.combined, method = 'hommel'),
                                         error = function(err){NA})
  }

  if(correction == "BH"){
    pvals_matrix$qval_BH <- tryCatch(p.adjust(pvals_matrix$pval.combined, method = 'BH'),
                                     error = function(err){NA})
  }

  if(correction == "BY"){
    pvals_matrix$qval_BY <- tryCatch(p.adjust(pvals_matrix$pval.combined, method = 'BY'),
                                     error = function(err){NA})
  }

  if(correction == "none"){
    pvals_matrix$qval_none <- tryCatch(p.adjust(pvals_matrix$pval.combined, method = 'none'),
                                       error = function(err){NA})
  }

  return(pvals_matrix)
}

# Internal Helper Function 2
#' Combine P-Values Using CCT
#'
#' Internal helper function to combine p-values using the CCT method.
#' @param pvals A vector of p-values.
#' @param weights Optional weights for the CCT method.
#' @keywords internal
CCT <- function(pvals, weights=NULL){
  #### check if there is NA
  pvals <- ifelse(is.na(pvals), 1, pvals)
  # if(sum(is.na(pvals)) > 0){
  #   stop("Cannot have NAs in the p-values!")
  # }

  # Check if there are any p-values left
  if (length(pvals) == 0) {
    return(NA)
  }


  #### check if all p-values are between 0 and 1
  if (any(pvals < 0 | pvals > 1)) {
    stop("All p-values must be between 0 and 1!")
  }

  #### check if there are p-values that are either exactly 0 or 1.
  is.zero <- (sum(pvals==0)>=1)
  is.one <- (sum(pvals==1)>=1)
  if(is.zero && is.one){
    stop("Cannot have both 0 and 1 p-values!")
  }
  if(is.zero){
    return(0)
  }
  if(is.one){
    warning("There are p-values that are exactly 1!")
    return(1)
  }

  #### check the validity of weights (default: equal weights) and standardize them.
  if(is.null(weights)){
    weights <- rep(1/length(pvals),length(pvals))
  }else if(length(weights)!=length(pvals)){
    stop("The length of weights should be the same as that of the p-values!")
  }else if(sum(weights < 0) > 0){
    stop("All the weights must be positive!")
  }else{
    weights <- weights/sum(weights)
  }

  #### check if there are very small non-zero p-values
  is.small <- (pvals < 1e-16)
  if (sum(is.small) == 0){
    cct.stat <- sum(weights*tan((0.5-pvals)*pi))
  }else{
    cct.stat <- sum((weights[is.small]/pvals[is.small])/pi)
    cct.stat <- cct.stat + sum(weights[!is.small]*tan((0.5-pvals[!is.small])*pi))
  }

  #### check if the test statistic is very large.
  if(cct.stat > 1e+15){
    pval <- (1/cct.stat)/pi
  }else{
    pval <- 1-pcauchy(cct.stat)
  }
  return(pval)
}

combinePvalCCT<-function(pvals_matrix){

  # Apply the CCT function to each row
  combined_p_value <- apply(pvals_matrix, 1, CCT)

  return(combined_p_value)
}
