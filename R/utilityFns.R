###########################
# DAssemble Normalization #
###########################

DAssemble_normalize_features <- function(features,
                                         metadata,
                                         method = c("scran", "TMM", "CLR", "TSS"),
                                         return_normalized = TRUE,
                                         pseudocount = 1) {
  
  ##########################
  # Match and Sanity Check #
  ##########################
  
  method <- match.arg(
    toupper(method),
    choices = c("SCRAN", "TMM", "CLR", "TSS")
  )
  
  if (is.null(rownames(features)) ||
      is.null(rownames(metadata)) ||
      !identical(rownames(features), rownames(metadata))) {
    stop("`features` and `metadata` must have identical rownames (samples).")
  }
  
  if (!is.data.frame(features)) {
    features <- as.data.frame(features)
  }
  
  ######################
  # Initial Extraction #
  ######################
  
  X_raw  <- as.matrix(features)     # preserve raw counts
  X      <- X_raw
  n_samp <- nrow(X)
  
  sf <- rowSums(X_raw)              # default: library size (used for TSS and CLR metadata)
  names(sf) <- rownames(X)
  
  #########################
  # Perform Normalization #
  #########################
  
  if (method == "SCRAN") {
    
    if (!requireNamespace("scran", quietly = TRUE)) {
      stop("The 'scran' package is required for method = 'scran'.")
    }
    
    sf <- scran::calculateSumFactors(t(X))
    
  } else if (method == "TMM") {
    
    if (!requireNamespace("edgeR", quietly = TRUE)) {
      stop("The 'edgeR' package is required for method = 'TMM'.")
    }
    
    dge <- edgeR::DGEList(counts = t(X))
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    eff_lib <- dge$samples$lib.size * dge$samples$norm.factors
    sf <- eff_lib / mean(eff_lib)
    
  } else if (method == "TSS") {
    
    # Total Sum Scaling: sf already equals library size (row sums of raw counts)
    # Keep sf as-is.
    
  } else if (method == "CLR") {
    
    ######################
    # CLR Transformation #
    ######################
    
    X_pc <- X + pseudocount
    logX <- log(X_pc)
    X    <- logX - rowMeans(logX, na.rm = TRUE)
    # sf remains original library size (rowSums of X_raw)
  }
  
  #################
  # Apply Scaling #
  #################
  
  if (return_normalized && method %in% c("SCRAN", "TMM", "TSS")) {
    X <- X / sf
  }
  
  ################################
  # Attach scale factor only     #
  # when returning UN-normalized #
  ################################
  
  if (!return_normalized) {
    metadata$scale_factor <- sf
  }
  
  ##########################
  # Return final output    #
  ##########################
  
  list(
    features = as.data.frame(X),
    metadata = metadata
  )
}

###########################
# DAssemble Call Function #
###########################

DA_do_call <- function(fun_or_name, args_list) {
  fn <- if (is.character(fun_or_name)) get(fun_or_name, mode = "function") else fun_or_name
  fm <- names(formals(fn))
  use <- args_list[intersect(names(args_list), fm)]
  missing_req <- setdiff(fm, names(use))
  if (length(missing_req)) for (m in missing_req) use[[m]] <- NULL
  do.call(fn, use)
}



# Vectorized row-wise CCT using your CCT() defined in your codebase.
DA_CCT_rows <- function(mat) {
  if (!exists("CCT")) stop("CCT not found in the environment.")
  apply(mat, 1L, function(pv) CCT(as.numeric(pv)))
}


##############################
## CCT: Cauchy Combination  ##
##############################

CCT <- function(pvals, weights = NULL) {
  # Replace NAs by 1 (no evidence)
  pvals <- ifelse(is.na(pvals), 1, pvals)
  
  # 1) Check all p between 0 and 1
  if ((sum(pvals < 0) + sum(pvals > 1)) > 0) {
    stop("All p-values must be between 0 and 1!")
  }
  
  # 2) Handle exact 0 or 1
  is.zero <- (sum(pvals == 0) >= 1)
  is.one  <- (sum(pvals == 1) >= 1)
  if (is.zero && is.one) {
    stop("Cannot have both 0 and 1 p-values!")
  }
  if (is.zero) {
    return(0)
  }
  if (is.one) {
    warning("There are p-values that are exactly 1!")
    return(1)
  }
  
  # 3) Weights: default = equal, then normalized
  if (is.null(weights)) {
    weights <- rep(1 / length(pvals), length(pvals))
  } else if (length(weights) != length(pvals)) {
    stop("The length of weights should be the same as that of the p-values!")
  } else if (sum(weights < 0) > 0) {
    stop("All the weights must be positive!")
  } else {
    weights <- weights / sum(weights)
  }
  
  # 4) Handle very small p-values to avoid overflow in tan()
  is.small <- (pvals < 1e-16)
  if (sum(is.small) == 0) {
    cct.stat <- sum(weights * tan((0.5 - pvals) * pi))
  } else {
    cct.stat <- sum((weights[is.small] / pvals[is.small]) / pi)
    cct.stat <- cct.stat + sum(weights[!is.small] * tan((0.5 - pvals[!is.small]) * pi))
  }
  
  # 5) Tail approximation when cct.stat is huge
  if (cct.stat > 1e+15) {
    pval <- (1 / cct.stat) / pi
  } else {
    pval <- 1 - pcauchy(cct.stat)
  }
  return(pval)
}

######################################
## Row-wise CCT with single-method  ##
## bypass when ncol(mat) == 1      ##
######################################

DA_CCT_rows <- function(mat) {
  if (!exists("CCT")) stop("CCT not found in the environment.")
  
  # If only one method contributes, return its p-values unchanged
  if (ncol(mat) == 1L) {
    return(as.numeric(mat[, 1L]))
  }
  
  # Otherwise, combine via CCT row-wise
  apply(mat, 1L, function(pv) CCT(as.numeric(pv)))
}

##########################################
## Two-vector helper (exactly 2 views)  ##
##########################################

combinePvalCCT <- function(pvalue1, pvalue2) {
  if (length(pvalue1) != length(pvalue2))
    stop("The lengths of two pvalue vectors should match up!")
  
  n <- length(pvalue1)
  combined_p_value <- numeric(n)
  
  for (i in seq_len(n)) {
    combined_p_value[i] <- CCT(c(pvalue1[i], pvalue2[i]))
  }
  combined_p_value
}

#######################################################
## Flexible helper for 1, 2, 3, ... p-value vectors ##
#######################################################

combinePvalCCT_multi <- function(...) {
  p_list <- list(...)
  if (length(p_list) == 0L) {
    stop("At least one p-value vector is required.")
  }
  
  # Check equal length
  lens <- vapply(p_list, length, integer(1))
  if (length(unique(lens)) != 1L) {
    stop("All p-value vectors must have the same length.")
  }
  
  # Build matrix and reuse DA_CCT_rows (which handles 1-column bypass)
  pmat <- do.call(cbind, p_list)
  DA_CCT_rows(pmat)
}

#########################################################
## Build subensembles using CCT (with single-method    ##
## bypass handled inside DA_CCT_rows)                  ##
#########################################################

DA_build_subensembles <- function(res, core_method, enh_list, p_adj = "BH") {
  if (!"feature" %in% names(res)) {
    stop("The input 'res' must contain a 'feature' column.")
  }
  
  method_to_pcol <- list()
  if (!is.null(core_method) && "Pval_core" %in% names(res)) {
    method_to_pcol$core <- "Pval_core"
  }
  if ("WLX" %in% enh_list && "Pval_WLX" %in% names(res)) {
    method_to_pcol$WLX <- "Pval_WLX"
  }
  if ("LR" %in% enh_list && "Pval_LR" %in% names(res)) {
    method_to_pcol$LR <- "Pval_LR"
  }
  if ("KS" %in% enh_list && "Pval_KS" %in% names(res)) {
    method_to_pcol$KS <- "Pval_KS"
  }
  
  labels <- names(method_to_pcol)
  if (length(labels) == 0L) {
    return(list())
  }
  
  # All non-empty subsets of available methods
  subsets <- list()
  for (k in seq_len(length(labels))) {
    cmb <- combn(labels, k, simplify = FALSE)
    subsets <- c(subsets, cmb)
  }
  
  ensembles <- list()
  feature <- res$feature
  
  for (lab_set in subsets) {
    # Build a human-readable key like "core+WLX+LR"
    name_parts <- vapply(
      lab_set,
      function(x) {
        if (x == "core") {
          if (is.null(core_method)) "core" else core_method
        } else {
          x
        }
      },
      character(1)
    )
    key <- paste(name_parts, collapse = "+")
    
    pcols <- unlist(method_to_pcol[lab_set], use.names = FALSE)
    pmat  <- as.matrix(res[, pcols, drop = FALSE])
    
    # Combine p-values row-wise with CCT (or identity if one column)
    p_comb <- DA_CCT_rows(pmat)
    
    ensembles[[key]] <- data.frame(
      feature = feature,
      Pval    = p_comb,
      adjPval = stats::p.adjust(p_comb, method = p_adj),
      stringsAsFactors = FALSE
    )
  }
  
  ensembles
}