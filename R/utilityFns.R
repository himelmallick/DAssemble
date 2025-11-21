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


DA_as_df <- function(x) {
  if (is.data.frame(x)) return(x)
  if (is.list(x) && !is.null(x$res) && is.data.frame(x$res)) return(x$res)
  stop("Unexpected fit output; expected a data.frame or a list with $res.")
}

DA_standardize_cols <- function(df, feature_pref = c("feature","feature"), p_pref, new_p_name) {
  df <- as.data.frame(df)
  gcol <- feature_pref[feature_pref %in% names(df)]
  if (length(gcol) == 0L) {
    if (!is.null(rownames(df))) df$feature <- rownames(df) else stop("No feature ID column found.")
  } else if (gcol[1] != "feature") {
    names(df)[names(df) == gcol[1]] <- "feature"
  }
  pcol <- p_pref[p_pref %in% names(df)]
  if (length(pcol) == 0L) stop("No p-value column found; candidates: ", paste(p_pref, collapse = ", "))
  names(df)[names(df) == pcol[1]] <- new_p_name
  df
}

# Safe join on "feature"
DA_join_by_feature <- function(x, y) {
  if (!("feature" %in% names(x)) || !("feature" %in% names(y))) stop("Both inputs must have 'feature'.")
  dplyr::left_join(x, y, by = "feature")
}

# Helper to ensure that feature matrix and metadata have matching rows (samples)
DA_align_features_metadata <- function(features, metadata) {
  if (!is.data.frame(features)) {
    features <- as.data.frame(features)
  }
  if (!is.data.frame(metadata)) {
    metadata <- as.data.frame(metadata)
  }
  if (is.null(rownames(features)) || is.null(rownames(metadata))) {
    stop("Both 'features' and 'metadata' must have rownames representing sample IDs.")
  }
  if (!setequal(rownames(features), rownames(metadata))) {
    stop("Rownames of 'features' and 'metadata' must match.")
  }
  common <- intersect(rownames(features), rownames(metadata))
  features <- features[common, , drop = FALSE]
  metadata <- metadata[common, , drop = FALSE]
  list(features = features, metadata = metadata)
}


## ======================================
## 1) Align counts & metadata (no norm)
## ======================================

DA_align_counts_and_md <- function(X, md,
                                   orientation = c("auto","features_by_samples","samples_by_features"),
                                   expVar, coVars = NULL) {
  orientation <- match.arg(orientation)
  X <- as.matrix(X)
  
  ## --- get sample IDs from metadata ---
  get_md_samples <- function(md) {
    rn <- rownames(md)
    if (!is.null(rn) && all(nzchar(rn))) return(rn)
    if ("sample" %in% names(md)) return(as.character(md$sample))
    rn
  }
  
  samp_ids <- get_md_samples(md)
  if (is.null(samp_ids) || !length(samp_ids) || any(!nzchar(samp_ids))) {
    stop("metadata must have rownames or a non-empty 'sample' column for matching.")
  }
  # ensure rownames are set
  if (is.null(rownames(md)) || any(!nzchar(rownames(md)))) {
    rownames(md) <- samp_ids
  }
  
  rn <- rownames(X)
  cn <- colnames(X)
  
  matches_row <- if (!is.null(rn)) sum(rn %in% samp_ids) else 0L
  matches_col <- if (!is.null(cn)) sum(cn %in% samp_ids) else 0L
  
  ## --- decide which dimension is samples, using names ---
  sample_dim <- NULL
  if (orientation == "samples_by_features") {
    sample_dim <- "row"
  } else if (orientation == "features_by_samples") {
    sample_dim <- "col"
  } else {  # auto
    if (matches_row > 0L || matches_col > 0L) {
      sample_dim <- if (matches_row >= matches_col) "row" else "col"
    } else {
      # fall back: assume rows are samples if we truly have no names
      sample_dim <- "row"
    }
  }
  
  ## --- subset & reorder by sample IDs; always end with samples x features ---
  if (sample_dim == "row") {
    if (is.null(rn)) stop("counts has no rownames to match metadata sample IDs.")
    common <- intersect(rn, samp_ids)
    if (!length(common)) {
      stop("No overlapping sample IDs between counts rownames and metadata.")
    }
    # align both counts and metadata to common sample set
    counts_sxf <- X[common, , drop = FALSE]  # samples x features
    md_s       <- md[common, , drop = FALSE]
    
  } else if (sample_dim == "col") {
    if (is.null(cn)) stop("counts has no colnames to match metadata sample IDs.")
    common <- intersect(cn, samp_ids)
    if (!length(common)) {
      stop("No overlapping sample IDs between counts colnames and metadata.")
    }
    X_sub <- X[, common, drop = FALSE]      # features x samples
    counts_sxf <- t(X_sub)                  # samples x features
    md_s       <- md[common, , drop = FALSE]
    
  } else {
    stop("Internal error: unknown sample_dim.")
  }
  
  # final check: rownames(counts) == rownames(metadata)
  stopifnot(identical(rownames(counts_sxf), rownames(md_s)))
  
  ## --- factorize expVar if binary (keeps numeric for continuous) ---
  if (!expVar %in% colnames(md_s)) {
    stop("expVar '", expVar, "' not found in metadata.")
  }
  if (is.character(md_s[[expVar]]) || is.logical(md_s[[expVar]]) ||
      (is.numeric(md_s[[expVar]]) && length(unique(md_s[[expVar]])) == 2L)) {
    md_s[[expVar]] <- droplevels(factor(md_s[[expVar]]))
  }
  
  storage.mode(counts_sxf) <- "integer"
  
  list(counts = counts_sxf, md = md_s)
}




## ===============================================
## 6) TPCC rows (use your CCT)
## ===============================================

# Vectorized row-wise CCT using your CCT() defined in your codebase.
DA_CCT_rows <- function(mat) {
  if (!exists("CCT")) stop("CCT not found in the environment.")
  apply(mat, 1L, function(pv) CCT(as.numeric(pv)))
}


## ===============================================
## 7) Metrics (fallback if you want them here)
## ===============================================

# If you want DAssemble to emit metrics immediately, this is a self-contained version.
# It does NOT change any modeling logic; it only uses the adjPval threshold you pass.
DA_compute_metrics <- function(res_df, truth, alpha, time_min, method_label) {
  res_df$adjPval[is.na(res_df$adjPval)] <- 1
  adj <- res_df$adjPval
  names(adj) <- res_df$feature
  
  if (is.null(truth)) {
    return(data.frame(Method = method_label,
                      Power = NA_real_, FPR = NA_real_, FDR = NA_real_,
                      F1 = NA_real_, MCC = NA_real_, ROC = NA_real_, pAUROC = NA_real_,
                      alpha = alpha, Time = time_min))
  }
  
  ## truth & prediction for confusion matrix
  truth_vec <- names(adj) %in% truth
  pred      <- adj <= alpha
  
  xtab <- table(
    factor(pred,      levels = c(TRUE, FALSE)),
    factor(truth_vec, levels = c(TRUE, FALSE))
  )
  cm <- caret::confusionMatrix(xtab, positive = "TRUE")
  
  Sens <- as.numeric(cm$byClass["Sensitivity"])
  Spec <- as.numeric(cm$byClass["Specificity"])
  FDR  <- 1 - as.numeric(cm$byClass["Precision"])
  F1   <- as.numeric(cm$byClass["F1"])
  MCC  <- mltools::mcc(preds = pred, actuals = truth_vec)
  
  ## ---- ROC / pAUROC ----
  response_roc <- factor(truth_vec, levels = c(FALSE, TRUE))
  
  # make safe scores: cap adj > 0, replace 0 with smallest positive
  adj_safe <- adj
  adj_safe[adj_safe <= 0] <- min(adj_safe[adj_safe > 0], na.rm = TRUE)
  score <- -log10(adj_safe)
  
  # drop non-finite scores for ROC only
  ok <- is.finite(score) & !is.na(score)
  response_roc_ok <- response_roc[ok]
  score_ok        <- score[ok]
  
  if (length(levels(response_roc_ok)) < 2 || length(unique(score_ok)) < 2) {
    ROC    <- NA_real_
    pAUROC <- NA_real_
  } else {
    # full ROC
    roc1 <- pROC::roc(
      response  = response_roc_ok,
      predictor = score_ok,
      direction = ">", quiet = TRUE
    )
    ROC <- as.numeric(roc1$auc)
    
    # partial AUC with correction, but safely
    pAUROC <- tryCatch({
      roc2 <- pROC::roc(
        response  = response_roc_ok,
        predictor = score_ok,
        direction = ">", quiet = TRUE,
        partial.auc = c(1, 0.85),
        partial.auc.correct = TRUE
      )
      as.numeric(roc2$auc)
    }, warning = function(w) {
      NA_real_
    }, error = function(e) {
      NA_real_
    })
    
    # enforce AUC >= 0.5 by flipping if needed
    if (!is.na(ROC) && ROC < 0.5) {
      ROC <- 1 - ROC
      if (!is.na(pAUROC)) pAUROC <- 1 - pAUROC
    }
  }
  
  data.frame(
    Method  = method_label,
    Power   = round(Sens, 3),
    FPR     = round(1 - Spec, 3),
    FDR     = round(FDR, 3),
    F1      = round(F1, 3),
    MCC     = round(MCC, 3),
    ROC     = round(ROC, 3),
    pAUROC  = round(pAUROC, 3),
    alpha   = alpha,
    Time    = time_min
  )
}


DA_norm_for_enhancers <- function(counts,
                                  domain = c("bulkrnaseq","singlecell","microbiome","none")) {
  domain <- match.arg(domain)
  counts <- as.matrix(counts)
  
  if (domain == "none") {
    # no normalization – just return input
    return(counts)
  }
  
  if (domain == "bulkrnaseq") {
    # TMM (edgeR) -> CPM (non-log)
    dge <- edgeR::DGEList(counts = counts)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    counts_norm <- edgeR::cpm(dge, log = FALSE)
    
  } else if (domain == "singlecell") {
    # SCRAN size factors -> divide counts by size factors
    if (!exists("SCRANnorm")) stop("SCRANnorm not found for domain = 'singlecell'.")
    sf <- SCRANnorm(counts)                   # length = ncol(counts)
    counts_norm <- t(t(counts) / sf)
    
  } else if (domain == "microbiome") {
    # CLR per sample (column)
    clr_fun <- function(x) {
      x <- as.numeric(x)
      x <- x + 1
      lx <- log(x)
      lx - mean(lx, na.rm = TRUE)
    }
    counts_norm <- apply(counts, 2L, clr_fun)
    if (is.vector(counts_norm)) {
      counts_norm <- matrix(counts_norm, nrow = nrow(counts))
    }
    rownames(counts_norm) <- rownames(counts)
    colnames(counts_norm) <- colnames(counts)
    
  } else {
    stop("Unknown domain: ", domain)
  }
  
  counts_norm
}




#### CCT P-value Combination ####

CCT <- function(pvals, weights=NULL){
  pvals <- ifelse(is.na(pvals), 1, pvals)
  
  #### check if all p-values are between 0 and 1
  
  if((sum(pvals<0) + sum(pvals>1)) > 0){
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
  
  #### check if there are very small non-zero p-value
  
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

#### Combine P-vals Function ####

combinePvalCCT<-function(pvalue1, pvalue2){
  
  if (length(pvalue1)!=length(pvalue2)){stop('the lengths of two pvalue vectors should match up!')}
  n<-length(pvalue1)
  combined_p_value<-c()
  for(i in 1:n){
    combined_p_value[i]<-CCT(c(pvalue1[i], pvalue2[i]))
  }
  return(combined_p_value)
}

#### Normalizing using Quantile ####

quan_norm <- function(features,metadata){
  norm.y <- preprocessCore::normalize.quantiles(as.matrix(features))
  norm.y <- data.frame(norm.y)
  names(norm.y) <- names(features)
  rownames(norm.y) <- rownames(features)
  return(norm.y)
}

#### Normalizing using Upper Quartile ####

uqrt_norm <- function(features, metadata){
  quant.exp <- apply(as.matrix(features), 2, function(x){quantile(x[x > 0], 0.75)})
  norm.y <- data.frame(t(t(as.matrix(features))/quant.exp))
  return(norm.y)
}

#### Normalizing using using CPM  ####

cpm_norm <- function(features, metadata, const.mult = 1e+06, prior.count = 1, log = TRUE){
  norm.factors <- edgeR::calcNormFactors(features)*colSums(features)
  cpm <- features %*% diag(1/(norm.factors)) * const.mult
  colnames(cpm) <- colnames(features)
  log(cpm + prior.count)
}

#### Normalizing using TMM (edgeR) ####

tmm_norm <- function(features, metadata){
  norm.y <- edgeR::DGEList(features)
  norm.y <- edgeR::calcNormFactors(norm.y, method = "TMM")
  norm.y <- as.data.frame(edgeR::cpm(norm.y, log = FALSE))
  return(norm.y)
}

#### Normalizing using Geometric Means (DESeq2) ####

rle_norm <- function(features, metadata, coVars = NULL, expVar = 'Exposure'){
  if(is.null(coVars)){
    metadata <- metadata[, c(expVar), drop = FALSE]
  }else{
    metadata <- metadata[, c(expVar, coVars)]
  }
  formula <- as.formula(paste('~', paste(colnames(metadata), collapse = "+"), sep = ''))
  x <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(features),
                                                       colData = metadata,
                                                       design = formula))
  gm_mean <- function(x, na.rm = TRUE){
    exp(sum(log(x[x > 0]), na.rm = na.rm)/length(x))
  }
  geoMeans <- apply(counts(x), 1, gm_mean)
  s <- DESeq2::estimateSizeFactors(x,geoMeans = geoMeans)
  s <- s$sizeFactor
  norm.y <- data.frame(t(apply(features, 1, function(x)x/s)))
  return(norm.y)
}


DA_build_subensembles <- function(res, core_method, enh_list, p_adj = "BH") {
  if (!("feature" %in% names(res))) {
    stop("`res` must contain a 'feature' column.")
  }
  
  # Map method labels to p-value column names in res
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
  
  # Generate all non-empty subsets of the available labels
  subsets <- list()
  for (k in seq_len(length(labels))) {
    cmb <- combn(labels, k, simplify = FALSE)
    subsets <- c(subsets, cmb)
  }
  
  ensembles <- list()
  feature <- res$feature
  
  for (lab_set in subsets) {
    # Define a convenient name: replace "core" by the actual core name if provided
    name_parts <- vapply(lab_set, function(x) {
      if (x == "core") {
        if (is.null(core_method)) "core" else core_method
      } else {
        x
      }
    }, character(1))
    key <- paste(name_parts, collapse = "+")
    
    pcols <- unlist(method_to_pcol[lab_set], use.names = FALSE)
    pmat  <- as.matrix(res[, pcols, drop = FALSE])
    
    # Combine p-values with CCT; if there is only one column this reduces to identity
    p_comb <- DA_CCT_rows(pmat)
    
    ensembles[[key]] <- data.frame(
      feature   = feature,
      Pval    = p_comb,
      adjPval = stats::p.adjust(p_comb, method = p_adj),
      stringsAsFactors = FALSE
    )
  }
  
  ensembles
}