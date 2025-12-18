#' Differential analysis using DAssemble
#'
#' The `DAssemble` function implements an ensemble strategy for differential
#' analysis. Users can select a single *core* model and optionally augment it
#' with up to three *enhancer* models. Enhancers are simple, often nonparametric
#' tests that provide complementary evidence. P-values from the core and
#' enhancer models are combined using the Cauchy Combination Test (CCT).
#'
#' @export
DAssemble <- function(features,
                      metadata,
                      core_method = NULL,
                      enhancers = NULL,
                      expVar = "group",
                      p_adj = "BY",
                      enhancer_norm = "TSS",
                      return_components   = TRUE,
                      return_subensembles = FALSE) {
  
  ##################################
  # 0) Track time and match domain #
  ##################################
  
  start_time  <- proc.time()
  norm_method <- match.arg(
    toupper(enhancer_norm),
    choices = c("SCRAN", "TMM", "CLR", "TSS")
  )
  norm_method <- tolower(norm_method)  
  
  ################################
  # 1) Validate core + enhancers #
  ################################
  
  valid_core <- c(
    "DESeq2", "edgeR", "limmaVOOM", "dearseq",
    "metagenomeSeq", "MAST", "Tweedieverse",
    "Maaslin2", "Maaslin3", "LOCOM",
    "LinDA", "ANCOMBC2", "Robseq", "ALDEx2"
  )
  
  if (is.null(core_method) || identical(core_method, "none")) {
    core_method <- NULL
  }
  if (!is.null(core_method) && !(core_method %in% valid_core)) {
    stop("Unknown core_method: ", core_method)
  }
  
  enh_list  <- if (is.null(enhancers)) character() else unique(enhancers)
  valid_enh <- c("WLX", "LR", "KS")
  
  if (length(enh_list) > length(valid_enh)) {
    stop("Too many enhancers specified. Available enhancers are: ",
         paste(valid_enh, collapse = ", "))
  }
  if (length(setdiff(enh_list, valid_enh)) > 0L) {
    stop("Unknown enhancers: ",
         paste(setdiff(enh_list, valid_enh), collapse = ", "))
  }
  if (is.null(core_method) && length(enh_list) == 0L) {
    stop("Specify core_method or at least one enhancer.")
  }
  
  #############################################
  # Warning for hurdle / zero-inflated models #
  #############################################
  
  zi_methods <- c("metagenomeSeq", "MAST", "Maaslin3")
  
  if (!is.null(core_method) &&
      core_method %in% zi_methods &&
      length(enh_list) > 0L) {
    warning(core_method,
            " is a hurdle / zero-inflated model; ensemble enhancers are not recommended ",
            "because these models already perform multi-part modeling. ",
            "Results may be harder to interpret.")
  }
  
  ##############################################
  # 2) Sanity check of input data and metadata #
  ##############################################
  
  if (!is.data.frame(features) || !is.data.frame(metadata)) {
    stop("`features` and `metadata` must both be data.frames.")
  }
  
  if (is.null(rownames(features)) || is.null(rownames(metadata)) ||
      !identical(rownames(features), rownames(metadata))) {
    stop("`features` and `metadata` must have identical rownames.")
  }
  
  ########################
  # 3) Check expVar type #
  ########################
  
  if (!expVar %in% colnames(metadata)) {
    stop("expVar '", expVar, "' not found in metadata.")
  }
  
  grp <- metadata[[expVar]]
  if (!is.factor(grp)) grp <- factor(grp)
  if (nlevels(grp) != 2L) {
    stop("expVar must be binary.")
  }
  metadata[[expVar]] <- droplevels(grp)
  
  ###############################################
  # 4) Tweedieverse: optional special norm      #
  ###############################################
  
  if (!is.null(core_method) && core_method == "Tweedieverse") {
    norm_out <- DAssemble_normalize_features(
      features          = features,
      metadata          = metadata,
      method            = 'scran',
      return_normalized = FALSE
    )
    features <- norm_out$features
    metadata <- norm_out$metadata
  }
  
  ##################
  ## 5) Core model #
  ##################
  
  components <- list()
  res        <- NULL
  
  if (!is.null(core_method)) {
    core_fun <- paste0("DA_fit_core_", core_method)
    
    res_core <- DA_do_call(core_fun, list(
      features = features,
      metadata = metadata,
      expVar   = expVar
    ))
    
    if (!"feature" %in% names(res_core)) {
      stop("Core result must contain a 'feature' column.")
    }
    if (!"pval_core" %in% names(res_core)) {
      stop("Core result must contain a 'pval_core' column.")
    }
    
    res <- res_core
    components[[core_method]] <- res_core
  }
  
  #################
  ## 6) Enhancers #
  #################
  
  if (length(enh_list) > 0L) {
    
    # Default normalization for enhancers (user choice)
    feat_norm_default <- DAssemble_normalize_features(
      features          = features,
      metadata          = metadata,
      method            = norm_method,
      return_normalized = TRUE
    )
    
    # Hard code for CLR/TSS + LR: nothing happens except appending the libsize to metadata
    feat_norm_for_LR <- NULL
    if (norm_method %in% c("clr", "tss") && "LR" %in% enh_list) {
      feat_norm_for_LR <- DAssemble_normalize_features(
        features          = features,
        metadata          = metadata,
        method            = toupper(norm_method),  # CLR or TSS
        return_normalized = FALSE
      )
    }
    
    for (e in enh_list) {
      fun <- paste0("DA_fit_enhancer_", e)
      
      # Force: LR uses NO CLR, while other enhancers use CLR (or whatever norm_method is)
      feat_use <- feat_norm_default
      if (e == "LR" && identical(norm_method, "clr") && !is.null(feat_norm_for_LR)) {
        feat_use <- feat_norm_for_LR
      }
      
      tmp <- do.call(fun, list(feat_use$features, feat_use$metadata, expVar))
      
      if (!"feature" %in% names(tmp)) {
        stop("Enhancer ", e, " result must contain a 'feature' column.")
      }
      
      pcol_expected <- paste0("pval_", e)
      if (!pcol_expected %in% names(tmp)) {
        stop("Enhancer ", e, " result must contain a '", pcol_expected, "' column.")
      }
      
      if (is.null(res)) {
        res <- tmp
      } else {
        res <- dplyr::left_join(res, tmp, by = c("feature", "metadata"))
      }
      
      components[[e]] <- tmp
    }
  }
  
  ###################################
  # 7) Combine p-values via CCT     #
  ###################################
  
  pv_cols <- grep("^pval_", names(res), value = TRUE)
  # If no core, this does nothing; kept for backwards compatibility
  if (is.null(core_method)) {
    pv_cols <- setdiff(pv_cols, "pval_core")
  }
  if (length(pv_cols) == 0L) {
    stop("No pval_ columns found for combination.")
  }
  
  res$pval_joint <- DA_CCT_rows(as.matrix(res[, pv_cols, drop = FALSE]))
  res$qval       <- stats::p.adjust(res$pval_joint, method = p_adj)
  
  label_parts <- character(0)
  if (!is.null(core_method)) label_parts <- c(label_parts, core_method)
  if (length(enh_list) > 0L) label_parts <- c(label_parts, enh_list)
  method_label <- paste0("CCT(", paste(label_parts, collapse = "+"), ")")
  res <- res[order(res$qval), , drop = FALSE]
  
  ####################
  # 8) Output object #
  ####################
  
  out <- list(
    Method   = method_label,
    res      = res,
    features = features,
    Time.min = (proc.time() - start_time)[["elapsed"]] / 60
  )
  
  if (return_components) {
    out$components <- components
  }
  
  if (return_subensembles && exists("DA_build_subensembles")) {
    res_sub <- res
    # DA_build_subensembles expects Pval_* columns
    names(res_sub) <- sub("^pval_", "Pval_", names(res_sub))
    out$ensembles <- DA_build_subensembles(res_sub, core_method, enh_list, p_adj)
  }
  
  out
}
