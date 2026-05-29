#' Differential analysis using DAssemble
#'
#' The `DAssemble` function implements an ensemble strategy for differential
#' analysis. Users can select a single *core* model and optionally augment it
#' with up to three *enhancer* models. Enhancers are simple, often nonparametric
#' tests that provide complementary evidence. P-values from the core and
#' enhancer models are combined using the Cauchy Combination Test (CCT).
#'
#' @param features A `MultiAssayExperiment` object, or a data.frame with samples
#'   in rows and features in columns.
#' @param metadata A data.frame with sample metadata. Leave as `NULL` when
#'   `features` is a `MultiAssayExperiment`; sample metadata is then taken from
#'   `colData(features)`.
#' @param assay_name Name of the assay/experiment to use when `features` is a
#'   `MultiAssayExperiment`. Required when the object contains multiple
#'   experiments.
#' @param core_method Character scalar naming the core method to run, or `NULL`
#'   / `"none"` for enhancer-only analysis.
#' @param enhancers Character vector containing any of `"WLX"`, `"LR"`, and
#'   `"KS"`.
#' @param expVar Character scalar naming the binary exposure variable in
#'   `metadata` or `colData(features)`.
#' @param p_adj Multiple-testing correction method passed to [stats::p.adjust].
#' @param enhancer_norm Normalization method for enhancer tests. One of
#'   `"TSS"`, `"CLR"`, `"TMM"`, or `"SCRAN"`.
#' @param return_components Logical; if `TRUE`, return per-method component
#'   result tables.
#' @param return_subensembles Logical; if `TRUE`, return CCT sub-ensemble
#'   results for all available method subsets.
#'
#' @return A list containing combined results in `res`, the analyzed feature
#'   table in `features`, elapsed time in `Time.min`, and optionally
#'   `components` and `ensembles`.
#'
#' @examples
#' features <- data.frame(
#'   taxon1 = c(12, 9, 11, 8, 25, 21, 23, 26),
#'   taxon2 = c(30, 28, 31, 29, 10, 12, 11, 9),
#'   taxon3 = c(5, 6, 4, 7, 5, 4, 6, 5),
#'   taxon4 = c(8, 7, 9, 6, 12, 13, 11, 14),
#'   row.names = paste0("sample", seq_len(8))
#' )
#' metadata <- data.frame(
#'   group = factor(rep(c("control", "case"), each = 4)),
#'   row.names = rownames(features)
#' )
#' result <- DAssemble(
#'   features = features,
#'   metadata = metadata,
#'   core_method = NULL,
#'   enhancers = "WLX",
#'   expVar = "group"
#' )
#' head(result$res)
#'
#' @export
DAssemble <- function(features,
                      metadata = NULL,
                      core_method = NULL,
                      enhancers = NULL,
                      expVar = "group",
                      assay_name = NULL,
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

  if (inherits(features, "MultiAssayExperiment")) {
    if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
      stop("The 'MultiAssayExperiment' package is required for MAE input.")
    }
    if (!is.null(metadata)) {
      stop("Do not provide `metadata` when `features` is a MultiAssayExperiment.")
    }
    mae <- features
    assay_names <- names(MultiAssayExperiment::experiments(mae))
    if (length(assay_names) == 0L) {
      stop("The MultiAssayExperiment contains no experiments.")
    }
    if (is.null(assay_name)) {
      if (length(assay_names) != 1L) {
        stop("Specify `assay_name`; the MultiAssayExperiment contains multiple experiments: ",
             paste(assay_names, collapse = ", "))
      }
      assay_name <- assay_names[[1L]]
    }
    if (!assay_name %in% assay_names) {
      stop("assay_name '", assay_name, "' not found in MultiAssayExperiment experiments.")
    }
    
    exp <- MultiAssayExperiment::experiments(mae)[[assay_name]]
    features <- as.data.frame(t(as.matrix(SummarizedExperiment::assay(exp))))
    metadata <- as.data.frame(SummarizedExperiment::colData(mae))
    
    shared_samples <- intersect(rownames(features), rownames(metadata))
    if (length(shared_samples) == 0L) {
      stop("No overlapping sample names between selected assay and MultiAssayExperiment colData.")
    }
    features <- features[shared_samples, , drop = FALSE]
    metadata <- metadata[shared_samples, , drop = FALSE]
  }
  
  if (!is.data.frame(features) || !is.data.frame(metadata)) {
    stop("`features` and `metadata` must both be data.frames, unless `features` is a MultiAssayExperiment.")
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
