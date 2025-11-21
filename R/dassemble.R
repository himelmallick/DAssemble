#' Differential analysis using DAssemble
#'
#' The `DAssemble` function implements an ensemble strategy for differential
#' analysis. Users can select a single *core* model and optionally augment it
#' with up to two *enhancer* models. Enhancers are simple, often nonparametric
#' tests that provide complementary evidence. P-values from the core and
#' enhancer models are combined using the Cauchy Combination Test (CCT).
#'
#' @param features Matrix or data.frame of features. Either
#'   features x samples or samples x features (see `orientation`).
#' @param metadata Data.frame of sample-level metadata. Row names or a
#'   column named "sample" should match sample IDs in `features`.
#' @param core_method Character: name of core method, one of
#'   "DESeq2", "edgeR", "limmaVOOM", "dearseq",
#'   "metagenomeSeq", "MAST", "Tweedieverse",
#'   "Maaslin2", "Maaslin3", "LOCOM",
#'   "LinDA", "ANCOMBC2", "Robseq", "ALDEx2",
#'   or NULL / "none" to omit the core.
#' @param enhancers Character vector of up to two enhancers chosen from
#'   "WLX", "LR", "KS". Use NULL for no enhancers.
#' @param expVar Character: name of primary exposure variable in `metadata`.
#'   Currently assumed binary for all methods.
#' @param coVars Optional character vector of covariate names in `metadata`
#'   (used by some core methods).
#' @param p_adj Multiple testing correction method passed to stats::p.adjust.
#' @param orientation "auto", "features_by_samples", or "samples_by_features".
#' @param contrast Optional contrast specification for certain core methods.
#' @param quiet Logical; if TRUE, core methods should be as quiet as possible.
#' @param domain "none","bulkrnaseq","singlecell","microbiome"; used by
#'   enhancer normalization.
#' @param return_components Logical; return per-method results in $components.
#' @param return_subensembles Logical; if TRUE, attempt to construct CCT-based
#'   subensembles via DA_build_subensembles (if available).
#'
#' @return A list with elements:
#'   \describe{
#'     \item{Method}{Label describing the ensemble, e.g. "CCT(DESeq2+WLX)".}
#'     \item{res}{Data.frame with columns including:
#'       `feature`, `pval_core` (if core used),
#'       `pval_WLX` / `pval_LR` / `pval_KS` (if enhancers used),
#'       `pval_joint` (CCT p-value), and `qval` (adjusted joint p-value).}
#'     \item{features}{Original features input.}
#'     \item{Time.min}{Elapsed time in minutes.}
#'     \item{components}{(Optional) list of per-method result tables.}
#'     \item{ensembles}{(Optional) sub-ensemble results if requested.}
#'   }
#'
#' @export
DAssemble <- function(features,
                      metadata,
                      core_method         = NULL,
                      enhancers           = NULL,
                      expVar,
                      coVars              = NULL,
                      p_adj               = "BH",
                      orientation         = c("auto",
                                              "features_by_samples",
                                              "samples_by_features"),
                      contrast            = NULL,
                      quiet               = TRUE,
                      domain              = c("none",
                                              "bulkrnaseq",
                                              "singlecell",
                                              "microbiome"),
                      return_components   = TRUE,
                      return_subensembles = FALSE) {
  
  start_time  <- proc.time()
  orientation <- match.arg(orientation)
  domain      <- match.arg(domain)
  
  ## 1) Validate core + enhancers -----------------------------------------
  
  valid_core <- c(
    "DESeq2", "edgeR", "limmaVOOM", "dearseq",
    "metagenomeSeq", "MAST", "Tweedieverse",
    "Maaslin2", "Maaslin3", "LOCOM",
    "LinDA", "ANCOMBC2", "Robseq", "ALDEx2"
  )
  
  if (is.null(core_method) || identical(core_method, "none")) {
    core_method <- NULL
  } else if (!core_method %in% valid_core) {
    stop("Unknown core_method: ", core_method,
         ". Supported core methods: ",
         paste(valid_core, collapse = ", "),
         " (or use NULL / 'none').")
  }
  
  if (is.null(enhancers)) {
    enh_list <- character()
  } else {
    enh_list <- unique(enhancers)
  }
  if (length(enh_list) > 2L) {
    stop("At most two enhancers are supported; you supplied: ",
         paste(enh_list, collapse = ", "))
  }
  
  valid_enh <- c("WLX", "LR", "KS")
  if (length(enh_list) > 0L && !all(enh_list %in% valid_enh)) {
    bad <- setdiff(enh_list, valid_enh)
    stop("Unknown enhancer(s): ", paste(bad, collapse = ", "),
         ". Supported enhancers: ", paste(valid_enh, collapse = ", "))
  }
  
  if (is.null(core_method) && length(enh_list) == 0L) {
    stop("You must specify either a core_method or at least one enhancer.")
  }
  
  ## 2) Align features and metadata (always get samples x features) -------
  
  aligned <- DA_align_counts_and_md(
    X           = features,
    md          = metadata,
    orientation = orientation,
    expVar      = expVar,
    coVars      = coVars
  )
  counts <- aligned$counts  # samples x features
  md     <- aligned$md
  
  if (!expVar %in% colnames(md)) {
    stop("expVar '", expVar, "' not found in metadata.")
  }
  
  grp <- md[[expVar]]
  if (!is.factor(grp)) grp <- factor(grp)
  grp <- droplevels(grp)
  if (nlevels(grp) != 2L) {
    stop("For now, DAssemble requires a binary expVar (2 levels).")
  }
  md[[expVar]] <- grp
  
  ## 3) Core model --------------------------------------------------------
  
  components <- list()
  
  if (!is.null(core_method)) {
    core_fun <- paste0("DA_fit_core_", core_method)
    if (!exists(core_fun, mode = "function")) {
      stop("Core method '", core_method,
           "' is not implemented. Please define ", core_fun, "().")
    }
    
    res_core <- DA_do_call(core_fun, list(
      features = counts,
      metadata = md,
      expVar   = expVar,
      coVars   = coVars,
      contrast = contrast,
      quiet    = quiet
    ))
    
    if (!("feature" %in% names(res_core))) {
      if ("Genes" %in% names(res_core)) {
        res_core$feature <- res_core$Genes
      } else if (!is.null(rownames(res_core))) {
        res_core$feature <- rownames(res_core)
      } else {
        stop("Core result must contain or be mappable to a 'feature' column.")
      }
    }
    
    if (!("pval_core" %in% names(res_core))) {
      if ("pvalue" %in% names(res_core)) {
        res_core$pval_core <- res_core$pvalue
        res_core$pvalue    <- NULL
      } else if ("Pval" %in% names(res_core)) {
        res_core$pval_core <- res_core$Pval
        res_core$Pval      <- NULL
      } else {
        stop("Core result must contain a 'pval_core' (or 'pvalue'/ 'Pval') column.")
      }
    }
    
    res_core <- res_core[, c("feature", "pval_core",
                             setdiff(names(res_core), c("feature", "pval_core")))]
    
    res <- res_core
    components$core <- res_core
  } else {
    genes <- colnames(counts)
    if (is.null(genes)) genes <- paste0("feat_", seq_len(ncol(counts)))
    res <- data.frame(
      feature   = genes,
      pval_core = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  
  ## 4) Enhancers ---------------------------------------------------------
  
  if (length(enh_list) > 0L) {
    counts_enh <- DA_norm_for_enhancers(counts, domain = domain)
    
    for (e in enh_list) {
      if (identical(e, "WLX")) {
        w <- DA_fit_enhancer_WLX(counts_enh, md, expVar, domain = domain)
        res <- dplyr::left_join(res, w, by = "feature")
        components$WLX <- w
        
      } else if (identical(e, "LR")) {
        lr <- DA_fit_enhancer_LR(counts_enh, md, expVar, domain = domain)
        res <- dplyr::left_join(res, lr, by = "feature")
        components$LR <- lr
        
      } else if (identical(e, "KS")) {
        ks <- DA_fit_enhancer_KS(counts_enh, md, expVar, domain = domain)
        res <- dplyr::left_join(res, ks, by = "feature")
        components$KS <- ks
      }
    }
  }
  
  ## 5) Combine p-values via CCT ------------------------------------------
  
  pv_cols <- grep("^pval_", names(res), value = TRUE)
  
  if (is.null(core_method)) {
    pv_cols <- setdiff(pv_cols, "pval_core")
  }
  
  if (length(pv_cols) == 0L) {
    stop("No per-method p-value columns found to combine (pval_*).")
  }
  
  pv_mat <- as.matrix(res[, pv_cols, drop = FALSE])
  res$pval_joint <- as.numeric(DA_CCT_rows(pv_mat))
  res$qval       <- as.numeric(stats::p.adjust(res$pval_joint, method = p_adj))
  
  if (!is.null(core_method)) {
    label_parts <- c(core_method, enh_list)
  } else {
    label_parts <- enh_list
  }
  method_label <- paste0("CCT(", paste(label_parts, collapse = "+"), ")")
  
  if ("pval_joint" %in% names(res)) {
    ord <- order(res$pval_joint, na.last = TRUE)
    res <- res[ord, , drop = FALSE]
  }
  
  ## 6) Output object -----------------------------------------------------
  
  elapsed <- proc.time() - start_time
  out <- list(
    Method   = method_label,
    res      = res,
    features = features,
    Time.min = as.numeric(elapsed["elapsed"]) / 60
  )
  
  if (return_components) {
    out$components <- components
  }
  
  if (return_subensembles) {
    if (!exists("DA_build_subensembles", mode = "function")) {
      warning("return_subensembles = TRUE, but DA_build_subensembles() is not defined; ",
              "no subensembles will be returned.")
    } else {
      # Temporary mapping for existing DA_build_subensembles that expects Pval_*:
      res_sub <- res
      if ("pval_core" %in% names(res_sub)) res_sub$Pval_core <- res_sub$pval_core
      if ("pval_WLX"  %in% names(res_sub)) res_sub$Pval_WLX  <- res_sub$pval_WLX
      if ("pval_LR"   %in% names(res_sub)) res_sub$Pval_LR   <- res_sub$pval_LR
      if ("pval_KS"   %in% names(res_sub)) res_sub$Pval_KS   <- res_sub$pval_KS
      
      out$ensembles <- DA_build_subensembles(
        res         = res_sub,
        core_method = core_method,
        enh_list    = enh_list,
        p_adj       = p_adj
      )
    }
  }
  
  out
}