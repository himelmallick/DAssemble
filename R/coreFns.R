##################################
# DAssemble Core DESeq2 Pipeline #
##################################

DA_fit_core_DESeq2 <- function(features, metadata, expVar) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("The DESeq2 package is required for core_method = 'DESeq2'.")
  }
  
  # Align samples (rows) between features and metadata
  aligned  <- DA_align_features_metadata(features, metadata)
  features <- aligned$features
  metadata <- aligned$metadata
  
  # Ensure binary factor exposure
  group <- metadata[[expVar]]
  if (is.data.frame(group)) group <- group[[1]]
  group <- droplevels(factor(group))
  if (nlevels(group) != 2L) {
    stop("DESeq2 core requires a binary expVar (2 levels).")
  }
  metadata[[expVar]] <- group
  lvl <- levels(group)
  
  # Fit DESeq2
  design_formula <- stats::as.formula(paste("~", expVar))
  countMat <- t(as.matrix(features))  # DESeq2 expects features x samples
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = countMat,
    colData   = metadata,
    design    = design_formula
  )
  
  dds <- DESeq2::DESeq(dds)
  
  # Use a proper contrast: level2 vs level1
  res <- DESeq2::results(dds, contrast = c(expVar, lvl[2L], lvl[1L]))
  res_df <- as.data.frame(res)
  
  # Add feature column
  res_df$feature <- rownames(res_df)
  
  # Standardize p-value column name
  if (!"pvalue" %in% names(res_df)) {
    stop("DESeq2::results output must contain a 'pvalue' column.")
  }
  res_df$pval_core <- res_df$pvalue
  res_df$pvalue    <- NULL
  
  # Reorder so feature + pval_core come first
  res_df <- res_df[, c("feature", "pval_core",
                       setdiff(names(res_df), c("feature", "pval_core")))]
  
  res_df
}


#########################################
# DAssemble Core metagenomeSeq Pipeline #
#########################################

DA_fit_core_metagenomeSeq <- function(features, metadata, expVar) {
  if (!requireNamespace("metagenomeSeq", quietly = TRUE)) {
    stop("The metagenomeSeq package is required for core_method = 'metagenomeSeq'.")
  }
  
  # Align samples (rows) between features and metadata
  aligned   <- DA_align_features_metadata(features, metadata)
  features  <- aligned$features
  metadata  <- aligned$metadata
  
  # Fit metagenomeSeq 
  design_formula <- stats::as.formula(paste("~", expVar))
  design <- stats::model.matrix(design_formula, data = metadata)
  count_table <- t(as.matrix(features))
  mgsdata <- metagenomeSeq::newMRexperiment(counts = count_table)
  mgsp    <- metagenomeSeq::cumNormStat(mgsdata)
  mgsdata <- metagenomeSeq::cumNorm(mgsdata, mgsp)
  fit <- metagenomeSeq::fitZig(obj = mgsdata, mod = design)
  
  # Gather metagenomeSeq results
  coef_mat <- fit@fit$coefficients
  pval_mat <- fit@eb$p.value
  
  # Drop intercept and any internal scaling column
  drop_cols <- intersect(colnames(coef_mat), c("(Intercept)", "scalingFactor"))
  keep_cols <- setdiff(colnames(coef_mat), drop_cols)
  coef_mat  <- coef_mat[, keep_cols, drop = FALSE]
  pval_mat  <- pval_mat[, keep_cols, drop = FALSE]
  
  # At this point we assume only one column corresponds to expVar
  # (no additional covariates were included).
  coef_vec <- coef_mat[, 1]
  pval_vec <- pval_mat[, 1]
  
  # Construct results table similar to DESeq2 core:
  # - feature column from rownames
  # - pval_core from raw p-values (no adjustment)
  res_df <- data.frame(
    coef   = as.numeric(coef_vec),
    pvalue = as.numeric(pval_vec),
    row.names = rownames(coef_mat),
    stringsAsFactors = FALSE
  )
  
  res_df$feature   <- rownames(res_df)
  res_df$pval_core <- res_df$pvalue
  res_df$pvalue    <- NULL
  
  res_df <- res_df[, c("feature", "pval_core",
                       setdiff(names(res_df), c("feature", "pval_core")))]
  
  res_df
}


#################################
# DAssemble Core edgeR Pipeline #
#################################

DA_fit_core_edgeR <- function(features, metadata, expVar) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("The edgeR package is required for core_method = 'edgeR'.")
  }
  
  # Align samples (rows) between features and metadata
  aligned  <- DA_align_features_metadata(features, metadata)
  features <- aligned$features
  metadata <- aligned$metadata
  
  # Fit edgeR
  design_formula <- stats::as.formula(paste("~", expVar))
  design <- stats::model.matrix(design_formula, data = metadata)
  countMat <- t(as.matrix(features))
  dgelist <- edgeR::DGEList(counts = countMat)
  dgelist <- edgeR::calcNormFactors(dgelist)
  dgelist <- edgeR::estimateDisp(dgelist, design)
  fit <- edgeR::glmFit(dgelist, design)
  lrt <- edgeR::glmLRT(fit, coef = 2)
  res_df <- as.data.frame(edgeR::topTags(lrt, n = nrow(countMat)))
  
  res_df$feature   <- rownames(res_df)
  res_df$pval_core <- res_df$PValue
  res_df$PValue    <- NULL
  
  res_df <- res_df[, c("feature", "pval_core",
                       setdiff(names(res_df), c("feature", "pval_core")))]
  
  res_df
}


##################################
# DAssemble Core Robseq Pipeline #
##################################

DA_fit_core_Robseq <- function(features, metadata, expVar) {
  if (!requireNamespace("Robseq", quietly = TRUE)) {
    stop("The Robseq package is required for core_method = 'Robseq'.")
  }
  
  # Align rownames / sample order (samples x features)
  aligned  <- DA_align_features_metadata(features, metadata)
  features <- aligned$features
  metadata <- aligned$metadata
  
  # Robseq expects genes in rows, samples in columns
  countMat <- t(as.matrix(features))
  
  # --- Prepare metadata for Robseq ---
  # Robseq's default expVar is "Exposure"; if we don't pass expVar,
  # it will look for metadata$Exposure. Make sure that exists.
  meta_rob <- as.data.frame(metadata)
  
  if (!expVar %in% names(meta_rob)) {
    stop("For Robseq core, metadata must contain column '", expVar, "'.")
  }
  
  # Guarantee an 'Exposure' column that Robseq can use by default
  if (!"Exposure" %in% names(meta_rob)) {
    meta_rob$Exposure <- meta_rob[[expVar]]
  }
  
  # Call robust.dge letting it use the default expVar = "Exposure"
  fit <- Robseq::robust.dge(
    features    = countMat,
    metadata    = meta_rob,
    norm.method = "RLE",
    # expVar omitted -> uses default "Exposure"
    coVars      = NULL,
    filter      = FALSE,
    parallel    = FALSE,
    ncores      = 1,
    verbose     = FALSE
  )
  
  if (!"res" %in% names(fit)) {
    stop("Robseq::robust.dge output does not contain a 'res' element.")
  }
  
  res_df <- as.data.frame(fit$res)
  
  # Expect columns: Genes, Log2FC, SE, LCI, UCI, Pval, adjPval
  if (!"Pval" %in% names(res_df)) {
    stop("Robseq results must contain a 'Pval' column.")
  }
  
  # Standardize to feature + pval_core for DAssemble
  if ("Genes" %in% names(res_df)) {
    res_df$feature <- res_df$Genes
  } else if (!is.null(rownames(res_df))) {
    res_df$feature <- rownames(res_df)
  } else {
    stop("Robseq results do not contain 'Genes' column and rownames are NULL.")
  }
  
  res_df$pval_core <- res_df$Pval
  
  # Drop originals if you don't want duplication
  res_df$Pval <- NULL
  if ("Genes" %in% names(res_df)) {
    res_df$Genes <- NULL
  }
  
  # Put feature + pval_core first
  res_df <- res_df[, c("feature", "pval_core",
                       setdiff(names(res_df), c("feature", "pval_core")))]
  
  res_df
}

#######################################
# DAssemble Core limma-voom Pipeline  #
#######################################

DA_fit_core_limmaVOOM <- function(features, metadata, expVar) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("The limma package is required for core_method = 'limmaVOOM'.")
  }
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("The edgeR package is required for core_method = 'limmaVOOM'.")
  }
  
  # Align samples (rows) between features and metadata
  aligned  <- DA_align_features_metadata(features, metadata)
  features <- aligned$features
  metadata <- aligned$metadata
  
  # Minimal design: intercept + expVar
  design_formula <- stats::as.formula(paste("~", expVar))
  design <- stats::model.matrix(design_formula, data = metadata)
  
  # Counts: samples x features
  countMat <- t(as.matrix(features))
  
  # Voom transform
  y <- edgeR::DGEList(counts = countMat)
  y <- edgeR::calcNormFactors(y)
  v <- limma::voom(y, design, plot = FALSE)
  
  fit <- limma::lmFit(v, design)
  fit <- limma::eBayes(fit)
  
  # Extract topTable for the coefficient of expVar
  tt <- limma::topTable(fit, coef = 2, number = nrow(countMat), sort.by = "none")
  
  res_df <- as.data.frame(tt)
  res_df$feature   <- rownames(res_df)
  res_df$pval_core <- res_df$P.Value
  
  res_df$P.Value <- NULL
  
  res_df <- res_df[, c("feature", "pval_core",
                       setdiff(names(res_df), c("feature", "pval_core")))]
  
  res_df
}


###################################
# DAssemble Core dearseq Pipeline #
###################################

DA_fit_core_dearseq <- function(features, metadata, expVar) {
  if (!requireNamespace("dearseq", quietly = TRUE)) {
    stop("The dearseq package is required for core_method = 'dearseq'.")
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("The SummarizedExperiment package is required for core_method = 'dearseq'.")
  }
  
  # Align samples (rows) between features and metadata
  aligned  <- DA_align_features_metadata(features, metadata)
  features <- aligned$features
  metadata <- aligned$metadata
  
  countMat <- t(as.matrix(features))
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = countMat),
    colData = metadata
  )
  
  which_test <- "asymptotic"
  fit <- dearseq::dear_seq(
    object         = se,
    variables2test = expVar,
    which_test     = which_test,
    preprocessed   = FALSE,
    progressbar    = FALSE,
    parallel_comp  = FALSE,
    which_weights  = "loclin"
  )
  
  # Raw p-values per feature
  pvals <- fit$pvals$rawPval
  genes <- rownames(fit$precision_weights)
  
  # Gather dearseq
  res_df <- data.frame(
    feature   = genes,
    pval_core = as.numeric(pvals),
    stringsAsFactors = FALSE
  )
  
  res_df
}


###############################
# DAssemble Core MAST Pipeline
###############################

DA_fit_core_MAST <- function(features, metadata, expVar) {
  if (!requireNamespace("MAST", quietly = TRUE)) {
    stop("The MAST package is required for core_method = 'MAST'.")
  }
  
  # Align samples (rows) between features and metadata
  aligned  <- DA_align_features_metadata(features, metadata)
  features <- aligned$features
  metadata <- aligned$metadata
  
  # Convert to SingleCellAssay
  sca <- MAST::FromMatrix(
    exprsArray = t(as.matrix(features)),
    cData      = metadata
  )
  
  # Fit hurdle model: \code{zlm(~ expVar, sca)}
  zlmFit <- MAST::zlm(
    stats::as.formula(paste("~", expVar)),
    sca
  )
  
  # Test the coefficient expVar
  summary_zlm <- MAST::summary(zlmFit, doLRT = expVar)
  
  # "Dat" slot includes per-feature statistics; we focus on the Hurdle test
  result <- summary_zlm$datatable
  res_hurdle <- result[result$contrast == expVar & result$component == "H", ]
  
  res_df <- as.data.frame(res_hurdle)
  
  # Ensure feature identifier
  if (!"primerid" %in% names(res_df)) {
    stop("MAST result does not contain 'primerid' column.")
  }
  res_df$feature <- res_df$primerid
  
  # Map p-value to pval_core
  if (!"Pr(>Chisq)" %in% names(res_df)) {
    stop("MAST result does not contain 'Pr(>Chisq)' column.")
  }
  res_df$pval_core <- res_df[["Pr(>Chisq)"]]
  
  # Drop original columns if desired
  res_df$primerid        <- NULL
  res_df[["Pr(>Chisq)"]] <- NULL
  
  res_df <- res_df[, c("feature", "pval_core",
                       setdiff(names(res_df), c("feature", "pval_core")))]
  
  res_df
}


###################################
# DAssemble Core ANCOMBC2        #
###################################

DA_fit_core_ANCOMBC2 <- function(features, metadata, expVar) {
  if (!requireNamespace("ANCOMBC", quietly = TRUE) ||
      !requireNamespace("TreeSummarizedExperiment", quietly = TRUE) ||
      !requireNamespace("S4Vectors", quietly = TRUE)) {
    stop("ANCOMBC, TreeSummarizedExperiment, and S4Vectors are required for core_method = 'ANCOMBC2'.")
  }
  
  # Align samples (rows) between features and metadata
  aligned  <- DA_align_features_metadata(features, metadata)
  features <- aligned$features
  metadata <- aligned$metadata
  
  assays <- S4Vectors::SimpleList(counts = t(as.matrix(features)))
  tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays  = assays,
    colData = metadata
  )
  
  res <- ANCOMBC::ancombc2(
    data   = tse,
    assay_name = "counts",
    tax_level  = NULL,
    fix_formula = stats::as.formula(paste("~", expVar)),
    rand_formula = NULL,
    p_adj_method = "BH",
    prv_cut   = 0,
    lib_cut   = 0,
    s0_perc   = 0.05,
    alpha     = 0.05,
    n_cl      = 1,
    global    = FALSE
  )
  
  # res$res contains per-feature statistics; we look for p-values for expVar
  res_df <- as.data.frame(res$res)
  
  # We expect at least one column with 'p_' or 'pval' referencing expVar
  p_candidates <- grep(paste0("(^p_|pval).*", expVar), names(res_df), value = TRUE)
  if (!length(p_candidates)) {
    # fallback: any column starting with "p_"
    p_candidates <- grep("^p_", names(res_df), value = TRUE)
  }
  if (!length(p_candidates)) {
    stop("ANCOMBC2 results do not contain a p-value column for expVar.")
  }
  pcol <- p_candidates[1L]
  
  res_df$feature   <- res_df$taxon
  res_df$pval_core <- res_df[[pcol]]
  
  res_df[[pcol]] <- NULL
  res_df$taxon   <- NULL
  
  res_df <- res_df[, c("feature", "pval_core",
                       setdiff(names(res_df), c("feature", "pval_core")))]
  res_df
}


#################################
# DAssemble Core ALDEx2 Pipeline
#################################

DA_fit_core_ALDEx2 <- function(features, metadata, expVar) {
  if (!requireNamespace("ALDEx2", quietly = TRUE)) {
    stop("The ALDEx2 package is required for core_method = 'ALDEx2'.")
  }
  
  # Align samples (rows) between features and metadata
  aligned  <- DA_align_features_metadata(features, metadata)
  features <- aligned$features
  metadata <- aligned$metadata
  
  # ALDEx2 expects features in rows and samples in columns
  countMat <- t(as.matrix(features))
  
  group <- metadata[[expVar]]
  if (is.data.frame(group)) group <- group[[1]]
  group <- droplevels(factor(group))
  if (nlevels(group) != 2L) {
    stop("ALDEx2 core requires a binary expVar.")
  }
  
  condition <- group
  aldex_obj <- ALDEx2::aldex(
    reads   = countMat,
    conds   = condition,
    test    = "t",
    effect  = TRUE,
    mc.samples = 8,
    denom   = "all"
  )
  
  res_df <- as.data.frame(aldex_obj)
  res_df$feature   <- rownames(res_df)
  
  # ALDEx2 has several p-values; we choose the Welch's t-test p-value (we.ep)
  if (!"we.ep" %in% names(res_df)) {
    stop("ALDEx2 result must contain 'we.ep' column.")
  }
  res_df$pval_core <- res_df$we.ep
  
  res_df <- res_df[, c("feature", "pval_core",
                       setdiff(names(res_df), c("feature", "pval_core")))]
  
  res_df
}


#################################
# DAssemble Core LinDA Pipeline #
#################################

DA_fit_core_LinDA <- function(features, metadata, expVar) {
  if (!requireNamespace("LinDA", quietly = TRUE)) {
    stop("The LinDA package is required for core_method = 'LinDA'.")
  }
  
  # Align samples (rows) between features and metadata
  aligned  <- DA_align_features_metadata(features, metadata)
  features <- aligned$features
  metadata <- aligned$metadata
  
  # LinDA expects otu_table features x samples
  otu_table <- t(as.matrix(features))
  meta      <- metadata
  
  formula_str <- paste("~", expVar)
  
  fit <- LinDA::linda(
    otu_tab  = otu_table,
    meta_tab = meta,
    formula  = formula_str,
    alpha    = 0.05
  )
  
  # Extract results for expVar
  if (!expVar %in% names(fit)) {
    stop("LinDA result does not contain a component for expVar '", expVar, "'.")
  }
  
  res_df <- as.data.frame(fit[[expVar]])
  
  # Ensure feature column
  if (!"feature" %in% names(res_df)) {
    if (!is.null(rownames(res_df))) {
      res_df$feature <- rownames(res_df)
    } else {
      stop("LinDA results do not contain 'feature' column and rownames are NULL.")
    }
  }
  
  # LinDA raw p-value is 'pvalue'
  if (!"pvalue" %in% names(res_df)) {
    stop("LinDA result must contain 'pvalue' column.")
  }
  res_df$pval_core <- res_df$pvalue
  res_df$pvalue    <- NULL
  
  res_df <- res_df[, c("feature", "pval_core",
                       setdiff(names(res_df), c("feature", "pval_core")))]
  
  res_df
}


#################################
# DAssemble Core LOCOM Pipeline #
#################################

DA_fit_core_LOCOM <- function(features, metadata, expVar) {
  if (!requireNamespace("LOCOM", quietly = TRUE)) {
    stop("The LOCOM package is required for core_method = 'LOCOM'.")
  }
  
  # Align samples (rows) between features and metadata
  aligned  <- DA_align_features_metadata(features, metadata)
  features <- aligned$features
  metadata <- aligned$metadata
  
  otu.table <- as.matrix(features)
  Y <- metadata[[expVar]]
  if (is.factor(Y)) Y <- as.numeric(Y) - 1
  if (!is.numeric(Y)) Y <- as.numeric(Y)
  
  fit <- LOCOM::locom(
    otu.table     = otu.table,
    Y             = Y,
    fdr.nominal   = 0.1,
    seed          = 1,
    filter.thresh = 0,
    n.perm.max    = 1000
  )
  
  # fit$feature.table should contain per-feature p-values
  res_df <- as.data.frame(fit$feature.table)
  
  # Try to find p-value column
  p_candidates <- grep("pval|p.value|p_value", names(res_df), ignore.case = TRUE, value = TRUE)
  if (!length(p_candidates)) {
    stop("LOCOM feature.table does not contain a recognizable p-value column.")
  }
  pcol <- p_candidates[1L]
  
  if (!"feature" %in% names(res_df)) {
    if (!is.null(rownames(res_df))) {
      res_df$feature <- rownames(res_df)
    } else {
      stop("LOCOM results do not contain 'feature' column and rownames are NULL.")
    }
  }
  
  res_df$pval_core <- res_df[[pcol]]
  res_df[[pcol]]   <- NULL
  
  res_df <- res_df[, c("feature", "pval_core",
                       setdiff(names(res_df), c("feature", "pval_core")))]
  
  res_df
}


###################################
# DAssemble Core MaAsLin2 Pipeline
###################################

DA_fit_core_Maaslin2 <- function(features, metadata, expVar) {
  if (!requireNamespace("Maaslin2", quietly = TRUE)) {
    stop("The Maaslin2 package is required for core_method = 'Maaslin2'.")
  }
  
  # Align samples (rows) between features and metadata
  aligned  <- DA_align_features_metadata(features, metadata)
  features <- aligned$features
  metadata <- aligned$metadata
  
  # Fit MaAsLin2
  tmp <- file.path(tempdir(), paste0("maaslin2_", as.integer(runif(1, 1e6, 1e7))))
  dir.create(tmp, showWarnings = FALSE)
  fit_out <- Maaslin2::Maaslin2(
    input_data    = t(as.matrix(features)),
    input_metadata = metadata,
    output        = tmp,
    fixed_effects = expVar
  )
  
  res_df <- fit_out$results
  
  unlink(tmp, recursive = TRUE)
  
  # Results must contain:
  #   feature, metadata, pval, qval, etc.
  if (!all(c("feature", "metadata") %in% names(res_df))) {
    stop("Maaslin2 results must contain 'feature' and 'metadata' columns.")
  }
  
  res_df <- res_df[res_df$metadata == expVar, , drop = FALSE]
  
  if (!"pval" %in% names(res_df)) {
    stop("Maaslin2 results must contain 'pval' column.")
  }
  
  res_df$pval_core <- res_df$pval
  res_df$pval      <- NULL
  
  # Maaslin2 already has a "feature" column
  res_df <- res_df[, c("feature", "pval_core",
                       setdiff(names(res_df), c("feature", "pval_core")))]
  res_df
}


####################################
# DAssemble Core MaAsLin3 Pipeline #
####################################

DA_fit_core_Maaslin3 <- function(features, metadata, expVar) {
  if (!requireNamespace("maaslin3", quietly = TRUE)) {
    stop("maaslin3 not installed.")
  }
  
  # Align samples (rows) between features and metadata
  aligned  <- DA_align_features_metadata(features, metadata)
  features <- aligned$features      # samples x features
  metadata <- aligned$metadata
  
  # Temporary output directory for maaslin3 side-effects
  tmp <- file.path(tempdir(), paste0("maaslin3_", as.integer(runif(1, 1e6, 1e7))))
  dir.create(tmp, showWarnings = FALSE)
  
  # Call maaslin3 with a minimal, broadly compatible set of arguments
  fit <- maaslin3::maaslin3(
    input_data     = t(as.matrix(features)),  # features x samples
    input_metadata = metadata,
    output         = tmp,
    fixed_effects  = expVar
  )
  
  # Clean up temporary directory
  unlink(tmp, recursive = TRUE)
  
  # MaAsLin3 returns a list; we pull the abundance results table
  res_ab <- fit$fit_data_abundance$results
  if (is.null(res_ab)) {
    stop("Maaslin3: fit_data_abundance$results is NULL; no abundance results produced.")
  }
  
  res_df <- as.data.frame(res_ab)
  
  # Filter to the exposure variable of interest
  if (!"metadata" %in% names(res_df)) {
    stop("Maaslin3 results must contain a 'metadata' column.")
  }
  res_df <- res_df[res_df$metadata == expVar, , drop = FALSE]
  
  # Ensure feature column
  if (!"feature" %in% names(res_df)) {
    if (!is.null(rownames(res_df))) {
      res_df$feature <- rownames(res_df)
    } else {
      stop("Maaslin3 results do not contain 'feature' column and rownames are NULL.")
    }
  }
  
  # Use the uncorrected per-association p-value as the core p-value
  if (!"pval_individual" %in% names(res_df)) {
    stop("Maaslin3 results must contain a 'pval_individual' column.")
  }
  res_df$pval_core <- res_df$pval_joint
  
  # ---- Avoid name collisions with DAssemble ----
  # 1) Drop the original pval_individual so only pval_core matches '^pval_'
  res_df$pval_individual <- NULL
  
  # 2) Rename Maaslin3's joint p-values so they don't override DAssemble's pval_joint
  if ("pval_joint" %in% names(res_df)) {
    names(res_df)[names(res_df) == "pval_joint"] <- "maaslin3_pval_joint"
  }
  if ("qval_joint" %in% names(res_df)) {
    names(res_df)[names(res_df) == "qval_joint"] <- "maaslin3_qval_joint"
  }
  
  # Reorder so feature + pval_core come first; keep everything else
  res_df <- res_df[, c("feature", "pval_core",
                       setdiff(names(res_df), c("feature", "pval_core")))]
  
  res_df
}


###################################
# DAssemble Core Tweedieverse     #
###################################

DA_fit_core_Tweedieverse <- function(features, metadata, expVar,
                                     transf = "NONE",
                                     multiple_qvalues = FALSE) {
  
  if (!requireNamespace("Tweedieverse", quietly = TRUE)) {
    stop("The Tweedieverse package is required for core_method = 'Tweedieverse'.")
  }
  
  ## 1) Align & coerce
  aligned  <- DA_align_features_metadata(features, metadata)
  features <- aligned$features
  metadata <- aligned$metadata
  
  ## Tweedieverse expects features x samples
  X <- if (nrow(features) == nrow(metadata)) t(features) else features
  X <- as.data.frame(as.matrix(X), check.names = FALSE)
  md <- as.data.frame(metadata, check.names = FALSE)
  
  ## Optional scran normalization
  scale_factor_arg <- NULL
  if (transf == "scran") {
    if (!requireNamespace("scran", quietly = TRUE) ||
        !requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      stop("scran normalization requires scran + SingleCellExperiment.")
    }
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = t(X)))
    clusters <- scran::quickCluster(sce)
    sce <- scran::computeSumFactors(sce, clusters = clusters)
    md$scranSize <- scran::sizeFactors(sce)
    scale_factor_arg <- "scranSize"
  }
  
  # Temporary output directory for maaslin3 side-effects
  tmp <- file.path(tempdir(), paste0("Tweedieverse_", as.integer(runif(1, 1e6, 1e7))))
  dir.create(tmp, showWarnings = FALSE)
  
  ## 2) Tweedieverse call
    # Temporary output directory for maaslin3 side-effects
  tmp <- file.path(tempdir(), paste0("maaslin3_", as.integer(runif(1, 1e6, 1e7))))
  dir.create(tmp, showWarnings = FALSE)
  tw <- Tweedieverse::Tweedieverse(
    input_features        = X,
    input_metadata        = md,
    output = tmp,
    correction      = "BH",
    max_significance = 1,
    prev_threshold  = 0,
    fixed_effects   = expVar,
    base_model      = "CPLM",
    adjust_offset   = TRUE,
    scale_factor    = scale_factor_arg,
    #median_comparison   = FALSE,
    #median_subtraction  = FALSE,
    na.action       = stats::na.pass
  )
  
  ## 3) Keep only main columns (just like DAssemble_Revision)
  df <- as.data.frame(tw)
  drop_cols <- c("value","qval","base.model","tweedie.index",
                 "N","N.not.zero","percent.zero")
  drop_cols <- drop_cols[drop_cols %in% names(df)]
  
  df <- df[, setdiff(names(df), drop_cols), drop = FALSE]
  
  ## 4) Standardize feature column
  if (!"feature" %in% names(df)) {
    if (!is.null(rownames(df))) {
      df$feature <- rownames(df)
    } else {
      stop("Tweedieverse results lack a 'feature' column.")
    }
  }
  
  ## 5) Extract p-value
  p_candidates <- grep("pval|p.value", names(df), ignore.case = TRUE, value = TRUE)
  if (!length(p_candidates)) {
    stop("Tweedieverse results contain no p-value column.")
  }
  pcol <- p_candidates[1]
  
  df$pval_core <- df[[pcol]]
  df[[pcol]]   <- NULL
  
  ## 6) Additional q-values (optional)
  if (multiple_qvalues) {
    stop("multiple_qvalues not implemented in this package version.")
  } else {
    df$adjPval <- p.adjust(df$pval_core, method = "BY")
  }
  
  ## 7) Reorder columns
  df <- df[order(df$adjPval), ]
  df <- df[, c("feature", "pval_core",
               setdiff(names(df), c("feature", "pval_core")))]
  
  df
}