#########################
# DAssemble Core DESeq2 #
#########################

DA_fit_core_DESeq2 <- function(features, metadata, expVar, coVars = NULL, ...) {
  
  ########################
  # Package sanity check #
  ########################
  
  if (!requireNamespace("DESeq2", quietly = TRUE))
    stop("DESeq2 required.")
  
  ############################
  # Standard DESeq2 pipeline #
  ############################
  
  design <- stats::as.formula(paste("~", build_rhs(expVar, coVars)))
  extra_args <- list(...)
  dds_args <- extra_args$DESeqDataSetFromMatrix %||% list()
  sizefactor_args <- extra_args$estimateSizeFactors %||% list()
  deseq_args <- extra_args$DESeq %||% list()
  results_args <- extra_args$results %||% list()

  x <- do.call(
    DESeq2::DESeqDataSetFromMatrix,
    c(list(countData = t(as.matrix(features)), colData = metadata, design = design), dds_args)
  )
  geoMeans <- apply(
    DESeq2::counts(x),
    1,
    DA_positive_geometric_mean
  )
  x <- do.call(
    DESeq2::estimateSizeFactors,
    c(list(object = x, geoMeans = geoMeans), sizefactor_args)
  )
  fit <- do.call(DESeq2::DESeq, c(list(object = x), deseq_args))
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  feature <- rownames(stats::coef(fit))
  result_names <- DESeq2::resultsNames(fit)
  coef_name <- grep(paste0("^", expVar), result_names, value = TRUE)[1]
  if (is.na(coef_name)) {
    stop("Could not identify the exposure coefficient in DESeq2 output.")
  }
  pval_core <- do.call(
    DESeq2::results,
    c(list(object = fit, name = coef_name), results_args)
  )$pvalue
  return(DA_format_result(feature, expVar, pval_core = pval_core))
}

########################
# DAssemble Core edgeR #
########################

DA_fit_core_edgeR <- function(features, metadata, expVar, coVars = NULL, ...) {
  
  ########################
  # Package sanity check #
  ########################
  
  if (!requireNamespace("edgeR", quietly = TRUE))
    stop("edgeR required.")
  
  ###########################
  # Standard edgeR pipeline #
  ###########################
  
  extra_args <- list(...)
  dgelist_args <- extra_args$DGEList %||% list()
  norm_args <- extra_args$calcNormFactors %||% list()
  commondisp_args <- extra_args$estimateGLMCommonDisp %||% list()
  trendeddisp_args <- extra_args$estimateGLMTrendedDisp %||% list()
  tagwisedisp_args <- extra_args$estimateGLMTagwiseDisp %||% list()
  glmfit_args <- extra_args$glmFit %||% list()
  glmlrt_args <- extra_args$glmLRT %||% list()

  d <- do.call(edgeR::DGEList, c(list(counts = t(features)), dgelist_args))
  d <- do.call(edgeR::calcNormFactors, c(list(object = d, method = "TMM"), norm_args))
  design <- stats::model.matrix(stats::as.formula(paste("~", build_rhs(expVar, coVars))), metadata)
  d <- do.call(edgeR::estimateGLMCommonDisp, c(list(y = d, design = design), commondisp_args))
  d <- do.call(edgeR::estimateGLMTrendedDisp, c(list(y = d, design = design), trendeddisp_args))
  d <- do.call(edgeR::estimateGLMTagwiseDisp, c(list(y = d, design = design), tagwisedisp_args))
  fit <- do.call(edgeR::glmFit, c(list(y = d, design = design), glmfit_args))
  coef_name <- get_exp_coef_name(metadata, expVar, coVars)
  fit <- do.call(
    edgeR::glmLRT,
    c(list(glmfit = fit, coef = which(colnames(design) == coef_name)), glmlrt_args)
  )
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  feature<-rownames(fit$table)
  pval_core<-fit$table$PValue
  return(DA_format_result(feature, expVar, pval_core = pval_core))
  
}

#############################
# DAssemble Core limmaVOOM  #
#############################

DA_fit_core_limmaVOOM <- function(features, metadata, expVar, coVars = NULL, ...) {
  
  ########################
  # Package sanity check #
  ########################
  
  if (!requireNamespace("limma", quietly = TRUE) ||
      !requireNamespace("edgeR", quietly = TRUE))
    stop("limma + edgeR required.")
  
  ###############################
  # Standard limmaVOOM pipeline #
  ###############################
  
  extra_args <- list(...)
  voom_args <- extra_args$voom %||% list()
  lmfit_args <- extra_args$lmFit %||% list()
  ebayes_args <- extra_args$eBayes %||% list()

  design <- stats::model.matrix(stats::as.formula(paste("~", build_rhs(expVar, coVars))), metadata)
  x<-t(as.matrix(features)+1) # Convert to matrix, round up to nearest integer, and transpose
  y <- do.call(limma::voom, c(list(counts = x, design = design, plot = FALSE), voom_args))
  fit <- do.call(limma::lmFit, c(list(object = y, design = design), lmfit_args))
  fit <- do.call(limma::eBayes, c(list(fit = fit), ebayes_args))
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  feature<-rownames(fit$coefficients)
  coef_name <- get_exp_coef_name(metadata, expVar, coVars)
  if (!coef_name %in% colnames(fit$p.value)) {
    stop("Could not find the exposure coefficient in limma output.")
  }
  pval_core <- fit$p.value[, coef_name]
  
  return(DA_format_result(feature, expVar, pval_core = pval_core))
}


###################################
# DAssemble Core metagenomeSeq    #
###################################

DA_fit_core_metagenomeSeq <- function(features, metadata, expVar, coVars = NULL, ...) {
  
  ########################
  # Package sanity check #
  ########################
  
  if (!requireNamespace("metagenomeSeq", quietly = TRUE))
    stop("metagenomeSeq required.")
  
  ###################################
  # Standard metagenomeSeq pipeline #
  ###################################
  
  extra_args <- list(...)
  mr_args <- extra_args$newMRexperiment %||% list()
  stat_args <- extra_args$cumNormStat %||% list()
  cumnorm_args <- extra_args$cumNorm %||% list()
  fitzig_args <- extra_args$fitZig %||% list()
  mrcoefs_args <- extra_args$MRcoefs %||% list()

  design <- stats::model.matrix(stats::as.formula(paste("~", build_rhs(expVar, coVars))), metadata)
  count_table <- t(features) 
  mgsdata <- do.call(metagenomeSeq::newMRexperiment, c(list(counts = count_table), mr_args))
  mgsp <- do.call(metagenomeSeq::cumNormStat, c(list(obj = mgsdata), stat_args))
  mgsdata <- do.call(metagenomeSeq::cumNorm, c(list(obj = mgsdata, p = mgsp), cumnorm_args))
  fit <- do.call(metagenomeSeq::fitZig, c(list(obj = mgsdata, mod = design), fitzig_args))
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  coef_name <- get_exp_coef_name(metadata, expVar, coVars)
  coef_table <- metagenomeSeq::MRcoefs(
    fit,
    coef = coef_name,
    number = Inf,
    group = 4,
    adjustMethod = "none"
  )
  if (length(mrcoefs_args) > 0L) {
    coef_table <- do.call(
      metagenomeSeq::MRcoefs,
      c(list(fit = fit, coef = coef_name, number = Inf, group = 4, adjustMethod = "none"), mrcoefs_args)
    )
  }
  feature <- rownames(coef_table)
  pval_core <- coef_table$pvalues
  
  return(DA_format_result(feature, expVar, pval_core = pval_core))
  
}

########################
# DAssemble Core MAST  #
########################

DA_fit_core_MAST <- function(features, metadata, expVar, coVars = NULL, ...) {
  
  ########################
  # Package sanity check #
  ########################
  
  if (!requireNamespace("MAST", quietly = TRUE))
    stop("MAST required.")
  
  ##########################
  # Standard MAST pipeline #
  ##########################
  
  countData <- t(features)
  expr <- log2(edgeR::cpm(countData)+1)
  
  sca <- MAST::FromMatrix(exprsArray = expr)
  cdr2 <- colSums(SummarizedExperiment::assay(sca)>0)
  cd <- SummarizedExperiment::colData(sca)
  cd <- SummarizedExperiment::colData(sca)
  for (nm in colnames(metadata)) {
    cd[[nm]] <- metadata[[nm]]
  }
  cd$cngeneson <- scale(cdr2)
  cd$condition <- droplevels(factor(metadata[[expVar]]))
  SummarizedExperiment::colData(sca) <- cd
  mast_covars <- c(coVars, "cngeneson")
  mast_formula <- stats::as.formula(
    paste("~", build_rhs("condition", mast_covars))
  )
  extra_args <- list(...)
  zlm_args <- extra_args$zlm %||% list()
  lrtest_args <- extra_args$lrTest %||% list()

  zlmCond <- do.call(MAST::zlm, c(list(formula = mast_formula, sca = sca), zlm_args))
  mast_design <- stats::model.matrix(mast_formula, data = as.data.frame(SummarizedExperiment::colData(sca)))
  test.cond <- get_exp_coef_name(
    as.data.frame(SummarizedExperiment::colData(sca)),
    "condition",
    mast_covars
  )
  mast_eval_env <- new.env(parent = asNamespace("MAST"))
  mast_eval_env$test.cond <- test.cond
  mast_eval_env$mast_terms <- colnames(mast_design)
  test.hypothesis <- eval(
    quote(CoefficientHypothesis(test.cond, mast_terms)),
    envir = mast_eval_env
  )
  lrt <- do.call(MAST::lrTest, c(list(object = zlmCond, hypothesis = test.hypothesis), lrtest_args))
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  feature <- dimnames(lrt)$primerid
  pval_core <- lrt[, "hurdle", "Pr(>Chisq)", drop = TRUE]
  return(DA_format_result(feature, expVar, pval_core = pval_core))
}

############################
# DAssemble Core dearseq   #
############################

DA_fit_core_dearseq <- function(features, metadata, expVar, coVars = NULL, ...) {
  
  ########################
  # Package sanity check #
  ########################
  
  if (!requireNamespace("dearseq", quietly = TRUE))
    stop("dearseq required.")
  
  #############################
  # Standard dearseq pipeline #
  #############################
  
  se_raw <- SummarizedExperiment::SummarizedExperiment(assays = as.matrix(t(features)), colData = metadata)
  which_test <- if (nrow(metadata) <= 20) "permutation" else "asymptotic"
  
  # Build covariate matrix if coVars provided
  cov_mat <- build_covariate_matrix(metadata, coVars)
  
  fit <- do.call(
    dearseq::dear_seq,
    c(list(
      object = se_raw,
      variables2test = expVar,
      covariates = cov_mat,
      which_test = which_test,
      preprocessed = FALSE,
      progressbar = FALSE,
      parallel_comp = FALSE,
      which_weights = "loclin"
    ), list(...))
  )
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  feature <- rownames(fit$pvals)
  pval_core <- fit$pvals$rawPval
  return(DA_format_result(feature, expVar, pval_core = pval_core))
}

#########################
# DAssemble Core Robseq #
#########################

DA_fit_core_Robseq <- function(features, metadata, expVar, coVars = NULL, ...) {
  
  ########################
  # Package sanity check #
  ########################
  
  required_pkgs <- c("DESeq2", "edgeR", "MASS", "dfadjust", "preprocessCore")
  missing_pkgs <- required_pkgs[
    !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_pkgs) > 0L) {
    stop("Robseq core requires: ", toString(missing_pkgs))
  }
  
  ############################
  # Standard Robseq pipeline #
  ############################
  
  countMat <- t(as.matrix(features))
  meta_rob <- metadata
  meta_rob$Exposure <- metadata[[expVar]]
  fit <- DAssemble_robust_dge(
    features = countMat,
    metadata = meta_rob,
    expVar = "Exposure",
    coVars = coVars,
    ...
  )
  res <- as.data.frame(fit$res)
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  feature <- rownames(res)
  pval_core <- res$Pval
  return(DA_format_result(feature, expVar, pval_core = pval_core))
}

##############################
# DAssemble Core Prevalence  #
##############################

DA_prevalence_augmentation_weight <- function(formula, df) {
  rhs_formula <- stats::delete.response(stats::terms(formula))
  n_predictors <- ncol(stats::model.matrix(rhs_formula, data = df))
  n_samples <- nrow(df)
  n_predictors / (2 * n_samples)
}

DA_prevalence_augment_data <- function(formula,
                                       df,
                                       response = "expr",
                                       weights = NULL) {
  n_samples <- nrow(df)
  if (is.null(weights)) {
    weights <- rep(1, n_samples)
  }
  if (length(weights) != n_samples) {
    stop("weights must have length equal to the number of samples.")
  }

  augmentation_weight <- DA_prevalence_augmentation_weight(formula, df)
  df_zero <- df
  df_one <- df
  df_zero[[response]] <- 0L
  df_one[[response]] <- 1L

  list(
    data = rbind(df, df_zero, df_one),
    weights = c(
      weights,
      rep(augmentation_weight, n_samples),
      rep(augmentation_weight, n_samples)
    )
  )
}

DA_prevalence_fit_augmented <- function(formula,
                                        df,
                                        has_random_effects = FALSE,
                                        ...) {
  extra_args <- list(...)
  user_weights <- extra_args$weights
  extra_args$weights <- NULL
  augmented <- DA_prevalence_augment_data(
    formula = formula,
    df = df,
    weights = user_weights
  )

  fit_fun <- if (isTRUE(has_random_effects)) glmmTMB::glmmTMB else stats::glm
  do.call(
    fit_fun,
    c(
      list(
        formula = formula,
        family = stats::binomial(),
        data = augmented$data,
        weights = augmented$weights
      ),
      extra_args
    )
  )
}

DA_prevalence_fit_firth <- function(formula, df, ...) {
  if (!requireNamespace("brglm2", quietly = TRUE)) {
    stop(
      "brglm2 is required for separation_method = 'firth' in the Prevalence core."
    )
  }

  do.call(
    stats::glm,
    c(
      list(
        formula = formula,
        family = stats::binomial(),
        data = df,
        method = brglm2::brglmFit,
        type = "AS_mean"
      ),
      list(...)
    )
  )
}

DA_fit_core_Prevalence <- function(features,
                                   metadata,
                                   expVar,
                                   coVars = NULL,
                                   random_effects = NULL,
                                   zero_threshold = 0,
                                   ...) {
  extra_args <- list(...)
  separation_method <- extra_args$separation_method %||% "augment"
  extra_args$separation_method <- NULL
  if (!separation_method %in% c("augment", "firth")) {
    stop("Prevalence separation_method must be one of 'augment' or 'firth'.")
  }
  formula <- build_model_formula(
    expVar = expVar,
    coVars = coVars,
    random_effects = random_effects,
    response = "expr"
  )
  coef_name <- get_exp_coef_name(metadata, expVar, coVars)
  has_random_effects <- !is.null(random_effects) && length(random_effects) > 0L

  if (has_random_effects && !requireNamespace("glmmTMB", quietly = TRUE)) {
    stop("glmmTMB is required for prevalence models with random_effects.")
  }
  if (has_random_effects && identical(separation_method, "firth")) {
    stop("Prevalence separation_method = 'firth' is only supported without random_effects.")
  }

  prevalence_stats <- vapply(seq_len(ncol(features)), function(j) {
    expr <- as.integer(features[, j] > zero_threshold)
    if (length(unique(expr)) < 2L) {
      return(c(coef = NA_real_, pval = NA_real_))
    }

    df <- cbind(metadata, expr = expr)
    fit <- if (has_random_effects || identical(separation_method, "augment")) {
      try(
        do.call(
          DA_prevalence_fit_augmented,
          c(
            list(
              formula = formula,
              df = df,
              has_random_effects = has_random_effects
            ),
            extra_args
          )
        ),
        silent = TRUE
      )
    } else {
      try(
        do.call(
          DA_prevalence_fit_firth,
          c(list(formula = formula, df = df), extra_args)
        ),
        silent = TRUE
      )
    }

    if (inherits(fit, "try-error")) {
      return(c(coef = NA_real_, pval = NA_real_))
    }

    sm <- try(summary(fit), silent = TRUE)
    if (inherits(sm, "try-error")) {
      return(c(coef = NA_real_, pval = NA_real_))
    }

    coef_table <- if (has_random_effects) sm$coefficients$cond else stats::coef(summary(fit))
    if (!coef_name %in% rownames(coef_table)) {
      return(c(coef = NA_real_, pval = NA_real_))
    }

    p_col <- intersect(c("Pr(>|z|)", "Pr(>|t|)"), colnames(coef_table))[1]
    if (is.na(p_col)) {
      return(c(coef = NA_real_, pval = NA_real_))
    }

    c(
      coef = coef_table[coef_name, "Estimate"],
      pval = coef_table[coef_name, p_col]
    )
  }, numeric(2))

  coef_prev <- prevalence_stats["coef", ]
  pval_prev <- prevalence_stats["pval", ]
  names(coef_prev) <- colnames(features)
  names(pval_prev) <- colnames(features)

  DA_format_result(
    feature = colnames(features),
    expVar = expVar,
    coef_core = coef_prev,
    pval_core = pval_prev
  )
}

###########################
# DAssemble Core MaAsLin2 #
###########################

DA_fit_core_Maaslin2 <- function(features,
                                 metadata,
                                 expVar,
                                 coVars = NULL,
                                 random_effects = NULL,
                                 ...) {
  
  ########################
  # Package sanity check #
  ########################
  
  if (!requireNamespace("Maaslin2", quietly = TRUE))
    stop("Maaslin2 required.")
  
  ##############################
  # Standard MaAsLin2 pipeline #
  ##############################
  
  tmp <- file.path(tempdir(), paste0("m2_", sample(1e8,1)))
  dir.create(tmp, showWarnings = FALSE)
  fit <- Maaslin2::Maaslin2(features, 
                            metadata, 
                            output = tmp, 
                            fixed_effects  = c(expVar, coVars),
                            random_effects = random_effects,
                            min_abundance = -Inf, # No additional filtering
                            save_scatter = FALSE, 
                            save_models = FALSE,
                            plot_heatmap = FALSE, 
                            plot_scatter = FALSE, 
                            max_significance = 1,
                            ...)
  res <- fit$results
  if (requireNamespace("logging", quietly = TRUE)) {
    try(logging::removeHandler("logging::writeToFile"), silent = TRUE)
  }
  unlink(tmp, recursive = TRUE)
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  res <- res[res$metadata == expVar, , drop = FALSE]
  feature   <- res$feature
  coef_core <- res$coef
  pval_core <- res$pval
  return(DA_format_result(
    feature,
    expVar,
    coef_core = coef_core,
    pval_core = pval_core
  ))
}

###########################
# DAssemble Core MaAsLin3 #
###########################

DA_fit_core_Maaslin3 <- function(features,
                                 metadata,
                                 expVar,
                                 coVars = NULL,
                                 random_effects = NULL,
                                 prevalence_only = FALSE,
                                 ...) {
  
  ########################
  # Package sanity check #
  ########################
  
  if (!requireNamespace("maaslin3", quietly = TRUE))
    stop("maaslin3 required.")
  
  ##############################
  # Standard MaAsLin3 pipeline #
  ##############################
  
  tmp <- file.path(tempdir(), paste0("m3_", sample(1e8,1)))
  dir.create(tmp, showWarnings = FALSE)
  fit <- maaslin3::maaslin3(
    features,
    metadata,
    output = tmp,
    fixed_effects  = c(expVar, coVars),
    random_effects = random_effects,
    min_prevalence = -Inf, # No additional filtering
    evaluate_only = if (isTRUE(prevalence_only)) "prevalence" else NULL,
    median_comparison_abundance = TRUE, 
    median_comparison_prevalence = TRUE,
    subtract_median = TRUE,
    plot_summary_plot = FALSE, 
    plot_associations = FALSE, 
    max_significance = 1,
    ...)
  if (requireNamespace("logging", quietly = TRUE)) {
    try(logging::removeHandler("logging::writeToFile"), silent = TRUE)
  }
  unlink(tmp, recursive = TRUE)
  res <- if (isTRUE(prevalence_only)) {
    as.data.frame(fit$fit_data_prevalence$results)
  } else {
    as.data.frame(fit$fit_data_abundance$results)
  }
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  res <- res[res$metadata == expVar, , drop = FALSE]
  feature<-res$feature
  pval_core <- if (isTRUE(prevalence_only)) res$pval else res$pval_joint
  return(DA_format_result(feature, expVar, pval_core = pval_core))
}



###########################
# DAssemble Core ANCOMBC2 #
###########################

DA_fit_core_ANCOMBC2 <- function(features, metadata, expVar, coVars = NULL, ...) {
  
  ########################
  # Package sanity check #
  ########################
  
  if (!requireNamespace("ANCOMBC", quietly = TRUE))
    stop("ANCOMBC required.")
  
  ##############################
  # Standard ANCOMBC2 pipeline #
  ##############################
  
  otu_data<-as.data.frame(t(features))
  metadata$sampleID<-rownames(metadata) # ANCOMBC2 does not work with single-column metadata
  
  fix_formula <- build_rhs(expVar, coVars)
  
  fit <- ANCOMBC::ancombc2(data = otu_data, 
                           meta_data = metadata, 
                           fix_formula = fix_formula,
                           prv_cut = 0, # No additional filtering
                           alpha = 1,
                           ...)
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  res <- as.data.frame(fit$res)
  p_col <- grep(paste0("^p_", expVar), names(res), value = TRUE)
  feature<-res$taxon
  pval_core <- res[[p_col]]
  return(DA_format_result(feature, expVar, pval_core = pval_core))
}

#########################
# DAssemble Core ALDEx2 #
#########################

DA_fit_core_ALDEx2 <- function(features, metadata, expVar, coVars = NULL, ...) {
  
  ########################
  # Package sanity check #
  ########################
  
  if (!requireNamespace("ALDEx2", quietly = TRUE))
    stop("ALDEx2 required.")
  
  ############################
  # Standard ALDEx2 pipeline #
  ############################
  
  group <- droplevels(factor(metadata[[expVar]]))
  
  if (is.null(coVars)) {
    
    # No covariates: use the standard aldex() convenience wrapper
    ald       <- do.call(
      ALDEx2::aldex,
      c(list(reads = t(as.matrix(features)), conditions = as.character(group)), list(...))
    )
    feature   <- rownames(ald)
    pval_core <- ald$we.ep
    
  } else {
    
    # Covariates present: use aldex.clr + aldex.glm
    clr_args <- list(...)
    glm_args <- clr_args$aldex.glm %||% list()
    clr_call_args <- clr_args$aldex.clr %||% list()
    clr_obj <- do.call(
      ALDEx2::aldex.clr,
      c(list(reads = t(as.matrix(features)), conds = as.character(group)), clr_call_args)
    )
    mm      <- model.matrix(as.formula(paste("~", build_rhs(expVar, coVars))), metadata)
    glm_res <- do.call(ALDEx2::aldex.glm, c(list(clr = clr_obj, mm = mm), glm_args))
    
    # Extract p-value column corresponding to expVar
    p_col   <- grep(paste0("^model\\.", expVar, ".*Pr"), colnames(glm_res), value = TRUE)[1]
    if (is.na(p_col))
      stop("Could not find a p-value column for expVar '", expVar, "' in aldex.glm output.")
    
    feature   <- rownames(glm_res)
    pval_core <- glm_res[[p_col]]
  }
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  return(DA_format_result(feature, expVar, pval_core = pval_core))
}



########################
# DAssemble Core LinDA #
########################

DA_fit_core_LinDA <- function(features, metadata, expVar, coVars = NULL, ...) {
  
  ########################
  # Package sanity check #
  ########################
  
  if (!requireNamespace("MicrobiomeStat", quietly = TRUE))
    stop("MicrobiomeStat required for LinDA.")
  
  ###########################
  # Standard LinDA pipeline #
  ###########################
  
  otu <- t(as.matrix(features))
  out <- MicrobiomeStat::linda(
    feature.dat = otu,
    meta.dat = metadata,
    formula = paste("~", build_rhs(expVar, coVars)),
    feature.dat.type = "count",
    prev.filter = 0,
    mean.abund.filter = 0,
    max.abund.filter = 0,
    is.winsor = FALSE,
    verbose = FALSE,
    ...
  )
  res <- as.data.frame(out$output[[1]])
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  feature <- rownames(res)
  pval_core <- res$pvalue
  return(DA_format_result(feature, expVar, pval_core = pval_core))
  
}


########################
# DAssemble Core LOCOM #
########################

DA_fit_core_LOCOM <- function(features, metadata, expVar, coVars = NULL, ...) {
  
  ########################
  # Package sanity check #
  ########################
  
  if (!requireNamespace("LOCOM2", quietly = TRUE))
    stop("LOCOM2 required.")
  
  ###########################
  # Standard LOCOM pipeline #
  ###########################
  
  otu <- as.matrix(features)
  Y <- metadata[[expVar]]
  if (is.factor(Y)) Y <- as.numeric(Y) - 1
  fit <- LOCOM2::locom2(
    otu.table = otu,
    Y = Y,
    seed = 1234,
    filter = FALSE,
    verbose = FALSE,
    ...
  )
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  feature <- names(fit$p.otu.Wald)
  if (is.null(feature)) {
    feature <- colnames(otu)
  }
  pval_core <- as.vector(fit$p.otu.Wald)
  return(DA_format_result(feature, expVar, pval_core = pval_core))
  
}

##############################
# DAssemble Core Tweedieverse #
##############################

DA_fit_core_Tweedieverse <- function(features,
                                     metadata,
                                     expVar,
                                     coVars = NULL,
                                     random_effects = NULL,
                                     ...) {
  
  ########################
  # Package sanity check #
  ########################
  
  required_pkgs <- c("cplm", "glmmTMB", "logging", "pbapply")
  missing_pkgs <- required_pkgs[
    !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_pkgs) > 0L) {
    stop("Tweedieverse core requires: ", toString(missing_pkgs))
  }
  
  ##################################
  # Standard Tweedieverse pipeline #
  ##################################
  
  res <- DAssemble_Tweedieverse(
    features,
    metadata,
    output = NULL,
    max_significance = 1,
    abd_threshold = -Inf, # No additional filtering
    fixed_effects = paste(c(expVar, coVars), collapse = ","),
    random_effects = if (length(random_effects) > 0L) paste(random_effects, collapse = ",") else NULL,
    adjust_offset = TRUE,
    median_comparison = FALSE,
    median_subtraction = FALSE,
    na.action = na.pass,
    ...
  )
  
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  res <- res[res$metadata == expVar, , drop = FALSE]
  feature<-res$feature
  pval_core <- res$pval
  return(DA_format_result(feature, expVar, pval_core = pval_core))
  
}
