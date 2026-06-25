#########################
# DAssemble Core DESeq2 #
#########################

DA_fit_core_DESeq2 <- function(features, metadata, expVar, coVars = NULL) {
  
  ########################
  # Package sanity check #
  ########################
  
  if (!requireNamespace("DESeq2", quietly = TRUE))
    stop("DESeq2 required.")
  
  ############################
  # Standard DESeq2 pipeline #
  ############################
  
  design <- stats::as.formula(paste("~", expVar))
  x <- DESeq2::DESeqDataSetFromMatrix(countData = t(as.matrix(features)), colData = metadata, design = design)
  gm_mean <- function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans <- apply(DESeq2::counts(x), 1, gm_mean)
  x <- DESeq2::estimateSizeFactors(x, geoMeans = geoMeans)
  fit <- DESeq2::DESeq(x)
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  feature<-rownames(stats::coef(fit))
  pval_core<-DESeq2::results(fit,name=DESeq2::resultsNames(fit)[2])$pvalue
  return(DA_format_result(feature, expVar, pval_core = pval_core))
}

########################
# DAssemble Core edgeR #
########################

DA_fit_core_edgeR <- function(features, metadata, expVar, coVars = NULL) {
  
  ########################
  # Package sanity check #
  ########################
  
  if (!requireNamespace("edgeR", quietly = TRUE))
    stop("edgeR required.")
  
  ###########################
  # Standard edgeR pipeline #
  ###########################
  
  d <- edgeR::DGEList(counts = t(features))
  d <- edgeR::calcNormFactors(d, method='TMM')
  design <- stats::model.matrix(stats::as.formula(paste("~", expVar)), metadata)
  d <- edgeR::estimateGLMCommonDisp(d, design)
  d <- edgeR::estimateGLMTrendedDisp(d, design)
  d <- edgeR::estimateGLMTagwiseDisp(d, design)
  fit <- edgeR::glmFit(d, design)
  fit<- edgeR::glmLRT(fit, 2)
  
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

DA_fit_core_limmaVOOM <- function(features, metadata, expVar, coVars = NULL) {
  
  ########################
  # Package sanity check #
  ########################
  
  if (!requireNamespace("limma", quietly = TRUE) ||
      !requireNamespace("edgeR", quietly = TRUE))
    stop("limma + edgeR required.")
  
  ###############################
  # Standard limmaVOOM pipeline #
  ###############################
  
  design <- stats::model.matrix(stats::as.formula(paste("~", expVar)), metadata)
  x<-t(as.matrix(features)+1) # Convert to matrix, round up to nearest integer, and transpose
  y <- limma::voom(x,design,plot=FALSE)
  fit <- limma::lmFit(y,design)
  fit <- limma::eBayes(fit)
  
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

DA_fit_core_metagenomeSeq <- function(features, metadata, expVar, coVars = NULL) {
  
  ########################
  # Package sanity check #
  ########################
  
  if (!requireNamespace("metagenomeSeq", quietly = TRUE))
    stop("metagenomeSeq required.")
  
  ###################################
  # Standard metagenomeSeq pipeline #
  ###################################
  
  design <- stats::model.matrix(stats::as.formula(paste("~", expVar)), metadata)
  count_table <- t(features) 
  mgsdata <- metagenomeSeq::newMRexperiment(counts = count_table)
  mgsp <- metagenomeSeq::cumNormStat(mgsdata)
  mgsdata <- metagenomeSeq::cumNorm(mgsdata, mgsp)
  fit <- metagenomeSeq::fitZig(obj=mgsdata,mod=design)
  
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
  feature <- rownames(coef_table)
  pval_core <- coef_table$pvalues
  
  return(DA_format_result(feature, expVar, pval_core = pval_core))
  
}

########################
# DAssemble Core MAST  #
########################

DA_fit_core_MAST <- function(features, metadata, expVar, coVars = NULL) {
  
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
  cd$cngeneson <- scale(cdr2)
  cd$condition <- droplevels(factor(metadata[[expVar]]))
  SummarizedExperiment::colData(sca) <- cd
  zlmCond <- MAST::zlm(~condition + cngeneson, sca)
  mast_design <- stats::model.matrix(
    ~condition + cngeneson,
    data = as.data.frame(SummarizedExperiment::colData(sca))
  )
  test.cond <- colnames(mast_design)[2L]
  mast_eval_env <- new.env(parent = asNamespace("MAST"))
  mast_eval_env$test.cond <- test.cond
  mast_eval_env$mast_terms <- colnames(mast_design)
  test.hypothesis <- eval(
    quote(CoefficientHypothesis(test.cond, mast_terms)),
    envir = mast_eval_env
  )
  lrt <- MAST::lrTest(zlmCond, test.hypothesis)
  
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

DA_fit_core_dearseq <- function(features, metadata, expVar, coVars = NULL) {
  
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
  
  fit <- dearseq::dear_seq(object = se_raw, variables2test = expVar, covariates = cov_mat,
                           which_test = which_test, preprocessed = FALSE,
                           progressbar = FALSE, parallel_comp = FALSE, 
                           which_weights = "loclin")
  
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

DA_fit_core_Robseq <- function(features, metadata, expVar, coVars = NULL) {
  
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
  fit <- DAssemble_robust_dge(features = countMat, metadata = meta_rob)
  res <- as.data.frame(fit$res)
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  feature <- rownames(res)
  pval_core <- res$Pval
  return(DA_format_result(feature, expVar, pval_core = pval_core))
}

###########################
# DAssemble Core MaAsLin2 #
###########################

DA_fit_core_Maaslin2 <- function(features, metadata, expVar, coVars = NULL) {
  
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
                            min_abundance = -Inf, # No additional filtering
                            save_scatter = FALSE, 
                            save_models = FALSE,
                            plot_heatmap = FALSE, 
                            plot_scatter = FALSE, 
                            max_significance = 1)
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

DA_fit_core_Maaslin3 <- function(features, metadata, expVar, coVars = NULL) {
  
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
    min_prevalence = -Inf, # No additional filtering
    median_comparison_abundance = TRUE, 
    median_comparison_prevalence = TRUE,
    subtract_median = TRUE,
    plot_summary_plot = FALSE, 
    plot_associations = FALSE, 
    max_significance = 1)
  if (requireNamespace("logging", quietly = TRUE)) {
    try(logging::removeHandler("logging::writeToFile"), silent = TRUE)
  }
  unlink(tmp, recursive = TRUE)
  res <- as.data.frame(fit$fit_data_abundance$results)
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  res <- res[res$metadata == expVar, , drop = FALSE]
  feature<-res$feature
  pval_core <- res$pval_joint
  return(DA_format_result(feature, expVar, pval_core = pval_core))
}



###########################
# DAssemble Core ANCOMBC2 #
###########################

DA_fit_core_ANCOMBC2 <- function(features, metadata, expVar, coVars = NULL) {
  
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
                           alpha = 1)
  
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

DA_fit_core_ALDEx2 <- function(features, metadata, expVar, coVars = NULL) {
  
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
    ald       <- ALDEx2::aldex(t(as.matrix(features)), conditions = as.character(group))
    feature   <- rownames(ald)
    pval_core <- ald$we.ep
    
  } else {
    
    # Covariates present: use aldex.clr + aldex.glm
    clr_obj <- ALDEx2::aldex.clr(t(as.matrix(features)),
                                 conds = as.character(group))
    mm      <- model.matrix(as.formula(paste("~", build_rhs(expVar, coVars))), metadata)
    glm_res <- ALDEx2::aldex.glm(clr_obj, mm)
    
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

DA_fit_core_LinDA <- function(features, metadata, expVar, coVars = NULL) {
  
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
    formula = paste("~", expVar),
    feature.dat.type = "count",
    prev.filter = 0,
    mean.abund.filter = 0,
    max.abund.filter = 0,
    is.winsor = FALSE,
    verbose = FALSE
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

DA_fit_core_LOCOM <- function(features, metadata, expVar, coVars = NULL) {
  
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
    verbose = FALSE
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

DA_fit_core_Tweedieverse <- function(features, metadata, expVar, coVars = NULL) {
  
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
    fixed_effects = expVar,
    adjust_offset = TRUE,
    median_comparison = FALSE,
    median_subtraction = FALSE,
    na.action = na.pass
  )
  
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  res <- res[res$metadata == expVar, , drop = FALSE]
  feature<-res$feature
  pval_core <- res$pval
  return(DA_format_result(feature, expVar, pval_core = pval_core))
  
}
