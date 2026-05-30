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
  
  design <- as.formula(paste("~", build_rhs(expVar, coVars)))
  x <- DESeq2::DESeqDataSetFromMatrix(countData = t(as.matrix(features)), colData = metadata, design = design)
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(x), 1, gm_mean)
  x = estimateSizeFactors(x, geoMeans = geoMeans)
  fit <- DESeq(x)
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  feature<-rownames(coef(fit))
  pval_core<-results(fit,name=resultsNames(fit)[2])$pvalue
  df<-cbind.data.frame(feature, pval_core)
  df$metadata<- expVar
  df<-dplyr::select(df, c('feature', 'metadata'), everything())
  return(df)
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
  
  d <- DGEList(counts = t(features))
  d <- suppressWarnings(edgeR::calcNormFactors(d, method='TMM'))
  design <- build_design_matrix(metadata, expVar, coVars)
  # d <- estimateDisp(d, design) -> This step is equivalent to the next three steps
  d <- estimateGLMCommonDisp(d, design)
  d <- estimateGLMTrendedDisp(d, design)
  d <- estimateGLMTagwiseDisp(d, design)
  
  fit <- glmFit(d, design)
  
  coef_name <- get_exp_coef_name(metadata, expVar, coVars)
  coef_idx  <- match(coef_name, colnames(design))
  if (is.na(coef_idx)) {
    stop("Could not match the exposure coefficient in edgeR design matrix.")
  }
  
  fit <- glmLRT(fit, coef = coef_idx)
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  feature<-rownames(fit$table)
  pval_core<-fit$table$PValue
  df<-cbind.data.frame(feature, pval_core)
  df$metadata<- expVar
  df<-dplyr::select(df, c('feature', 'metadata'), everything())
  return(df)
  
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
  
  design <- build_design_matrix(metadata, expVar, coVars)
  x<-t(as.matrix(features)+1) # Convert to matrix, round up to nearest integer, and transpose
  y <- voom(x,design,plot=FALSE)
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
  
  df<-cbind.data.frame(feature, pval_core)
  df$metadata<- expVar
  df<-dplyr::select(df, c('feature', 'metadata'), everything())
  return(df)
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
  
  design <- build_design_matrix(metadata, expVar, coVars)
  count_table <- t(features) 
  mgsdata <- newMRexperiment(counts = count_table)
  mgsp <- cumNormStat(mgsdata)
  mgsdata <- cumNorm(mgsdata, mgsp)
  fit <- fitZig(obj=mgsdata,mod=design)
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  feature<-rownames(fit@fit$coefficients)
  coef_name <- get_exp_coef_name(metadata, expVar, coVars)
  if (!coef_name %in% colnames(fit@eb$p.value)) {
    stop("Could not find the exposure coefficient in metagenomeSeq output.")
  }
  pval_core <- fit@eb$p.value[, coef_name]
  
  df<-cbind.data.frame(feature, pval_core)
  df$metadata<- expVar
  df<-dplyr::select(df, c('feature', 'metadata'), everything())
  return(df)
  
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
  cdr2 <- colSums(assay(sca) > 0)
  
  colData(sca)$cngeneson <- as.numeric(scale(cdr2))
  colData(sca)$condition <- droplevels(factor(metadata[[expVar]]))
  
  if (!is.null(coVars) && length(coVars) > 0L) {
    for (cv in coVars) {
      colData(sca)[[cv]] <- metadata[[cv]]
    }
  }
  
  cov_terms <- if (!is.null(coVars) && length(coVars) > 0L) {
    paste("+", paste(coVars, collapse = " + "))
  } else {
    ""
  }
  
  zlm_formula <- as.formula(paste("~ condition + cngeneson", cov_terms))
  zlmCond <- MAST::zlm(formula = zlm_formula, sca = sca)
  
  test.cond <- grep("^condition", colnames(zlmCond@coefC), value = TRUE)[1]
  if (is.na(test.cond)) {
    stop("Could not identify the exposure coefficient in MAST.")
  }
  
  summaryCond <- summary(zlmCond, doLRT = test.cond)
  summaryDt <- summaryCond$datatable
  
  fcHurdle <- merge(
    summaryDt[contrast == test.cond & component == "H",
              .(primerid, `Pr(>Chisq)`)],
    summaryDt[contrast == test.cond & component == "logFC",
              .(primerid, coef, ci.hi, ci.lo)],
    by = "primerid"
  )
  
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  feature<-fcHurdle$primerid
  pval_core<-fcHurdle[,2]
  df<-cbind.data.frame(feature, pval_core)
  names(df)[2] <- c('pval_core')
  df$metadata<- expVar
  df<-dplyr::select(df, c('feature', 'metadata'), everything())
  return(df)
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
  
  se_raw <- SummarizedExperiment(assays = as.matrix(t(features)), colData = metadata)
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
  df<-cbind.data.frame(feature, pval_core)
  df$metadata<- expVar
  df<-dplyr::select(df, c('feature', 'metadata'), everything())
  return(df)
}

#########################
# DAssemble Core Robseq #
#########################

DA_fit_core_Robseq <- function(features, metadata, expVar, coVars = NULL) {
  
  ########################
  # Package sanity check #
  ########################
  
  if (!requireNamespace("Robseq", quietly = TRUE))
    stop("Robseq required.")
  
  ############################
  # Standard Robseq pipeline #
  ############################
  
  countMat <- t(as.matrix(features))
  meta_rob <- metadata
  meta_rob$Exposure <- metadata[[expVar]]
  fit <- Robseq::robust.dge(features = countMat, metadata = meta_rob, coVars = coVars)
  res <- as.data.frame(fit$res)
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  feature <- rownames(res)
  pval_core <- res$Pval
  df<-cbind.data.frame(feature, pval_core)
  df$metadata<- expVar
  df<-dplyr::select(df, c('feature', 'metadata'), everything())
  return(df)
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
  unlink(tmp, recursive = TRUE)
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  res <- res[res$metadata == expVar, , drop = FALSE]
  feature   <- res$feature
  coef_core <- res$coef
  pval_core <- res$pval
  df <- cbind.data.frame(feature, coef_core, pval_core)
  
  df$metadata<- expVar
  df<-dplyr::select(df, c('feature', 'metadata'), everything())
  return(df)
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
  unlink(tmp, recursive = TRUE)
  res <- as.data.frame(fit$fit_data_abundance$results)
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  res <- res[res$metadata == expVar, , drop = FALSE]
  feature<-res$feature
  pval_core <- res$pval_joint
  df<-cbind.data.frame(feature, pval_core)
  df$metadata<- expVar
  df<-dplyr::select(df, c('feature', 'metadata'), everything())
  return(df)
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
  df<-cbind.data.frame(feature, pval_core)
  df$metadata<- expVar
  df<-dplyr::select(df, c('feature', 'metadata'), everything())
  return(df)
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
  
  df <- cbind.data.frame(feature, pval_core)
  df$metadata <- expVar
  df <- dplyr::select(df, c('feature', 'metadata'), everything())
  return(df)
}



########################
# DAssemble Core LinDA #
########################

DA_fit_core_LinDA <- function(features, metadata, expVar, coVars = NULL) {
  
  ########################
  # Package sanity check #
  ########################
  
  if (!requireNamespace("LinDA", quietly = TRUE))
    stop("LinDA required.")
  
  ###########################
  # Standard LinDA pipeline #
  ###########################
  
  otu <- t(as.matrix(features))
  out <- LinDA::linda(otu.tab = otu, meta = metadata, formula = paste("~", build_rhs(expVar, coVars)))
  res <- as.data.frame(out$output[[1]])
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  feature <- rownames(res)
  pval_core <- res$pvalue
  df<-cbind.data.frame(feature, pval_core)
  df$metadata<- expVar
  df<-dplyr::select(df, c('feature', 'metadata'), everything())
  return(df)
  
}


########################
# DAssemble Core LOCOM #
########################

DA_fit_core_LOCOM <- function(features, metadata, expVar, coVars = NULL) {
  
  ########################
  # Package sanity check #
  ########################
  
  if (!requireNamespace("LOCOM", quietly = TRUE))
    stop("LOCOM required.")
  
  ###########################
  # Standard LOCOM pipeline #
  ###########################
  
  otu <- as.matrix(features)
  Y <- metadata[[expVar]]
  if (is.factor(Y)) Y <- as.numeric(Y) - 1
  
  C_mat <- build_covariate_matrix(metadata, coVars)
  
  fit <- LOCOM::locom(otu.table = otu, Y = Y, C = C_mat, seed = 1234, filter.thresh = 0) # No additional filtering
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  feature <- colnames(fit$p.otu)
  pval_core <- as.vector(fit$p.otu)
  df<-cbind.data.frame(feature, pval_core)
  df$metadata<- expVar
  df<-dplyr::select(df, c('feature', 'metadata'), everything())
  return(df)
  
}

##############################
# DAssemble Core Tweedieverse #
##############################

DA_fit_core_Tweedieverse <- function(features, metadata, expVar, coVars = NULL) {
  
  ########################
  # Package sanity check #
  ########################
  
  if (!requireNamespace("Tweedieverse", quietly = TRUE))
    stop("Tweedieverse required.")
  
  ##################################
  # Standard Tweedieverse pipeline #
  ##################################
  
  tmp <- file.path(tempdir(), paste0("tw_", sample(1e8,1)))
  dir.create(tmp, showWarnings = FALSE)
  res <- Tweedieverse(features,
                      metadata,
                      output = tmp,
                      max_significance = 1,
                      abd_threshold = -Inf, # No additional filtering
                      fixed_effects  = c(expVar, coVars),
                      adjust_offset = TRUE, # Expects a variable named scale_factor
                      median_comparison = FALSE,
                      median_subtraction = FALSE,
                      na.action = na.pass)
  
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  res <- res[res$metadata == expVar, , drop = FALSE]
  feature<-res$feature
  pval_core <- res$pval
  df<-cbind.data.frame(feature, pval_core)
  df$metadata<- expVar
  df<-dplyr::select(df, c('feature', 'metadata'), everything())
  return(df)
  
}
