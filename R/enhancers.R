###############################################
# DAssemble Enhancer Logistic Regression (LR) #
###############################################

DA_fit_enhancer_LR <- function(features, metadata, expVar, coVars = NULL) {
  
  ########################
  # Standard LR pipeline #
  ########################
  
  # Offset: scale_factor if present, otherwise library size
  if ("scale_factor" %in% colnames(metadata)) {
    offset_raw <- metadata[, "scale_factor"]
  } else {
    offset_raw <- rowSums(features)
  }
  log_offset <- log(offset_raw)
  
  # Build formula with expVar only
  formula <- as.formula(paste("expr ~", build_rhs(expVar, coVars)))
  coef_name <- get_exp_coef_name(metadata, expVar, coVars)
  
  ##################
  # Per-feature LR #
  ##################
  
  lr_stats <- vapply(seq_len(ncol(features)), function(j) {
    expr <- as.integer(features[, j] > 0)
    df   <- cbind(metadata, expr = expr)
    fit  <- try(glm(formula = formula, family = binomial(),
                    data = df, offset = log_offset), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      sm <- coef(summary(fit))
      if (coef_name %in% rownames(sm)) {
        return(c(coef = sm[coef_name, 1], pval = sm[coef_name, 4]))
      }
    }
    c(coef = NA_real_, pval = NA_real_)
  }, numeric(2))
  coef_LR <- lr_stats["coef", ]
  pval_LR <- lr_stats["pval", ]
  names(coef_LR) <- colnames(features)
  names(pval_LR) <- colnames(features)
  
  ###########################################
  # Standardized output - Feature + Pvalues #
  ###########################################
  
  feature <- names(pval_LR)
  return(DA_format_result(
    feature,
    expVar,
    coef_LR = coef_LR,
    pval_LR = pval_LR
  ))
}


##############################################
# DAssemble Enhancer Kolmogorov–Smirnov (KS) #
##############################################

DA_fit_enhancer_KS <- function(features, metadata, expVar, coVars = NULL) {

  ########################
  # Standard KS pipeline #
  ########################
  if (!is.null(coVars)) warning("KS is a nonparametric two-group test; coVars are ignored.")
  
  group <- metadata[[expVar]]
  g1 <- which(group == levels(group)[1])
  g2 <- which(group == levels(group)[2])
  
  ##################
  # Per-feature KS #
  ##################
  
  pval_KS <- vapply(seq_len(ncol(features)), function(j) {
    x1 <- features[g1, j]
    x2 <- features[g2, j]
    tst <- try(stats::ks.test(x1, x2), silent = TRUE)
    if (inherits(tst, "try-error")) {
      NA_real_
    } else {
      tst$p.value
    }
  }, numeric(1))
  names(pval_KS)<-colnames(features)
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  feature <- names(pval_KS)
  return(DA_format_result(feature, expVar, pval_KS = pval_KS))
}


##############################################
# DAssemble Enhancer Wilcoxon Rank Sum (WLX) #
##############################################

DA_fit_enhancer_WLX <- function(features, metadata, expVar, coVars = NULL) {
 
  #########################
  # Standard WLX pipeline #
  #########################
  if (!is.null(coVars)) warning("WLX is a nonparametric two-group test; coVars are ignored.")

  group <- metadata[[expVar]]
  g1 <- which(group == levels(group)[1])
  g2 <- which(group == levels(group)[2])
  
  ###################
  # Per-feature WLX #
  ###################
  
  pval_WLX <- vapply(seq_len(ncol(features)), function(j) {
    x1 <- features[g1, j]
    x2 <- features[g2, j]
    tst <- try(stats::wilcox.test(x1, x2), silent = TRUE)
    if (inherits(tst, "try-error")) {
      NA_real_
    } else {
      tst$p.value
    }
  }, numeric(1))
  names(pval_WLX)<-colnames(features)
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  feature <- names(pval_WLX)
  return(DA_format_result(feature, expVar, pval_WLX = pval_WLX))
}
