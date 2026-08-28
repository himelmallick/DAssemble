###############################################
# DAssemble Enhancer Logistic Regression (LR) #
###############################################

DA_fit_enhancer_LR <- function(features,
                               metadata,
                               expVar,
                               coVars = NULL,
                               random_effects = NULL,
                               ...) {
  extra_args <- list(...)
  separation_method <- extra_args$separation_method %||% "augment"
  extra_args$separation_method <- NULL
  if (!separation_method %in% c("augment", "firth")) {
    stop("LR separation_method must be one of 'augment' or 'firth'.")
  }
  
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
  formula <- build_model_formula(
    expVar = expVar,
    coVars = coVars,
    random_effects = random_effects,
    response = "expr"
  )
  coef_name <- get_exp_coef_name(metadata, expVar, coVars)
  has_random_effects <- !is.null(random_effects) && length(random_effects) > 0L

  if (has_random_effects && !requireNamespace("glmmTMB", quietly = TRUE)) {
    stop("glmmTMB is required for LR with random_effects.")
  }
  if (has_random_effects && identical(separation_method, "firth")) {
    stop("LR separation_method = 'firth' is only supported without random_effects.")
  }
  
  ##################
  # Per-feature LR #
  ##################
  
  lr_stats <- vapply(seq_len(ncol(features)), function(j) {
    expr <- as.integer(features[, j] > 0)
    df   <- cbind(metadata, expr = expr)
    if (length(unique(expr)) < 2L) {
      return(c(coef = NA_real_, pval = NA_real_))
    }
    fit  <- try(
      if (has_random_effects || identical(separation_method, "augment")) {
        do.call(
          DA_prevalence_fit_augmented,
          c(
            list(
              formula = formula,
              df = df,
              has_random_effects = has_random_effects,
              offset = log_offset
            ),
            extra_args
          )
        )
      } else {
        do.call(
          DA_prevalence_fit_firth,
          c(
            list(
              formula = formula,
              df = df,
              offset = log_offset
            ),
            extra_args
          )
        )
      },
      silent = TRUE
    )
    if (!inherits(fit, "try-error")) {
      sm <- if (has_random_effects) {
        summary(fit)$coefficients$cond
      } else {
        coef(summary(fit))
      }
      if (coef_name %in% rownames(sm)) {
        p_col <- intersect(c("Pr(>|z|)", "Pr(>|t|)"), colnames(sm))[1]
        if (!is.na(p_col)) {
          return(c(coef = sm[coef_name, "Estimate"], pval = sm[coef_name, p_col]))
        }
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

DA_fit_enhancer_KS <- function(features, metadata, expVar, coVars = NULL, ...) {

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
    tst <- try(do.call(stats::ks.test, c(list(x = x1, y = x2), list(...))), silent = TRUE)
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

DA_fit_enhancer_WLX <- function(features, metadata, expVar, coVars = NULL, ...) {
 
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
    tst <- try(do.call(stats::wilcox.test, c(list(x = x1, y = x2), list(...))), silent = TRUE)
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
