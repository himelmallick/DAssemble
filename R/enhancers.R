###############################################
# DAssemble Enhancer Logistic Regression (LR) #
###############################################

DA_fit_enhancer_LR <- function(features, metadata, expVar) {
  
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
  formula <- as.formula(paste("expr ~", expVar))
  
  ##################
  # Per-feature LR #
  ##################
  
  pval_LR <- rep(NA_real_, ncol(features))
  names(pval_LR) <- colnames(features)
  
  for (j in seq_len(ncol(features))) {
    expr <- as.integer(features[, j] > 0)
    df   <- cbind(metadata, expr = expr)
    fit <- try(glm(formula = formula, family = binomial(), data = df, offset  = log_offset), silent = TRUE)
    
    if (!inherits(fit, "try-error")) {
      pval_LR[j] <- coef(summary(fit))[2, 4]   # p-value for expVar
    }
  }
  
  ###########################################
  # Standardized output - Feature + Pvalues #
  ###########################################
  
  feature<-names(pval_LR)
  df<-cbind.data.frame(feature, pval_LR)
  df$metadata<- expVar
  df<-dplyr::select(df, c('feature', 'metadata'), everything())
  return(df)
}


##############################################
# DAssemble Enhancer Kolmogorov–Smirnov (KS) #
##############################################

DA_fit_enhancer_KS <- function(features, metadata, expVar) {

  ########################
  # Standard KS pipeline #
  ########################
  
  group <- metadata[[expVar]]
  g1 <- which(group == levels(group)[1])
  g2 <- which(group == levels(group)[2])
  
  ##################
  # Per-feature KS #
  ##################
  
  pval_KS <- rep(NA_real_, ncol(features))
  names(pval_KS)<-colnames(features)
  for (j in seq_len(ncol(features))) {
    x1 <- features[g1, j]
    x2 <- features[g2, j]
    tst <- try(stats::ks.test(x1, x2), silent = TRUE)
    if (inherits(tst, "try-error")) {
      pval_KS[j] <- NA_real_
    } else {
      pval_KS[j] <- tst$p.value
    }
  }
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  feature<-names(pval_KS)
  df<-cbind.data.frame(feature, pval_KS)
  df$metadata<- expVar
  df<-dplyr::select(df, c('feature', 'metadata'), everything())
  return(df)
}


##############################################
# DAssemble Enhancer Wilcoxon Rank Sum (WLX) #
##############################################

DA_fit_enhancer_WLX <- function(features, metadata, expVar) {
 
  #########################
  # Standard WLX pipeline #
  #########################
  
  group <- metadata[[expVar]]
  g1 <- which(group == levels(group)[1])
  g2 <- which(group == levels(group)[2])
  
  ###################
  # Per-feature WLX #
  ###################
  
  pval_WLX <- rep(NA_real_, ncol(features))
  names(pval_WLX)<-colnames(features)
  for (j in seq_len(ncol(features))) {
    x1 <- features[g1, j]
    x2 <- features[g2, j]
    tst <- try(stats::wilcox.test(x1, x2), silent = TRUE)
    if (inherits(tst, "try-error")) {
      pvals[j] <- NA_real_
    } else {
      pval_WLX[j] <- tst$p.value
    }
  }
  
  #####################################################
  # Standardized output - Feature + Metadata + Pvalue #
  #####################################################
  
  feature<-names(pval_WLX)
  df<-cbind.data.frame(feature, pval_WLX)
  df$metadata<- expVar
  df<-dplyr::select(df, c('feature', 'metadata'), everything())
  return(df)
}