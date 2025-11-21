########################################################
# DAssemble Enhancer Logistic Regression (LR) Pipeline #
########################################################

DA_fit_enhancer_LR <- function(features, metadata, expVar,
                               domain = c("bulkrnaseq",
                                          "singlecell",
                                          "microbiome",
                                          "none")) {
  domain  <- match.arg(domain)
  
  aligned  <- DA_align_features_metadata(features, metadata)
  features <- aligned$features
  metadata <- aligned$metadata
  
  group <- metadata[[expVar]]
  if (is.data.frame(group)) group <- group[[1]]
  group <- droplevels(factor(group))
  if (nlevels(group) != 2L) {
    stop("LR enhancer requires a binary expVar.")
  }
  
  mat <- as.matrix(features)
  if (nrow(mat) < ncol(mat)) mat <- t(mat)
  feat_names <- colnames(mat)
  if (is.null(feat_names)) feat_names <- paste0("feat_", seq_len(ncol(mat)))
  
  y <- group
  X <- mat
  
  pvals <- rep(NA_real_, ncol(X))
  for (j in seq_len(ncol(X))) {
    xj <- X[, j]
    df <- data.frame(y = y, x = xj)
    fit <- try(stats::glm(y ~ x, data = df, family = stats::binomial()),
               silent = TRUE)
    if (inherits(fit, "try-error")) {
      pvals[j] <- NA_real_
    } else {
      coefs <- summary(fit)$coefficients
      if ("x" %in% rownames(coefs)) {
        pvals[j] <- coefs["x", "Pr(>|z|)"]
      } else {
        pvals[j] <- NA_real_
      }
    }
  }
  
  data.frame(
    feature  = feat_names,
    pval_LR  = as.numeric(pvals),
    stringsAsFactors = FALSE
  )
}


##########################################################
# DAssemble Enhancer Kolmogorov–Smirnov (KS) Pipeline    #
##########################################################

DA_fit_enhancer_KS <- function(features, metadata, expVar,
                               domain = c("bulkrnaseq",
                                          "singlecell",
                                          "microbiome",
                                          "none")) {
  domain  <- match.arg(domain)
  
  aligned  <- DA_align_features_metadata(features, metadata)
  features <- aligned$features
  metadata <- aligned$metadata
  
  group <- metadata[[expVar]]
  if (is.data.frame(group)) group <- group[[1]]
  group <- droplevels(factor(group))
  if (nlevels(group) != 2L) {
    stop("KS enhancer requires a binary expVar.")
  }
  
  mat <- as.matrix(features)
  if (nrow(mat) < ncol(mat)) mat <- t(mat)
  feat_names <- colnames(mat)
  if (is.null(feat_names)) feat_names <- paste0("feat_", seq_len(ncol(mat)))
  
  g1 <- which(group == levels(group)[1L])
  g2 <- which(group == levels(group)[2L])
  
  pvals <- rep(NA_real_, ncol(mat))
  for (j in seq_len(ncol(mat))) {
    x1 <- mat[g1, j]
    x2 <- mat[g2, j]
    tst <- try(stats::ks.test(x1, x2), silent = TRUE)
    if (inherits(tst, "try-error")) {
      pvals[j] <- NA_real_
    } else {
      pvals[j] <- tst$p.value
    }
  }
  
  data.frame(
    feature  = feat_names,
    pval_KS  = as.numeric(pvals),
    stringsAsFactors = FALSE
  )
}


#########################################################
# DAssemble Enhancer Wilcoxon Rank-Sum (WLX) Pipeline   #
#########################################################

DA_fit_enhancer_WLX <- function(features, metadata, expVar,
                                domain = c("bulkrnaseq",
                                           "singlecell",
                                           "microbiome",
                                           "none")) {
  domain  <- match.arg(domain)
  
  aligned  <- DA_align_features_metadata(features, metadata)
  features <- aligned$features
  metadata <- aligned$metadata
  
  group <- metadata[[expVar]]
  if (is.data.frame(group)) group <- group[[1]]
  group <- droplevels(factor(group))
  if (nlevels(group) != 2L) {
    stop("WLX enhancer requires a binary expVar.")
  }
  
  mat <- as.matrix(features)
  if (nrow(mat) < ncol(mat)) mat <- t(mat)
  feat_names <- colnames(mat)
  if (is.null(feat_names)) feat_names <- paste0("feat_", seq_len(ncol(mat)))
  
  g1 <- which(group == levels(group)[1L])
  g2 <- which(group == levels(group)[2L])
  
  pvals <- rep(NA_real_, ncol(mat))
  for (j in seq_len(ncol(mat))) {
    x1 <- mat[g1, j]
    x2 <- mat[g2, j]
    tst <- try(stats::wilcox.test(x1, x2), silent = TRUE)
    if (inherits(tst, "try-error")) {
      pvals[j] <- NA_real_
    } else {
      pvals[j] <- tst$p.value
    }
  }
  
  data.frame(
    feature   = feat_names,
    pval_WLX  = as.numeric(pvals),
    stringsAsFactors = FALSE
  )
}