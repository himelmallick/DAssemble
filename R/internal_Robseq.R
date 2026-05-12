#### Internal Robseq implementation ####

DAssemble_robseq_quantile_norm <- function(features) {
  norm_y <- preprocessCore::normalize.quantiles(as.matrix(features))
  norm_y <- data.frame(norm_y)
  names(norm_y) <- names(features)
  rownames(norm_y) <- rownames(features)
  norm_y
}

DAssemble_robseq_upper_quartile_norm <- function(features) {
  quant_exp <- apply(
    as.matrix(features),
    2,
    function(x) stats::quantile(x[x > 0], 0.75)
  )
  data.frame(t(t(as.matrix(features)) / quant_exp))
}

DAssemble_robseq_tmm_norm <- function(features) {
  norm_y <- edgeR::DGEList(features)
  norm_y <- edgeR::calcNormFactors(norm_y, method = "TMM")
  as.data.frame(edgeR::cpm(norm_y, log = FALSE))
}

DAssemble_robseq_rle_norm <- function(features,
                                      metadata,
                                      coVars = NULL,
                                      expVar = "Exposure") {
  if (is.null(coVars)) {
    metadata <- metadata[, c(expVar), drop = FALSE]
  } else {
    metadata <- metadata[, c(expVar, coVars)]
  }
  formula <- stats::as.formula(
    paste("~", paste(colnames(metadata), collapse = "+"), sep = "")
  )
  x <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(
    countData = as.matrix(features),
    colData = metadata,
    design = formula
  ))
  gm_mean <- function(x, na.rm = TRUE) {
    exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
  }
  geoMeans <- apply(DESeq2::counts(x), 1, gm_mean)
  s <- DESeq2::estimateSizeFactors(x, geoMeans = geoMeans)$sizeFactor
  data.frame(t(apply(features, 1, function(x) x / s)))
}

DAssemble_robseq_per_gene_mod <- function(expr, formula, regData, expVar) {
  regData <- data.frame(expr, regData)
  fit_rlm <- tryCatch(
    suppressWarnings(MASS::rlm(formula, data = regData)),
    error = function(err) suppressWarnings(stats::lm(formula, data = regData))
  )

  fit_se <- tryCatch(
    dfadjust::dfadjustSE(fit_rlm),
    error = function(err) {
      fit_lm <- suppressWarnings(stats::lm(formula, data = regData))
      dfadjust::dfadjustSE(fit_lm)
    }
  )

  if (is.list(fit_se)) {
    rDF <- as.data.frame(cbind(
      fit_se$coefficients,
      "p-value" = 2 * stats::pt(
        -abs(fit_se$coefficients[, "Estimate"] /
               fit_se$coefficients[, "HC2 se"]),
        df = fit_se$coefficients[, "df"]
      )
    ))
    exp_rows <- grep(expVar, rownames(rDF))
    log2FC <- round(rDF[exp_rows, "Estimate"], 3)
    se <- round(rDF[exp_rows, "HC2 se"], 3)
    U.CI <- log2FC + stats::qnorm(.025) * se
    L.CI <- log2FC - stats::qnorm(.025) * se
    pval <- rDF[exp_rows, "p-value"]
  } else {
    log2FC <- NA
    se <- NA
    U.CI <- NA
    L.CI <- NA
    pval <- NA
  }

  data.frame(
    log2FC = log2FC,
    SE = se,
    L.CI = L.CI,
    U.CI = U.CI,
    Pval = pval
  )
}

DAssemble_robust_dge <- function(features,
                                 metadata,
                                 norm.method = "RLE",
                                 expVar = "Exposure",
                                 coVars = NULL,
                                 filter = FALSE,
                                 verbose = FALSE) {
  start_time <- Sys.time()
  if (is.null(coVars)) {
    regData <- metadata[, c(expVar), drop = FALSE]
  } else {
    regData <- metadata[, c(expVar, coVars)]
  }
  regData[sapply(regData, is.character)] <-
    lapply(regData[sapply(regData, is.character)], as.factor)
  formula <- stats::as.formula(
    paste("expr ~ ", paste(colnames(regData), collapse = "+"))
  )

  if (norm.method == "TMM") {
    norm_y <- suppressMessages(DAssemble_robseq_tmm_norm(features))
    norm_y <- log2(norm_y + 0.5)
  } else if (norm.method == "RLE") {
    norm_y <- suppressMessages(
      DAssemble_robseq_rle_norm(features, metadata, coVars, expVar)
    )
    norm_y <- log2(norm_y + 0.5)
  } else if (norm.method == "CPM") {
    norm_y <- data.frame(edgeR::cpm(features, log = TRUE, prior.count = 1))
  } else if (norm.method == "Quantile") {
    norm_y <- suppressMessages(DAssemble_robseq_quantile_norm(features))
    norm_y <- log2(norm_y + 0.5)
  } else if (norm.method == "UQuantile") {
    norm_y <- suppressMessages(DAssemble_robseq_upper_quartile_norm(features))
    norm_y <- log2(norm_y + 0.5)
  } else {
    stop("Unknown Robseq normalization method: ", norm.method)
  }

  if (verbose && requireNamespace("pbapply", quietly = TRUE)) {
    res <- pbapply::pbapply(norm_y, 1, function(x) {
      DAssemble_robseq_per_gene_mod(
        expr = x,
        formula = formula,
        regData = regData,
        expVar = expVar
      )
    })
  } else {
    res <- apply(norm_y, 1, function(x) {
      DAssemble_robseq_per_gene_mod(
        expr = x,
        formula = formula,
        regData = regData,
        expVar = expVar
      )
    })
  }
  res <- do.call("rbind", res)

  Genes <- rownames(norm_y)
  output <- data.frame(Genes = Genes, res)
  output$adjPval <- stats::p.adjust(output$Pval, method = "BH")
  output <- output[order(output$adjPval), ]
  stop_time <- Sys.time()
  time_min <- round(difftime(stop_time, start_time, units = "mins")[[1]], 3)
  list(
    Method = "Robseq",
    res = output,
    features = features,
    Time.min = time_min
  )
}
