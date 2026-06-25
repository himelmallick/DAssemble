make_internal_toy_data <- function() {
  features <- data.frame(
    taxon1 = c(20, 22, 19, 21, 2, 3, 1, 2),
    taxon2 = c(1, 2, 1, 2, 18, 20, 21, 19),
    taxon3 = c(8, 8, 7, 9, 8, 7, 9, 8),
    taxon4 = c(0, 1, 0, 1, 0, 0, 1, 0),
    check.names = FALSE
  )
  rownames(features) <- paste0("sample", seq_len(nrow(features)))
  metadata <- data.frame(
    group = factor(rep(c("A", "B"), each = 4)),
    row.names = rownames(features)
  )
  list(features = features, metadata = metadata)
}

test_that("result formatting and feature mapping helpers are stable", {
  formatted <- DAssemble:::DA_format_result(
    feature = c("tax a", "tax-b"),
    expVar = "group",
    pval_core = c(0.01, 0.2),
    coef_core = c(1.5, -0.5)
  )

  expect_named(formatted, c("feature", "metadata", "pval_core", "coef_core"))
  expect_equal(formatted$metadata, c("group", "group"))
  expect_equal(formatted$feature, c("tax a", "tax-b"))

  features <- data.frame(
    "tax a" = 1:2,
    "tax-b" = 3:4,
    check.names = FALSE
  )
  feature_map <- DAssemble:::make_feature_join_map(features)
  normalized <- DAssemble:::normalize_result_features(
    data.frame(feature = c("tax.a", "tax.b"), pval_core = c(0.1, 0.2)),
    feature_map
  )
  expect_equal(normalized$feature, c("tax a", "tax-b"))

  duplicate_features <- data.frame(
    "tax-a" = 1:2,
    "tax.a" = 3:4,
    check.names = FALSE
  )
  expect_error(
    DAssemble:::make_feature_join_map(duplicate_features),
    "non-unique"
  )

  expect_error(
    DAssemble:::normalize_result_features(data.frame(p = 0.1), feature_map),
    "feature"
  )
  expect_error(
    DAssemble:::normalize_result_features(
      data.frame(feature = "missing", pval_core = 0.1),
      feature_map
    ),
    "Failed to map"
  )
})

test_that("formula and dispatch helpers handle covariates and missing arguments", {
  metadata <- data.frame(
    group = factor(c("A", "A", "B", "B")),
    batch = factor(c("x", "y", "x", "y")),
    age = c(30, 40, 35, 45)
  )

  expect_equal(DAssemble:::build_rhs("group"), "group")
  expect_equal(DAssemble:::build_rhs("group", c("batch", "age")), "group + batch + age")

  design <- DAssemble:::build_design_matrix(metadata, "group", c("batch", "age"))
  expect_true("(Intercept)" %in% colnames(design))
  expect_true(any(grepl("^group", colnames(design))))
  expect_true(any(grepl("^batch", colnames(design))))

  expect_equal(DAssemble:::get_exp_coef_name(metadata, "group", "batch"), colnames(design)[2])

  covariates <- DAssemble:::build_covariate_matrix(metadata, c("batch", "age"))
  expect_false("(Intercept)" %in% colnames(covariates))
  expect_equal(nrow(covariates), nrow(metadata))
  expect_null(DAssemble:::build_covariate_matrix(metadata, NULL))

  fn <- function(a, b = NULL, c = NULL) {
    list(a = a, b = b, c = c)
  }
  out <- DAssemble:::DA_do_call(fn, list(a = 1, extra = 99))
  expect_equal(out, list(a = 1, b = NULL, c = NULL))
})

test_that("normalization helper validates samples and optional methods", {
  toy <- make_internal_toy_data()
  bad_metadata <- toy$metadata
  rownames(bad_metadata) <- paste0("bad", seq_len(nrow(bad_metadata)))

  expect_error(
    DAssemble:::DAssemble_normalize_features(toy$features, bad_metadata, method = "TSS"),
    "identical rownames"
  )
  expect_error(
    DAssemble:::DAssemble_normalize_features(toy$features, toy$metadata, method = "bad"),
    "should be one of"
  )

  skip_if_not_installed("edgeR")
  tmm <- DAssemble:::DAssemble_normalize_features(
    toy$features,
    toy$metadata,
    method = "TMM",
    return_normalized = TRUE
  )
  expect_s3_class(tmm$features, "data.frame")
  expect_equal(dim(tmm$features), dim(toy$features))
})

test_that("CCT and subensemble helpers cover error and combination paths", {
  expect_error(DAssemble:::CCT(c(0.1, 0.2), weights = 1), "length of weights")
  expect_error(DAssemble:::CCT(c(0.1, 0.2), weights = c(1, -1)), "positive")
  expect_warning(
    expect_true(DAssemble:::CCT(c(NA, 0.2)) >= 0),
    "exactly 1"
  )

  res <- data.frame(
    feature = c("tax1", "tax2"),
    Pval_core = c(0.01, 0.8),
    Pval_WLX = c(0.02, 0.7),
    Pval_LR = c(0.03, 0.6),
    Pval_KS = c(0.04, 0.5)
  )
  ensembles <- DAssemble:::DA_build_subensembles(
    res,
    core_method = "DESeq2",
    enh_list = c("WLX", "LR", "KS"),
    p_adj = "BH"
  )
  expect_length(ensembles, 15)
  expect_true(all(c("DESeq2", "WLX", "DESeq2+WLX+LR+KS") %in% names(ensembles)))
  expect_equal(ensembles$WLX$Pval, res$Pval_WLX)
  expect_true(all(ensembles$`DESeq2+WLX`$Pval >= 0 & ensembles$`DESeq2+WLX`$Pval <= 1))

  expect_equal(
    DAssemble:::DA_build_subensembles(data.frame(feature = "tax1"), NULL, character()),
    list()
  )
  expect_error(
    DAssemble:::DA_build_subensembles(data.frame(Pval_WLX = 0.1), NULL, "WLX"),
    "feature"
  )
})

test_that("SummarizedExperiment assay extraction returns selected assay", {
  skip_if_not_installed("SummarizedExperiment")

  counts <- matrix(1:6, nrow = 2)
  rownames(counts) <- c("gene1", "gene2")
  colnames(counts) <- paste0("sample", 1:3)
  shifted <- counts + 10
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, shifted = shifted)
  )

  expect_equal(DAssemble:::extractAssay(se, "counts"), as.data.frame(counts))
  expect_equal(DAssemble:::extractAssay(se, "shifted"), as.data.frame(shifted))
  expect_null(DAssemble:::extractAssay(se, "missing"))
})

test_that("Tweedieverse utility functions handle small deterministic inputs", {
  expect_equal(DAssemble:::entropy(c("A", "A", "A")), 0)
  expect_equal(round(DAssemble:::entropy(c("A", "A", "B", "B")), 6), 1)

  set.seed(1)
  df <- data.frame(
    taxon = paste0("tax", 1:5),
    metadata = c("group", "group", "group", "batch", "batch"),
    effect_size = c(0.5, 0.2, -0.1, 1.0, 0.8),
    pval = c(0.01, 0.2, 0.8, 0.99, 0.99),
    stderr = rep(0.2, 5),
    qval = c(0.02, 0.3, 0.9, 0.99, 0.99)
  )
  adjusted <- DAssemble:::median_comparison_tweedie(
    df,
    p_cutoff = 0.95,
    subtract_median = TRUE,
    n_sims = 25
  )

  expect_named(adjusted, c(names(df), "pval_median", "coef_median"))
  expect_equal(adjusted$taxon, df$taxon)
  expect_true(all(!is.na(adjusted$pval_median[adjusted$metadata == "group"])))
  expect_true(all(is.na(adjusted$pval_median[adjusted$metadata == "batch"])))
  expect_equal(
    adjusted$coef_median[adjusted$metadata == "group"],
    df$effect_size[df$metadata == "group"] - stats::median(df$effect_size[1:3])
  )
})

test_that("Robseq lightweight helpers validate normalization and model fallback paths", {
  features <- data.frame(
    sample1 = c(10, 5, 0),
    sample2 = c(20, 10, 2),
    sample3 = c(30, 15, 3),
    sample4 = c(40, 20, 4),
    row.names = paste0("gene", 1:3)
  )
  metadata <- data.frame(
    Exposure = factor(c("A", "A", "B", "B")),
    row.names = colnames(features)
  )

  upper <- DAssemble:::DAssemble_robseq_upper_quartile_norm(features)
  expect_s3_class(upper, "data.frame")
  expect_equal(dim(upper), dim(features))

  expect_error(
    DAssemble:::DAssemble_robust_dge(features, metadata, norm.method = "bad"),
    "Unknown Robseq normalization method"
  )

  skip_if_not_installed("MASS")
  skip_if_not_installed("dfadjust")
  model_out <- DAssemble:::DAssemble_robseq_per_gene_mod(
    expr = c(1, 2, 3, 4),
    formula = stats::as.formula("expr ~ Exposure"),
    regData = metadata,
    expVar = "Exposure"
  )
  expect_named(model_out, c("log2FC", "SE", "L.CI", "U.CI", "Pval"))
})

test_that("DAssemble validates MultiAssayExperiment metadata and feature names", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")
  toy <- make_internal_toy_data()

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = t(as.matrix(toy$features)))
  )
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(counts = se),
    colData = S4Vectors::DataFrame(toy$metadata)
  )

  expect_error(
    DAssemble(
      features = mae,
      metadata = toy$metadata,
      enhancers = "WLX",
      expVar = "group"
    ),
    "Do not provide `metadata`"
  )
  expect_error(
    DAssemble(
      features = mae,
      assay_name = "missing",
      enhancers = "WLX",
      expVar = "group"
    ),
    "not found"
  )

  duplicate_features <- toy$features
  colnames(duplicate_features)[1:2] <- c("tax-a", "tax.a")
  expect_error(
    DAssemble(
      features = duplicate_features,
      metadata = toy$metadata,
      enhancers = "WLX",
      expVar = "group"
    ),
    "non-unique"
  )
})
