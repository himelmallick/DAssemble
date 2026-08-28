make_toy_data <- function() {
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

test_that("DAssemble validates requested methods and inputs", {
  toy <- make_toy_data()

  expect_error(
    DAssemble(toy$features, toy$metadata),
    "Specify core_method or at least one enhancer"
  )
  expect_error(
    DAssemble(toy$features, toy$metadata, core_method = "not_a_method"),
    "Unknown core_method"
  )
  expect_error(
    DAssemble(toy$features, toy$metadata, enhancers = "NOPE"),
    "Unknown enhancers"
  )

  bad_metadata <- toy$metadata
  rownames(bad_metadata) <- paste0("other", seq_len(nrow(bad_metadata)))
  expect_error(
    DAssemble(toy$features, bad_metadata, enhancers = "WLX"),
    "identical rownames"
  )

  non_binary <- toy$metadata
  non_binary$group <- factor(c("A", "A", "B", "B", "C", "C", "C", "C"))
  expect_error(
    DAssemble(toy$features, non_binary, enhancers = "WLX"),
    "expVar must be binary"
  )

  toy$metadata$subject <- factor(rep(seq_len(4), each = 2))
  expect_error(
    DAssemble(
      toy$features,
      toy$metadata,
      core_method = "edgeR",
      random_effects = "subject"
    ),
    "longitudinal-compatible core methods"
  )
})

test_that("normalization returns expected shapes and scale factors", {
  toy <- make_toy_data()

  tss <- DAssemble:::DAssemble_normalize_features(
    toy$features,
    toy$metadata,
    method = "TSS",
    return_normalized = TRUE
  )
  expect_s3_class(tss$features, "data.frame")
  expect_equal(dim(tss$features), dim(toy$features))
  expect_equal(unname(rowSums(tss$features)), rep(1, nrow(toy$features)))

  tss_raw <- DAssemble:::DAssemble_normalize_features(
    toy$features,
    toy$metadata,
    method = "TSS",
    return_normalized = FALSE
  )
  expect_equal(tss_raw$features, toy$features)
  expect_equal(unname(tss_raw$metadata$scale_factor), unname(rowSums(toy$features)))

  clr <- DAssemble:::DAssemble_normalize_features(
    toy$features,
    toy$metadata,
    method = "CLR",
    return_normalized = TRUE
  )
  expect_equal(unname(rowMeans(as.matrix(clr$features))), rep(0, nrow(toy$features)))
})

test_that("CCT helpers handle one or multiple p-value columns", {
  expect_equal(DAssemble:::CCT(c(0, 0.5)), 0)
  expect_warning(expect_equal(DAssemble:::CCT(c(1, 0.5)), 1))
  expect_error(DAssemble:::CCT(c(-0.1, 0.5)), "between 0 and 1")
  expect_error(DAssemble:::CCT(c(0, 1)), "Cannot have both")

  one_col <- matrix(c(0.01, 0.2, 0.8), ncol = 1)
  expect_equal(DAssemble:::DA_CCT_rows(one_col), as.numeric(one_col[, 1]))

  two_col <- cbind(c(0.01, 0.2, 0.8), c(0.02, 0.3, 0.7))
  combined <- DAssemble:::DA_CCT_rows(two_col)
  expect_length(combined, nrow(two_col))
  expect_true(all(combined >= 0 & combined <= 1))
})

test_that("enhancers return standardized p-value tables", {
  toy <- make_toy_data()

  wlx <- DAssemble:::DA_fit_enhancer_WLX(toy$features, toy$metadata, "group")
  ks <- DAssemble:::DA_fit_enhancer_KS(toy$features, toy$metadata, "group")

  expect_named(wlx, c("feature", "metadata", "pval_WLX"))
  expect_named(ks, c("feature", "metadata", "pval_KS"))
  expect_equal(wlx$feature, colnames(toy$features))
  expect_equal(ks$feature, colnames(toy$features))
  expect_true(all(wlx$pval_WLX >= 0 & wlx$pval_WLX <= 1, na.rm = TRUE))
  expect_true(all(ks$pval_KS >= 0 & ks$pval_KS <= 1, na.rm = TRUE))

  lr_meta <- DAssemble:::DAssemble_normalize_features(
    toy$features,
    toy$metadata,
    method = "TSS",
    return_normalized = FALSE
  )$metadata
  lr <- DAssemble:::DA_fit_enhancer_LR(toy$features, lr_meta, "group")
  expect_true(all(c("feature", "metadata", "pval_LR") %in% names(lr)))
  expect_equal(lr$feature, colnames(toy$features))
})

test_that("enhancer-only DAssemble returns combined results and subensembles", {
  toy <- make_toy_data()

  res_captured <- capture_warnings(DAssemble(
    features = toy$features,
    metadata = toy$metadata,
    core_method = NULL,
    enhancers = c("WLX", "KS"),
    expVar = "group",
    p_adj = "BH",
    enhancer_norm = "TSS",
    return_components = TRUE,
    return_subensembles = TRUE
  ))
  expect_captured_warning(res_captured, "exactly 1")
  res <- res_captured$value

  expect_type(res, "list")
  expect_equal(res$Method, "CCT(WLX+KS)")
  expect_true(all(c("feature", "metadata", "pval_WLX", "pval_KS",
                    "pval_joint", "qval") %in% names(res$res)))
  expect_equal(nrow(res$res), ncol(toy$features))
  expect_true(all(res$res$pval_joint >= 0 & res$res$pval_joint <= 1, na.rm = TRUE))
  expect_true(all(res$res$qval >= 0 & res$res$qval <= 1, na.rm = TRUE))
  expect_named(res$components, c("WLX", "KS"))
  expect_true(all(c("WLX", "KS", "WLX+KS") %in% names(res$ensembles)))
})

test_that("single enhancer DAssemble bypasses CCT and keeps p-values", {
  toy <- make_toy_data()

  res <- DAssemble(
    features = toy$features,
    metadata = toy$metadata,
    enhancers = "WLX",
    expVar = "group",
    p_adj = "BH",
    enhancer_norm = "TSS"
  )

  expect_equal(res$res$pval_joint, res$res$pval_WLX)
  expect_equal(res$res$qval, stats::p.adjust(res$res$pval_WLX, method = "BH"))
})

test_that("DAssemble supports longitudinal LR through random_effects", {
  toy <- make_toy_data()
  toy$metadata$subject <- factor(rep(seq_len(4), each = 2))

  testthat::local_mocked_bindings(
    DA_fit_enhancer_LR = function(features,
                                  metadata,
                                  expVar,
                                  coVars = NULL,
                                  random_effects = NULL) {
      expect_equal(random_effects, "subject")
      data.frame(
        feature = colnames(features),
        metadata = expVar,
        pval_LR = rep(0.25, ncol(features)),
        stringsAsFactors = FALSE
      )
    },
    .package = "DAssemble"
  )

  res <- DAssemble(
    features = toy$features,
    metadata = toy$metadata,
    core_method = NULL,
    enhancers = "LR",
    expVar = "group",
    random_effects = "subject",
    p_adj = "BH"
  )

  expect_equal(res$res$pval_joint, res$res$pval_LR)
})

test_that("DAssemble supports longitudinal MaAsLin3 core through random_effects", {
  toy <- make_toy_data()
  toy$metadata$subject <- factor(rep(seq_len(4), each = 2))

  testthat::local_mocked_bindings(
    DA_fit_core_Maaslin3 = function(features,
                                    metadata,
                                    expVar,
                                    coVars = NULL,
                                    random_effects = NULL,
                                    prevalence_only = FALSE) {
      expect_equal(random_effects, "subject")
      expect_true(is.null(prevalence_only) || identical(prevalence_only, FALSE))
      data.frame(
        feature = colnames(features),
        metadata = expVar,
        pval_core = rep(0.15, ncol(features)),
        stringsAsFactors = FALSE
      )
    },
    .package = "DAssemble"
  )

  res <- DAssemble(
    features = toy$features,
    metadata = toy$metadata,
    core_method = "Maaslin3",
    expVar = "group",
    random_effects = "subject",
    p_adj = "BH"
  )

  expect_equal(res$res$pval_joint, res$res$pval_core)
})

test_that("DAssemble forwards method_args to core and enhancer wrappers", {
  toy <- make_toy_data()

  testthat::local_mocked_bindings(
    DA_fit_core_Maaslin2 = function(features,
                                    metadata,
                                    expVar,
                                    coVars = NULL,
                                    random_effects = NULL,
                                    normalization = "TSS",
                                    transform = "LOG",
                                    ...) {
      expect_equal(normalization, "CLR")
      expect_equal(transform, "NONE")
      data.frame(
        feature = colnames(features),
        metadata = expVar,
        pval_core = rep(0.2, ncol(features)),
        stringsAsFactors = FALSE
      )
    },
    DA_fit_enhancer_LR = function(features,
                                  metadata,
                                  expVar,
                                  coVars = NULL,
                                  random_effects = NULL,
                                  family = "default",
                                  ...) {
      expect_equal(family, "binomial")
      data.frame(
        feature = colnames(features),
        metadata = expVar,
        pval_LR = rep(0.3, ncol(features)),
        stringsAsFactors = FALSE
      )
    },
    .package = "DAssemble"
  )

  res <- DAssemble(
    features = toy$features,
    metadata = toy$metadata,
    core_method = "Maaslin2",
    enhancers = "LR",
    expVar = "group",
    method_args = list(
      core = list(transform = "LOG"),
      Maaslin2 = list(normalization = "CLR", transform = "NONE"),
      enhancer = list(family = "binomial")
    )
  )

  expect_true(all(c("pval_core", "pval_LR", "pval_joint") %in% names(res$res)))
})

test_that("LR enhancer uses augmented logistic fitting by default", {
  toy <- make_toy_data()
  called <- FALSE

  testthat::local_mocked_bindings(
    DA_prevalence_fit_augmented = function(formula, df, has_random_effects = FALSE, ...) {
      called <<- TRUE
      fit <- stats::glm(formula = formula, family = stats::binomial(), data = df)
      fit
    },
    DA_prevalence_fit_firth = function(formula, df, offset = NULL, ...) {
      stop("should not use firth by default")
    },
    .package = "DAssemble"
  )

  res <- DAssemble:::DA_fit_enhancer_LR(
    features = toy$features,
    metadata = toy$metadata,
    expVar = "group"
  )

  expect_true(called)
  expect_true("pval_LR" %in% names(res))
  expect_equal(length(res$pval_LR), ncol(toy$features))
})

test_that("LR enhancer supports opt-in firth fitting without random effects", {
  toy <- make_toy_data()
  called <- FALSE

  testthat::local_mocked_bindings(
    DA_prevalence_fit_augmented = function(formula, df, has_random_effects = FALSE, ...) {
      stop("should not use augmentation when firth is requested")
    },
    DA_prevalence_fit_firth = function(formula, df, offset = NULL, ...) {
      called <<- TRUE
      stats::glm(formula = formula, family = stats::binomial(), data = df, offset = offset)
    },
    .package = "DAssemble"
  )

  res <- DAssemble:::DA_fit_enhancer_LR(
    features = toy$features,
    metadata = toy$metadata,
    expVar = "group",
    separation_method = "firth"
  )

  expect_true(called)
  expect_true("pval_LR" %in% names(res))
})

test_that("LR firth fitting is rejected for longitudinal random effects", {
  toy <- make_toy_data()
  toy$metadata$subject <- factor(rep(seq_len(4), each = 2))

  expect_error(
    DAssemble:::DA_fit_enhancer_LR(
      features = toy$features,
      metadata = toy$metadata,
      expVar = "group",
      random_effects = "subject",
      separation_method = "firth"
    ),
    "only supported without random_effects"
  )
})

test_that("DAssemble accepts MultiAssayExperiment input", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")
  toy <- make_toy_data()

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = t(as.matrix(toy$features)))
  )
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(counts = se),
    colData = S4Vectors::DataFrame(toy$metadata)
  )

  res <- DAssemble(
    features = mae,
    enhancers = "WLX",
    expVar = "group",
    p_adj = "BH",
    enhancer_norm = "TSS"
  )

  expect_equal(nrow(res$res), ncol(toy$features))
  expect_equal(res$res$pval_joint, res$res$pval_WLX)
  expect_equal(rownames(res$features), rownames(toy$features))
})

test_that("MultiAssayExperiment input requires assay_name for multiple experiments", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")
  toy <- make_toy_data()

  se_counts <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = t(as.matrix(toy$features)))
  )
  se_other <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = t(as.matrix(toy$features + 1)))
  )
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(counts = se_counts, shifted = se_other),
    colData = S4Vectors::DataFrame(toy$metadata)
  )

  expect_error(
    DAssemble(features = mae, enhancers = "WLX", expVar = "group"),
    "Specify `assay_name`"
  )

  res <- DAssemble(
    features = mae,
    assay_name = "shifted",
    enhancers = "WLX",
    expVar = "group",
    enhancer_norm = "TSS"
  )
  expect_equal(nrow(res$res), ncol(toy$features))
})


test_that("LOCOM wrapper uses LOCOM2 and returns standardized output", {
  skip_if_not_installed("LOCOM2")
  toy <- make_toy_data()

  res <- DAssemble:::DA_fit_core_LOCOM(toy$features, toy$metadata, "group")

  expect_named(res, c("feature", "metadata", "pval_core"))
  expect_setequal(res$feature, colnames(toy$features))
  expect_equal(res$metadata, rep("group", ncol(toy$features)))
  expect_true(all(res$pval_core >= 0 & res$pval_core <= 1, na.rm = TRUE))
})

test_that("Robseq wrapper is internally available", {
  skip_if_not_installed("dfadjust")
  skip_if_not_installed("preprocessCore")
  toy <- make_toy_data()

  res <- DAssemble:::DA_fit_core_Robseq(toy$features, toy$metadata, "group")

  expect_named(res, c("feature", "metadata", "pval_core"))
  expect_setequal(res$feature, colnames(toy$features))
  expect_equal(res$metadata, rep("group", ncol(toy$features)))
})
