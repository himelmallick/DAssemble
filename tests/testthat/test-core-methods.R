make_core_toy_data <- function() {
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
    batch = factor(rep(c("x", "y"), 4)),
    row.names = rownames(features)
  )
  list(features = features, metadata = metadata)
}

expect_core_result <- function(res, feature_names) {
  expect_s3_class(res, "data.frame")
  expect_true(all(c("feature", "metadata", "pval_core") %in% names(res)))
  expect_setequal(res$feature, feature_names)
  expect_true(all(res$pval_core >= 0 & res$pval_core <= 1, na.rm = TRUE))
}

cleanup_logging_handlers <- function() {
  if (requireNamespace("logging", quietly = TRUE)) {
    for (handler in c("logging::writeToFile", "writeToFile")) {
      try(logging::removeHandler(handler), silent = TRUE)
    }
  }
}

test_that("lightweight core wrappers return standardized outputs when dependencies are installed", {
  toy <- make_core_toy_data()

  skip_if_not_installed("DESeq2")
  set.seed(1)
  deseq_features <- as.data.frame(matrix(
    stats::rnbinom(240, mu = rep(c(20, 80), each = 120), size = 5),
    nrow = 12
  ))
  colnames(deseq_features) <- paste0("gene", seq_len(ncol(deseq_features)))
  rownames(deseq_features) <- paste0("sample", seq_len(nrow(deseq_features)))
  deseq_metadata <- data.frame(
    group = factor(rep(c("A", "B"), each = 6)),
    row.names = rownames(deseq_features)
  )
  deseq <- DAssemble:::DA_fit_core_DESeq2(deseq_features, deseq_metadata, "group")
  expect_core_result(deseq, colnames(deseq_features))

  skip_if_not_installed("edgeR")
  edger <- DAssemble:::DA_fit_core_edgeR(toy$features, toy$metadata, "group")
  expect_core_result(edger, colnames(toy$features))

  skip_if_not_installed("limma")
  voom <- DAssemble:::DA_fit_core_limmaVOOM(toy$features, toy$metadata, "group")
  expect_core_result(voom, colnames(toy$features))
})

test_that("Bioconductor core wrappers return standardized outputs when dependencies are installed", {
  toy <- make_core_toy_data()

  skip_if_not_installed("ALDEx2")
  aldex <- DAssemble:::DA_fit_core_ALDEx2(toy$features, toy$metadata, "group")
  expect_core_result(aldex, colnames(toy$features))

  skip_if_not_installed("dearseq")
  dear <- DAssemble:::DA_fit_core_dearseq(toy$features, toy$metadata, "group")
  expect_core_result(dear, colnames(toy$features))

  skip_if_not_installed("MicrobiomeStat")
  linda <- DAssemble:::DA_fit_core_LinDA(toy$features, toy$metadata, "group")
  expect_core_result(linda, colnames(toy$features))
})

test_that("additional microbiome core wrappers return standardized outputs", {
  toy <- make_core_toy_data()

  skip_if_not_installed("metagenomeSeq")
  metaseq <- DAssemble:::DA_fit_core_metagenomeSeq(toy$features, toy$metadata, "group")
  expect_core_result(metaseq, colnames(toy$features))

  skip_if_not_installed("MAST")
  mast <- DAssemble:::DA_fit_core_MAST(toy$features, toy$metadata, "group")
  expect_core_result(mast, colnames(toy$features))

  skip_if_not_installed("ANCOMBC")
  ancombc_captured <- capture_warnings(
    DAssemble:::DA_fit_core_ANCOMBC2(toy$features, toy$metadata, "group")
  )
  expect_captured_warning(ancombc_captured, "large number of taxa")
  ancombc <- ancombc_captured$value
  expect_core_result(ancombc, colnames(toy$features))
})

test_that("MaAsLin core wrappers return standardized outputs when installed", {
  toy <- make_core_toy_data()

  skip_if_not_installed("Maaslin2")
  maaslin2 <- DAssemble:::DA_fit_core_Maaslin2(toy$features, toy$metadata, "group")
  expect_core_result(maaslin2, colnames(toy$features))

  cleanup_logging_handlers()

  skip_if_not_installed("maaslin3")
  maaslin3_captured <- capture_warnings(
    DAssemble:::DA_fit_core_Maaslin3(toy$features, toy$metadata, "group")
  )
  expect_captured_warning(maaslin3_captured, "perfect fit")
  maaslin3 <- maaslin3_captured$value
  expect_core_result(maaslin3, colnames(toy$features))

  cleanup_logging_handlers()
})

test_that("Tweedieverse core wrapper standardizes mocked internal output", {
  toy <- make_core_toy_data()

  testthat::local_mocked_bindings(
    DAssemble_Tweedieverse = function(features, metadata, fixed_effects, ...) {
      expect_equal(features, toy$features)
      expect_equal(metadata, toy$metadata)
      expect_equal(fixed_effects, "group")
      data.frame(
        feature = colnames(features),
        metadata = "group",
        pval = seq_along(features) / 10,
        stringsAsFactors = FALSE
      )
    },
    .package = "DAssemble"
  )

  tw <- DAssemble:::DA_fit_core_Tweedieverse(toy$features, toy$metadata, "group")
  expect_core_result(tw, colnames(toy$features))
})

test_that("model-specific wrappers validate optional dependency availability", {
  toy <- make_core_toy_data()
  checked_missing_dependency_path <- FALSE

  if (!requireNamespace("maaslin3", quietly = TRUE)) {
    checked_missing_dependency_path <- TRUE
    expect_error(
      DAssemble:::DA_fit_core_Maaslin3(toy$features, toy$metadata, "group"),
      "maaslin3 required"
    )
  }
  if (!requireNamespace("Maaslin2", quietly = TRUE)) {
    checked_missing_dependency_path <- TRUE
    expect_error(
      DAssemble:::DA_fit_core_Maaslin2(toy$features, toy$metadata, "group"),
      "Maaslin2 required"
    )
  }
  if (!requireNamespace("MAST", quietly = TRUE)) {
    checked_missing_dependency_path <- TRUE
    expect_error(
      DAssemble:::DA_fit_core_MAST(toy$features, toy$metadata, "group"),
      "MAST required"
    )
  }
  if (!requireNamespace("metagenomeSeq", quietly = TRUE)) {
    checked_missing_dependency_path <- TRUE
    expect_error(
      DAssemble:::DA_fit_core_metagenomeSeq(toy$features, toy$metadata, "group"),
      "metagenomeSeq required"
    )
  }
  expect_type(checked_missing_dependency_path, "logical")
})
