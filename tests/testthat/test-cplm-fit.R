make_cplm_toy_data <- function(n = 8) {
  features <- data.frame(
    taxon1 = seq_len(n),
    taxon2 = c(2, 1, 2, 3, 5, 8, 13, 21)[seq_len(n)],
    check.names = FALSE
  )
  rownames(features) <- paste0("sample", seq_len(n))
  metadata <- data.frame(
    group = factor(rep(c("A", "B"), each = n / 2)),
    row.names = rownames(features)
  )
  metadata$offset <- rowSums(features)
  list(features = features, metadata = metadata)
}

expect_cplm_results <- function(fit, feature_names) {
  expect_type(fit, "list")
  expect_named(fit, "results")
  expect_s3_class(fit$results, "data.frame")
  expect_true(all(c(
    "feature", "metadata", "value", "coef", "stderr", "pval",
    "base.model", "tweedie.index", "qval"
  ) %in% names(fit$results)))
  expect_setequal(fit$results$feature, feature_names)
  expect_equal(fit$results$base.model, rep("CPLM", length(feature_names)))
}

test_that("fit.CPLM covers fixed-effect Tweedie model output", {
  skip_if_not_installed("cplm")
  skip_if_not_installed("pbapply")
  toy <- make_cplm_toy_data()

  fit <- DAssemble:::fit.CPLM(
    features = toy$features,
    metadata = toy$metadata,
    formula = expr ~ group,
    cores = 1
  )

  expect_cplm_results(fit, colnames(toy$features))
  expect_true(all(fit$results$pval >= 0 & fit$results$pval <= 1, na.rm = TRUE))
  expect_equal(fit$results$metadata, rep("group", ncol(toy$features)))
})

test_that("fit.CPLM covers random-effect Tweedie model output", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("pbapply")
  toy <- make_cplm_toy_data(n = 8)
  toy$metadata$subject <- factor(rep(seq_len(4), each = 2))

  fit_captured <- capture_warnings(DAssemble:::fit.CPLM(
    features = toy$features,
    metadata = toy$metadata,
    formula = expr ~ group,
    random_effects_formula = expr ~ (1 | subject),
    cores = 1
  ))
  fit <- fit_captured$value

  expect_cplm_results(fit, colnames(toy$features))
  expect_equal(fit$results$metadata, rep("group", ncol(toy$features)))
})

test_that("fit.CPLM returns NA rows when per-feature fitting fails", {
  skip_if_not_installed("cplm")
  skip_if_not_installed("pbapply")
  toy <- make_cplm_toy_data()

  fit <- DAssemble:::fit.CPLM(
    features = toy$features["taxon1"],
    metadata = toy$metadata,
    formula = expr ~ group,
    link = "bad",
    cores = 1
  )

  expect_equal(fit$results$feature, "taxon1")
  expect_true(all(is.na(fit$results[c("coef", "stderr", "pval", "qval")])))
})
