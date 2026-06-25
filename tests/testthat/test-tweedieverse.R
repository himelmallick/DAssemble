make_tweedie_toy_data <- function() {
  features <- data.frame(
    taxon1 = c(10, 12, 9, 11, 2, 1),
    taxon2 = c(1, 2, 1, 8, 9, 10),
    taxon3 = c(4, 4, 5, 4, 5, 4),
    check.names = FALSE
  )
  rownames(features) <- paste0("sample", seq_len(nrow(features)))
  metadata <- data.frame(
    group = factor(rep(c("A", "B"), each = 3)),
    batch = factor(c("x", "y", "x", "y", "x", "y")),
    age = c(30, 31, 32, 40, 41, 42),
    row.names = rownames(features)
  )
  list(features = features, metadata = metadata)
}

mock_tweedie_fit <- function(features, metadata, ...) {
  list(results = data.frame(
    feature = colnames(features),
    metadata = "group",
    value = "B",
    coef = seq_len(ncol(features)) / 10,
    stderr = rep(0.1, ncol(features)),
    pval = seq_len(ncol(features)) / 10,
    qval = seq_len(ncol(features)) / 10,
    stringsAsFactors = FALSE
  ))
}

test_that("fit.Tweedieverse dispatches only the CPLM internal model", {
  toy <- make_tweedie_toy_data()
  metadata <- transform(toy$metadata, offset = rowSums(toy$features))

  testthat::local_mocked_bindings(
    fit.CPLM = function(features, metadata, formula, adjust_offset = TRUE, ...) {
      expect_equal(colnames(features), colnames(toy$features))
      expect_true("offset" %in% colnames(metadata))
      expect_true(inherits(formula, "formula"))
      list(results = data.frame(feature = colnames(features), qval = 1))
    },
    .package = "DAssemble"
  )

  fit <- DAssemble:::fit.Tweedieverse(
    features = toy$features,
    metadata = metadata,
    base_model = "CPLM",
    adjust_offset = TRUE
  )
  expect_named(fit, "results")

  expect_error(
    DAssemble:::fit.Tweedieverse(
      features = toy$features,
      metadata = metadata,
      base_model = "ZICP"
    ),
    "Only base_model = 'CPLM'"
  )
})

test_that("DAssemble_Tweedieverse validates inputs before model fitting", {
  toy <- make_tweedie_toy_data()

  expect_error(
    DAssemble:::DAssemble_Tweedieverse(input_features = matrix(1:4, nrow = 2)),
    "Input data of class"
  )
  expect_error(
    DAssemble:::DAssemble_Tweedieverse(input_features = toy$features),
    "input_metadata must be provided"
  )
  expect_error(
    DAssemble:::DAssemble_Tweedieverse(
      input_features = toy$features,
      input_metadata = toy$metadata,
      link = "bad"
    ),
    "Option not valid"
  )
  expect_error(
    DAssemble:::DAssemble_Tweedieverse(
      input_features = toy$features,
      input_metadata = toy$metadata,
      base_model = "ZICP"
    ),
    "Option not valid"
  )
  expect_error(
    DAssemble:::DAssemble_Tweedieverse(
      input_features = toy$features,
      input_metadata = toy$metadata,
      correction = "bad"
    ),
    "Option not valid"
  )
  expect_error(
    DAssemble:::DAssemble_Tweedieverse(
      input_features = toy$features,
      input_metadata = toy$metadata,
      prev_threshold = 2
    ),
    "outside \\[0, 1\\]"
  )
})

test_that("DAssemble_Tweedieverse preprocesses data and output with mocked fitting", {
  toy <- make_tweedie_toy_data()

  testthat::local_mocked_bindings(
    fit.Tweedieverse = function(features, metadata, formula, random_effects_formula, ...) {
      expect_equal(colnames(features), colnames(toy$features))
      expect_true("offset" %in% colnames(metadata))
      expect_true(inherits(formula, "formula"))
      expect_null(random_effects_formula)
      mock_tweedie_fit(features, metadata)
    },
    .package = "DAssemble"
  )

  res <- DAssemble:::DAssemble_Tweedieverse(
    input_features = toy$features,
    input_metadata = toy$metadata,
    output = NULL,
    fixed_effects = "group",
    max_significance = 1,
    median_comparison = FALSE,
    standardize = FALSE
  )

  expect_s3_class(res, "data.frame")
  expect_named(
    res,
    c("feature", "metadata", "value", "coef", "stderr", "pval", "qval",
      "N", "N.not.zero", "percent.zero")
  )
  expect_equal(res$feature, colnames(toy$features))
  expect_equal(res$metadata, rep("group", ncol(toy$features)))
})

test_that("DAssemble_Tweedieverse handles transposed data and random effects setup", {
  toy <- make_tweedie_toy_data()

  testthat::local_mocked_bindings(
    fit.Tweedieverse = function(features, metadata, formula, random_effects_formula, ...) {
      expect_equal(rownames(features), rownames(toy$metadata))
      expect_equal(colnames(features), colnames(toy$features))
      expect_true(inherits(random_effects_formula, "formula"))
      mock_tweedie_fit(features, metadata)
    },
    .package = "DAssemble"
  )

  res_captured <- capture_warnings(
    DAssemble:::DAssemble_Tweedieverse(
      input_features = as.data.frame(t(toy$features)),
      input_metadata = toy$metadata,
      output = NULL,
      fixed_effects = "group,batch",
      random_effects = "batch",
      reference = "group,A",
      max_significance = 1,
      median_comparison = FALSE
    )
  )
  expect_captured_warning(res_captured, "as.is")
  res <- res_captured$value
  expect_equal(nrow(res), ncol(toy$features))
})
