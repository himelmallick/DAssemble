library(DAssemble)

# Locate the data file in the installed package, or next to this script when
# running from a source checkout.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
data_candidates <- c(
  system.file("extdata", "hmp2_adult.RData", package = "DAssemble"),
  if (length(script_arg)) {
    file.path(
      dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]))),
      "hmp2_adult.RData"
    )
  },
  file.path(getwd(), "inst", "extdata", "hmp2_adult.RData")
)
data_file <- data_candidates[file.exists(data_candidates)][1L]
if (is.na(data_file)) {
  stop("Could not find hmp2_adult.RData.")
}

load(data_file)

if (!exists("adult_taxa", inherits = FALSE) ||
    !exists("adult_meta", inherits = FALSE)) {
  stop("hmp2_adult.RData must contain adult_taxa and adult_meta.")
}
if (!identical(rownames(adult_taxa), rownames(adult_meta))) {
  stop("The taxa and metadata sample identifiers must be identical and ordered.")
}

binary_factor <- function(x) {
  factor(as.integer(x), levels = c(0, 1))
}

# Prepare the variables used by both the UC and CD models.
adult_meta$diagnosis_CD <- binary_factor(adult_meta$diagnosis == "CD")
adult_meta$diagnosis_UC <- binary_factor(adult_meta$diagnosis == "UC")
adult_meta$dysbiosis_CD <- binary_factor(
  adult_meta$dysbiosis_state == "dysbiosis_CD"
)
adult_meta$dysbiosis_UC <- binary_factor(
  adult_meta$dysbiosis_state == "dysbiosis_UC"
)
adult_meta$dysbiosis_nonIBD <- binary_factor(
  adult_meta$dysbiosis_state == "dysbiosis_nonIBD"
)
adult_meta$abx_used <- binary_factor(adult_meta$Antibiotics == "Yes")
adult_meta$participant_id <- factor(adult_meta$`Participant ID`)

.hmp2_model_spec <- function(disease = c("UC", "CD")) {
  disease <- match.arg(disease)
  other_disease <- if (disease == "UC") "CD" else "UC"

  list(
    exposure = paste0("dysbiosis_", disease),
    covariates = c(
      "diagnosis_CD",
      "diagnosis_UC",
      paste0("dysbiosis_", other_disease),
      "dysbiosis_nonIBD",
      "abx_used",
      "reads_filtered"
    )
  )
}

run_hmp2_dassemble <- function(disease = c("UC", "CD")) {
  spec <- .hmp2_model_spec(disease)

  suppressWarnings(DAssemble::DAssemble(
    features = adult_taxa,
    metadata = adult_meta,
    expVar = spec$exposure,
    coVars = spec$covariates,
    random_effects = "participant_id",
    core_method = "Maaslin2",
    enhancers = "LR",
    method_args = list(
      Maaslin2 = list(median_comparison = TRUE)
    ),
    enhancer_norm = "TSS",
    p_adj = "BH",
    return_components = TRUE
  ))
}

# Change this value to "CD" to run the dysbiosis CD versus none analysis.
# The complete 1,667-taxon model can take substantial time to finish.
if (sys.nframe() == 0L) {
  disease <- "UC"
  fit <- run_hmp2_dassemble(disease)

  ensemble_results <- fit$res
  maaslin2_results <- fit$components$Maaslin2
  lr_results <- fit$components$LR

  print(utils::head(ensemble_results))
  print(utils::head(maaslin2_results))
  print(utils::head(lr_results))
}
