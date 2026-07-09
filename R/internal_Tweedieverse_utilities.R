get_AICtab<-function(fit){
  
  ########################
  # Flag invalid options #
  ########################
  
  if (!inherits(fit, c("cpglm", "zcpglm", "glmmTMB"))) {
    stop('Not supported. Valid options are cplm , zcpglm, and glmmTMB')
  }
  
  ######################
  #  Initialize AICtab #
  ######################
  
  AICtab<-rep(NA, 5)
  
  ###########################
  # Case-by-Case Extraction #
  ###########################
  
  if (inherits(fit, "cpglm")) {
    
    ##########################################
    # Back calculate logLik and BIC from AIC #
    ##########################################
    
    AIC<-fit$aic
    AIC_multiplier<-length(fit$y) - fit$df.residual
    logLik<-(AIC - 2*AIC_multiplier)/2
    BIC_multiplier<-AIC_multiplier*log(length(fit$y))
    BIC<-BIC_multiplier + 2*logLik
    deviance<-fit$deviance
    df.resid<-fit$df.residual
    
    # Coherent output
    AICtab<-c(AIC, BIC, logLik, deviance, df.resid)
  }
  
  if (inherits(fit, "zcpglm")) {
    
    ##########################################
    # Back calculate AIC and BIC from logLik #
    ##########################################
    
    logLik<--fit$llik
    AIC_multiplier<-length(fit$y) - fit$df.residual
    BIC_multiplier<-AIC_multiplier*log(length(fit$y))
    AIC<-2*AIC_multiplier + 2*logLik
    BIC<-BIC_multiplier + 2*logLik
    deviance<-NA
    df.resid<-fit$df.residual
    
    # Coherent output
    AICtab<-c(AIC, BIC, logLik, deviance, df.resid)
    
  }
  
  if (inherits(fit, "glmmTMB")) {
    
    #######################################
    # Extract AICtab from glmmTMB objects #
    #######################################
    
    AICtab<-summary(fit)["AICtab"]$AICtab
    
  }
  
  ##########
  # Return #
  ##########
  
  names(AICtab)<-c('AIC', 'BIC', 'logLik', 'deviance', 'df.resid')
  return(AICtab)
}




# Adapted form: https://rstudio-pubs-static.s3.amazonaws.com/455435_30729e265f7a4d049400d03a18e218db.html

#' Entropy of a Vector
#'
#' Compute Shannon entropy for a vector.
#'
#' @param target A vector of values.
#' @return A numeric entropy value.
#' @examples
#' entropy(c("A", "A", "B", "B", "C"))
#' @noRd
entropy <- function(target) {
  freq <- table(target)/length(target)
  # vectorize
  vec <- as.data.frame(freq)[,2]
  #drop 0 to avoid NaN resulting from log2
  vec<-vec[vec>0]
  #compute entropy
  -sum(vec * log2(vec))
}

# Written by Grace
read_input_table <- function(path) {
  utils::read.delim(
    path,
    header = TRUE,
    sep = "\t",
    fill = TRUE,
    comment.char = "",
    check.names = FALSE,
    row.names = 1
  )
}

zcpglm_fit <- function(...) {
  zcpglm_internal <- get("zcpglm", asNamespace("cplm"))
  zcpglm_internal(...)
}

extractAssay <- function(input, assay_name = "counts") {
  
  # Extract assay name based on the user input
  if (assay_name %in% SummarizedExperiment::assayNames(input)) {
    counts_data <- SummarizedExperiment::assay(input, assay_name)
    return(as.data.frame(as.matrix(counts_data)))
  } else {
    return(NULL)
  }
}


#' Median Comparison for Compositionality Adjustment
#'
#' Adjust Tweedieverse(or any other differential analysis methods) coefficient estimates and p-values by testing each taxon
#' against the *median* effect for the same metadata variable - a simple
#' post-hoc strategy to curb false discoveries driven by the compositional
#' nature of microbiome count data (after the approach adopted in **MaAsLin 3**).
#'
#' @param df A `data.frame` returned by **Tweedieverse** or any other differential analysis methods containing (at
#'   minimum) the columns  
#'   `taxon`, `metadata`, `effect_size`, `pval`, `stderr`, and `qval`.
#' @param p_cutoff Numeric.  Upper bound on the original p-value to include an
#'   effect in the median calculation (default `0.95` = all non-missing).
#' @param subtract_median Logical.  If `TRUE`, subtracts the group
#'   median from every coefficient before returning it in `coef_median`;
#'   otherwise the original `effect_size` is copied unchanged
#'   (default `FALSE`).
#' @param n_sims Integer.  Number of Monte-Carlo simulations used to estimate
#'   the covariance between each coefficient and the group median
#'   (default `10 000`).
#' @param median_threshold Numeric.  Absolute difference below which a
#'   coefficient is considered effectively equal to the median and assigned
#'   `pval_median = 1` (default `0`).
#'
#' @details
#' For every distinct value in `metadata` the algorithm
#' \enumerate{
#'   \item keeps coefficients with `pval < p_cutoff` and computes their median;
#'   \item simulates `n_sims` draws of coefficients using a normal
#'         approximation (`N(effect_size, stderr^2)`) and records the empirical
#'         distribution of the simulated medians;
#'   \item derives a variance-inflated *offset* that accounts for the
#'         covariance between each coefficient and the group median;
#'   \item performs a two-sided Z-test of `H0 : beta_i = offset_i`, returning the
#'         resulting p-value in `pval_median`.
#' }
#'
#' @return
#' The input `df` with two new columns:
#' \describe{
#'   \item{`coef_median`}{Median-centred coefficient (or the original
#'   `effect_size` if `subtract_median = FALSE`).}
#'   \item{`pval_median`}{Two-sided p-value from the median-comparison test.}
#' }
#'
#' @section Warning:
#' The procedure is heuristic and relies on normal approximations.  Results
#' may be unstable for very small sample sizes or when `stderr` values are
#' zero or missing.
#'
#' @seealso [MaAsLin 3 GitHub](https://github.com/biobakery/maaslin3)
#'
#' @examples
#' toy <- data.frame(
#'   taxon = c("tax1", "tax2"),
#'   metadata = c("grp", "grp"),
#'   effect_size = c(0.4, 0.2),
#'   pval = c(0.01, 0.2),
#'   stderr = c(0.1, 0.1),
#'   qval = c(0.02, 0.25)
#' )
#' median_comparison_tweedie(toy, n_sims = 100)
#'
#' \dontrun{
#' 
#' ######################
#' # HMP2 input_features Analysis #
#' ######################
#'
#' #############
#' # Load input_features #
#' #############
#' 
#' library(data.table)
#' input_features <- fread("https://raw.githubusercontent.com/biobakery/Maaslin2/master/inst/extdata/HMP2_taxonomy.tsv", sep ="\t")
#' input_metadata <-fread("https://raw.githubusercontent.com/biobakery/Maaslin2/master/inst/extdata/HMP2_metadata.tsv", sep ="\t")
#'
#' ###############
#' # Format data #
#' ###############
#'
#' library(tibble)
#' features<- column_to_rownames(input_features, 'ID')
#' metadata<- column_to_rownames(input_metadata, 'ID')
#'
#' #############
#' # Fit Model #
#' #############
#'
#' HMP2 <- DAssemble_Tweedieverse(
#' features,
#' metadata,
#' output = './demo_output/HMP2', # Assuming demo_output exists
#' fixed_effects = c('diagnosis', 'dysbiosisnonIBD','dysbiosisUC','dysbiosisCD', 'antibiotics', 'age'),
#' random_effects = c('site', 'subject'),
#' base_model = 'CPLM',
#' adjust_offset = FALSE, # No offset as the values are relative abundances
#' cores = 8, # Make sure your computer has the capability
#' median_comparison = TRUE,
#' median_subtraction = TRUE,
#' standardize = FALSE,
#' reference = c('diagnosis,nonIBD'))
#' 
#' HMP2_adj <- median_comparison_tweedie(HMP2,
#'                                         p_cutoff = 0.95,
#'                                         subtract_median = TRUE,
#'                                         n_sims = 10000,
#'                                         median_threshold = 0)
#'
#' head(HMP2_adj[, c("taxon", "metadata", "coef_median", "pval_median")])
#' 
#' }
#'
#' @noRd
median_comparison_tweedie <- function(df,
                                      p_cutoff = 0.95,
                                      subtract_median = FALSE,
                                      n_sims = 10000,
                                      median_threshold = 0) {
  # df is your Tweedieverse output data.frame with columns:
  #   taxon, metadata, effect_size, pval, stderr, qval
  #
  # We'll store results here:
  df$pval_median <- NA_real_
  df$coef_median <- df$effect_size  # By default, same as effect_size
  
  # Process each metadata variable separately
  df$.orig_order <- seq_len(nrow(df))
  split_indices <- split(seq_len(nrow(df)), df$metadata)
  processed_groups <- lapply(split_indices, function(sub_idx) {
    # 1) Subset to just this metadata predictor
    sub_df  <- df[sub_idx, , drop = FALSE]

    # 2) Filter out obviously "bad" or huge p-values before computing the median
    use_idx <- which(!is.na(sub_df$pval) & sub_df$pval < p_cutoff)
    if (length(use_idx) == 0) {
      # If none are usable, move on
      return(sub_df)
    }

    # 3) Compute the "group-wide" median of the usable coefficients
    cur_median <- stats::median(sub_df$effect_size[use_idx], na.rm = TRUE)
    if (is.na(cur_median)) {
      # If no valid median, skip
      return(sub_df)
    }
    
    # 4) Optionally shift each coefficient by the median
    if (subtract_median) {
      sub_df$coef_median <- sub_df$effect_size - cur_median
    } else {
      sub_df$coef_median <- sub_df$effect_size
    }
    
    coefs    <- sub_df$effect_size
    ses      <- sub_df$stderr
    n_coefs  <- length(coefs)
    
    # Identify which coefficients were used for the median
    use_bool <- rep(FALSE, n_coefs)
    use_bool[use_idx] <- TRUE
    
    # 5) Simulate draws to approximate correlation of each coefficient w/ median
    #    sim_results has columns = draws, row 1 = simulated median, next rows = coefs
    sim_results <- replicate(n_sims, {
      sim_coefs   <- stats::rnorm(n_coefs, mean = coefs, sd = ses)
      sim_median  <- stats::median(sim_coefs[use_bool])
      c(sim_median, sim_coefs)
    })
    
    sim_medians <- sim_results[1, ]
    all_sims    <- sim_results[-1, , drop = FALSE]  # row i => draws for coef i
    
    # Covariance of each coefficient with the median, across draws
    cov_adjust <- apply(all_sims, 1, function(x) stats::cov(x, sim_medians))
    
    # 6) "offset to test" for each coefficient, per the MaAsLin 3 approach:
    #    offset_i = coefs[i] +/- ... depends on difference from median & correlation
    median_sd <- sd(sim_medians)  # the empirical SD of the simulated median
    offsets_to_test <- abs(cur_median - coefs) *
      sqrt( (ses^2) / ( ses^2 + median_sd^2 - 2*cov_adjust ) ) + coefs
    
    # 7) For each coefficient, finalize the p-value vs. the median
    #    - If difference from median < threshold => p=1
    #    - Else do test:  H0: coef[i] = offsets_to_test[i]
    #        =>  z = (coef - offset) / SE
    pvals_median <- vapply(seq_len(n_coefs), function(i) {
      if (abs(coefs[i] - cur_median) < median_threshold) {
        # If difference is trivially small => p=1
        1
      } else if (is.na(coefs[i]) || is.na(ses[i]) || ses[i] == 0) {
        NA_real_
      } else {
        # Normal approx. test for H0: coefs[i] == offsets_to_test[i]
        z_stat <- (coefs[i] - offsets_to_test[i]) / ses[i]
        2 * stats::pnorm(abs(z_stat), lower.tail = FALSE)
      }
    }, numeric(1))

    # Save the results back into the subset
    sub_df$pval_median <- pvals_median
    sub_df
  })
  df <- do.call(rbind, processed_groups)
  df <- df[order(df$.orig_order), , drop = FALSE]
  df$.orig_order <- NULL
  rownames(df) <- NULL

  # Return the augmented data
  return(df)
}
