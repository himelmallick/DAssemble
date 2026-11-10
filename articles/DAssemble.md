# DAssemble: Ensemble Differential Analysis

## Introduction

DAssemble is an R package for ensemble-based differential association
analysis of high-throughput omics data. It provides a flexible framework
for combining a user-specified **core method** with one or more
complementary **enhancer** methods, allowing users to integrate
statistical evidence across different modeling assumptions rather than
relying on a single differential analysis approach. DAssemble aggregates
model-specific p-values using the **Cauchy Combination Tests (CCT)** and
reports both ensemble-level and method-specific results.

The package was developed to support robust and interpretable
differential analysis across multiple omics domains, including bulk
RNA-seq, single-cell RNA-seq, and microbiome data. We submitted
DAssemble to Bioconductor to make the method easier to install,
document, test, and integrate with established Bioconductor workflows
and data structures, thereby supporting reproducible and extensible
analysis of high-throughput biological data.

## Overview

`DAssemble` implements an ensemble framework for differential-abundance
and differential-expression analysis. Users can choose a single core
model and optionally combine it with enhancer tests using the Cauchy
Combination Test.

The package accepts a `MultiAssayExperiment` object directly. When
multiple experiments are present, use `assay_name` to select the
experiment to analyze. For lower-level workflows, `features` can also be
a data frame with samples in rows and features in columns, paired with a
`metadata` data frame with matching sample row names. The exposure
variable supplied through `expVar` must be binary.

The main
[`DAssemble()`](https://himelmallick.github.io/DAssemble/reference/DAssemble.md)
function also supports:

- multiple adjustment covariates through `coVars`
- repeated-measures / longitudinal designs through `random_effects`
- open-ended method-specific controls through `method_args`

These arguments are available only for methods whose wrappers can
express the corresponding model directly. For unsupported methods,
DAssemble stops with a clear error rather than fitting a mis-specified
model.

## Main API

The primary entry point is:

``` r

DAssemble(
  features,
  metadata = NULL,
  core_method = NULL,
  enhancers = NULL,
  expVar = "group",
  coVars = NULL,
  random_effects = NULL,
  method_args = NULL,
  assay_name = NULL,
  p_adj = "BY",
  enhancer_norm = "TSS",
  return_components = TRUE,
  return_subensembles = FALSE
)
```

`method_args` is a named control list used to forward additional
arguments to underlying methods. When an argument overlaps with a
DAssemble wrapper default, the value in `method_args` takes precedence.
The usual pattern is:

- `method_args$core` for arguments shared by the selected core method
- `method_args$enhancer` for arguments shared by the selected enhancers
- `method_args$<MethodName>` for method-specific overrides

For the `LR` enhancer, the wrapper additionally recognizes
`method_args$LR$separation_method`. The default is `"augment"`, which
follows the MaAsLin3-style augmented logistic fit. For cross-sectional
presence/absence logistic analyses, users can instead request `"firth"`
to use bias-reduced logistic regression through `brglm2`.

For the `Maaslin2` core, the wrapper also recognizes
`method_args$Maaslin2$median_comparison = TRUE`. When requested,
DAssemble uses raw `fit$results`, applies the internal
`median_comparison_tweedie()` adjustment there, and only then reduces
the output back to the standardized exposure-specific result table.

Wrappers with multi-stage internals can also accept subcall-specific
entries. For example, `method_args$DESeq2` may contain items such as
`DESeq`, `results`, or `estimateSizeFactors`, while `method_args$edgeR`
may contain items such as `glmFit` or `glmLRT`.

``` r

res <- DAssemble(
  features = features,
  metadata = metadata,
  enhancers = "LR",
  expVar = "group",
  method_args = list(
    LR = list(separation_method = "firth")
  )
)
```

For example, MaAsLin3’s prevalence median comparison can be enabled
with:

``` r

res <- DAssemble(
  features = features,
  metadata = metadata,
  core_method = "Maaslin3",
  expVar = "group",
  method_args = list(
    Maaslin3 = list(median_comparison_prevalence = TRUE)
  )
)
```

## Supported methods

Core methods currently available in
[`DAssemble()`](https://himelmallick.github.io/DAssemble/reference/DAssemble.md)
include:

- `DESeq2`
- `edgeR`
- `limmaVOOM`
- `metagenomeSeq`
- `MAST`
- `dearseq`
- `ALDEx2`
- `LinDA`
- `LOCOM`
- `Maaslin2`
- `Maaslin3`
- `Tweedieverse`
- `Robseq`
- `ANCOMBC2`

Enhancers currently available are:

- `WLX`
- `LR`
- `KS`

## Covariate and longitudinal support

The current implementation supports multiple covariates (`coVars`) for
the following methods:

| Method        | Multiple covariates |
|:--------------|:--------------------|
| DESeq2        | Yes                 |
| edgeR         | Yes                 |
| limmaVOOM     | Yes                 |
| metagenomeSeq | Yes                 |
| MAST          | Yes                 |
| dearseq       | Yes                 |
| ALDEx2        | Yes                 |
| LinDA         | Yes                 |
| Maaslin2      | Yes                 |
| Maaslin3      | Yes                 |
| Tweedieverse  | Yes                 |
| Robseq        | Yes                 |
| ANCOMBC2      | Yes                 |
| LR            | Yes                 |
| WLX           | No                  |
| KS            | No                  |
| LOCOM         | No                  |

When `random_effects` is supplied, the following methods are currently
treated as longitudinal-compatible:

| Method       | Longitudinal support |
|:-------------|:---------------------|
| Maaslin2     | Yes                  |
| Maaslin3     | Yes                  |
| Tweedieverse | Yes                  |
| LR           | Yes                  |

All other methods should be treated as cross-sectional only in the
current package implementation.

For the `LR` enhancer, the default fitting path uses MaAsLin3-style
augmentation in both cross-sectional and longitudinal settings. The
alternative `method_args$LR$separation_method = "firth"` is available
only for the non-longitudinal model.

`Maaslin3` is used through its standard package pathway in
[`DAssemble()`](https://himelmallick.github.io/DAssemble/reference/DAssemble.md)
and is not forced into a prevalence-only mode.

## Installation

Install `DAssemble` from Bioconductor with `BiocManager`:

``` r

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("DAssemble")
```

## Bulk RNA-seq Example

The `airway` dataset from Bioconductor contains RNA-seq counts from
airway smooth muscle cells treated with dexamethasone. The example below
runs `DAssemble` with `DESeq2` as the core method and two enhancers.

``` r

data("airway", package = "airway")
counts <- SummarizedExperiment::assay(airway, "counts")
keep_genes <- order(rowSums(counts), decreasing = TRUE)[seq_len(500L)]
counts <- counts[keep_genes, , drop = FALSE]
metadata <- as.data.frame(SummarizedExperiment::colData(airway))
metadata <- metadata[colnames(counts), , drop = FALSE]

se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = counts)
)

mae <- MultiAssayExperiment::MultiAssayExperiment(
  experiments = list(rnaseq = se),
  colData = S4Vectors::DataFrame(metadata)
)

res <- DAssemble::DAssemble(
  features = mae,
  assay_name = "rnaseq",
  core_method = "DESeq2",
  enhancers = c("WLX", "LR"),
  expVar = "dex",
  coVars = "cell",
  p_adj = "BH",
  enhancer_norm = "tmm",
  method_args = list(
    DESeq2 = list(
      DESeq = list(fitType = "local")
    ),
    LR = list(
      control = glm.control(maxit = 100)
    )
  ),
  return_components = TRUE,
  return_subensembles = TRUE
)
#> using pre-existing size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
#> calcNormFactors has been renamed to normLibSizes

head(res$res)
#>           feature metadata    pval_core   pval_WLX coef_LR pval_LR pval_joint
#> 1 ENSG00000115414      dex 9.388328e-01 1.00000000      NA      NA          1
#> 2 ENSG00000011465      dex 1.523315e-03 0.48571429      NA      NA          1
#> 3 ENSG00000198804      dex 8.842089e-01 0.88571429      NA      NA          1
#> 4 ENSG00000156508      dex 4.568638e-04 0.11428571      NA      NA          1
#> 5 ENSG00000164692      dex 1.759305e-08 0.02857143      NA      NA          1
#> 6 ENSG00000116260      dex 2.836045e-01 0.68571429      NA      NA          1
#>   qval
#> 1    1
#> 2    1
#> 3    1
#> 4    1
#> 5    1
#> 6    1
names(res$components)
#> [1] "DESeq2" "WLX"    "LR"
```

This analysis adjusts for `cell`, combines DESeq2 p-values with Wilcoxon
and logistic-regression enhancer p-values, and forwards example control
arguments through `method_args`. The combined results are returned in
`res$res`, and the per-method outputs are available in `res$components`.

## Longitudinal HMP2 Example

The Human Microbiome Project 2 adult example included with the package
shows a longitudinal microbiome analysis with multiple covariates. The
feature table and metadata are stored in `inst/extdata` and can be
accessed after installation with
[`system.file()`](https://rdrr.io/r/base/system.file.html). To keep the
vignette build lightweight, the evaluated example below uses the most
abundant taxa and runs the Dysbiosis_UC scenario.

``` r

hmp2_file <- system.file(
  "extdata",
  "hmp2_adult.RData",
  package = "DAssemble",
  mustWork = TRUE
)
load(hmp2_file)

stopifnot(
  exists("adult_taxa"),
  exists("adult_meta"),
  identical(rownames(adult_taxa), rownames(adult_meta))
)

binary_factor <- function(x) {
  factor(as.integer(x), levels = c(0, 1))
}

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

disease <- "UC"
other_disease <- "CD"
hmp2_taxa <- adult_taxa[
  ,
  order(colSums(adult_taxa), decreasing = TRUE)[seq_len(20L)],
  drop = FALSE
]

res <- DAssemble::DAssemble(
  features = hmp2_taxa,
  metadata = adult_meta,
  expVar = paste0("dysbiosis_", disease),
  coVars = c(
    "diagnosis_CD",
    "diagnosis_UC",
    paste0("dysbiosis_", other_disease),
    "dysbiosis_nonIBD",
    "abx_used",
    "reads_filtered"
  ),
  random_effects = "participant_id",
  core_method = "Maaslin2",
  enhancers = "LR",
  method_args = list(
    Maaslin2 = list(median_comparison = TRUE)
  ),
  enhancer_norm = "TSS",
  p_adj = "BH",
  return_components = TRUE
)
#> [1] "Creating output feature tables folder"
#> [1] "Creating output fits folder"
#> 2026-11-10 04:17:16.902989 INFO::Writing function arguments to log file
#> 2026-11-10 04:17:16.914932 INFO::Verifying options selected are valid
#> 2026-11-10 04:17:16.938054 INFO::Determining format of input files
#> 2026-11-10 04:17:16.938847 INFO::Input format is data samples as rows and metadata samples as rows
#> 2026-11-10 04:17:16.944594 INFO::Formula for random effects: expr ~ (1 | participant_id)
#> 2026-11-10 04:17:16.945222 INFO::Formula for fixed effects: expr ~  dysbiosis_UC + diagnosis_CD + diagnosis_UC + dysbiosis_CD + dysbiosis_nonIBD + abx_used + reads_filtered
#> 2026-11-10 04:17:16.945689 INFO::Factor detected for categorial metadata 'dysbiosis_UC'. Provide a reference argument or manually set factor ordering to change reference level.
#> 2026-11-10 04:17:16.94608 INFO::Factor detected for categorial metadata 'diagnosis_CD'. Provide a reference argument or manually set factor ordering to change reference level.
#> 2026-11-10 04:17:16.946433 INFO::Factor detected for categorial metadata 'diagnosis_UC'. Provide a reference argument or manually set factor ordering to change reference level.
#> 2026-11-10 04:17:16.946768 INFO::Factor detected for categorial metadata 'dysbiosis_CD'. Provide a reference argument or manually set factor ordering to change reference level.
#> 2026-11-10 04:17:16.94715 INFO::Factor detected for categorial metadata 'dysbiosis_nonIBD'. Provide a reference argument or manually set factor ordering to change reference level.
#> 2026-11-10 04:17:16.947489 INFO::Factor detected for categorial metadata 'abx_used'. Provide a reference argument or manually set factor ordering to change reference level.
#> 2026-11-10 04:17:16.947815 INFO::Filter data based on min abundance and min prevalence
#> 2026-11-10 04:17:16.948172 INFO::Total samples in data: 1078
#> 2026-11-10 04:17:16.948502 INFO::Min samples required with min abundance for a feature not to be filtered: 107.800000
#> 2026-11-10 04:17:16.949663 INFO::Total filtered features: 0
#> 2026-11-10 04:17:16.950168 INFO::Filtered feature names from abundance and prevalence filtering:
#> 2026-11-10 04:17:16.951415 INFO::Total filtered features with variance filtering: 0
#> 2026-11-10 04:17:16.951823 INFO::Filtered feature names from variance filtering:
#> 2026-11-10 04:17:16.952192 INFO::Running selected normalization method: TSS
#> 2026-11-10 04:17:17.00817 INFO::Applying z-score to standardize continuous metadata
#> 2026-11-10 04:17:17.017607 INFO::Running selected transform method: LOG
#> 2026-11-10 04:17:17.019854 INFO::Running selected analysis method: LM
#> 2026-11-10 04:17:17.21761 INFO::Fitting model to feature number 1, UNCLASSIFIED
#> 2026-11-10 04:17:17.305012 INFO::Fitting model to feature number 2, k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae.g__Phocaeicola.s__Phocaeicola_vulgatus.t__SGB1814
#> 2026-11-10 04:17:17.351648 INFO::Fitting model to feature number 3, k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae.g__Bacteroides.s__Bacteroides_uniformis.t__SGB1836
#> 2026-11-10 04:17:17.398811 INFO::Fitting model to feature number 4, k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae.g__Phocaeicola.s__Phocaeicola_dorei.t__SGB1815
#> 2026-11-10 04:17:17.446687 INFO::Fitting model to feature number 5, k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae.g__Bacteroides.s__Bacteroides_stercoris.t__SGB1830
#> 2026-11-10 04:17:17.591343 INFO::Fitting model to feature number 6, k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Prevotellaceae.g__Prevotella.s__Prevotella_copri_clade_A.t__SGB1626
#> 2026-11-10 04:17:17.634575 INFO::Fitting model to feature number 7, k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae.g__Bacteroides.s__Bacteroides_ovatus.t__SGB1871
#> 2026-11-10 04:17:17.679068 INFO::Fitting model to feature number 8, k__Bacteria.p__Firmicutes.c__Clostridia.o__Eubacteriales.f__Oscillospiraceae.g__Faecalibacterium.s__Faecalibacterium_prausnitzii.t__SGB15342
#> 2026-11-10 04:17:17.723837 INFO::Fitting model to feature number 9, k__Bacteria.p__Firmicutes.c__Clostridia.o__Eubacteriales.f__Lachnospiraceae.g__Lachnospiraceae_unclassified.s__Eubacterium_rectale.t__SGB4933
#> 2026-11-10 04:17:17.769191 INFO::Fitting model to feature number 10, k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Rikenellaceae.g__Alistipes.s__Alistipes_putredinis.t__SGB2318
#> 2026-11-10 04:17:17.82326 INFO::Fitting model to feature number 11, k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Tannerellaceae.g__Parabacteroides.s__Parabacteroides_distasonis.t__SGB1934
#> 2026-11-10 04:17:18.560122 INFO::Fitting model to feature number 12, k__Bacteria.p__Firmicutes.c__Clostridia.o__Eubacteriales.f__Lachnospiraceae.g__Roseburia.s__Roseburia_intestinalis.t__SGB4951
#> 2026-11-10 04:17:18.604475 INFO::Fitting model to feature number 13, k__Bacteria.p__Firmicutes.c__Clostridia.o__Eubacteriales.f__Oscillospiraceae.g__Faecalibacterium.s__Faecalibacterium_prausnitzii.t__SGB15316
#> 2026-11-10 04:17:18.648463 INFO::Fitting model to feature number 14, k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae.g__Bacteroides.s__Bacteroides_fragilis.t__SGB1855
#> 2026-11-10 04:17:18.69172 INFO::Fitting model to feature number 15, k__Bacteria.p__Firmicutes.c__Clostridia.o__Eubacteriales.f__Lachnospiraceae.g__Roseburia.s__Roseburia_faecis.t__SGB4925
#> 2026-11-10 04:17:18.736925 INFO::Fitting model to feature number 16, k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae.g__Bacteroides.s__Bacteroides_caccae.t__SGB1877
#> 2026-11-10 04:17:18.782524 INFO::Fitting model to feature number 17, k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae.g__Phocaeicola.s__Phocaeicola_massiliensis.t__SGB1812
#> 2026-11-10 04:17:18.828215 INFO::Fitting model to feature number 18, k__Bacteria.p__Firmicutes.c__Clostridia.o__Eubacteriales.f__Oscillospiraceae.g__Faecalibacterium.s__Faecalibacterium_prausnitzii.t__SGB15332
#> 2026-11-10 04:17:18.874717 INFO::Fitting model to feature number 19, k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae.g__Bacteroides.s__Bacteroides_thetaiotaomicron.t__SGB1861
#> 2026-11-10 04:17:18.9202 INFO::Fitting model to feature number 20, k__Bacteria.p__Firmicutes.c__Clostridia.o__Eubacteriales.f__Oscillospiraceae.g__Faecalibacterium.s__Faecalibacterium_prausnitzii.t__SGB15318
#> 2026-11-10 04:17:18.972611 INFO::Counting total values for each feature
#> 2026-11-10 04:17:18.977453 INFO::Writing filtered data to file /tmp/RtmpQhqX3C/m2_78435735/features/filtered_data.tsv
#> 2026-11-10 04:17:18.98647 INFO::Writing filtered, normalized data to file /tmp/RtmpQhqX3C/m2_78435735/features/filtered_data_norm.tsv
#> 2026-11-10 04:17:18.995394 INFO::Writing filtered, normalized, transformed data to file /tmp/RtmpQhqX3C/m2_78435735/features/filtered_data_norm_transformed.tsv
#> 2026-11-10 04:17:19.006048 INFO::Writing residuals to file /tmp/RtmpQhqX3C/m2_78435735/fits/residuals.rds
#> 2026-11-10 04:17:19.013204 INFO::Writing fitted values to file /tmp/RtmpQhqX3C/m2_78435735/fits/fitted.rds
#> 2026-11-10 04:17:19.019592 INFO::Writing extracted random effects to file /tmp/RtmpQhqX3C/m2_78435735/fits/ranef.rds
#> 2026-11-10 04:17:19.020493 INFO::Writing all results to file (ordered by increasing q-values): /tmp/RtmpQhqX3C/m2_78435735/all_results.tsv
#> 2026-11-10 04:17:19.02188 INFO::Writing the significant results (those which are less than or equal to the threshold of 1.000000 ) to file (ordered by increasing q-values): /tmp/RtmpQhqX3C/m2_78435735/significant_results.tsv

head(res$res)
#>                                                                                                                                          feature
#> 16                                                                                                                                  UNCLASSIFIED
#> 15          k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Phocaeicola|s__Phocaeicola_massiliensis|t__SGB1812
#> 3  k__Bacteria|p__Firmicutes|c__Clostridia|o__Eubacteriales|f__Lachnospiraceae|g__Lachnospiraceae_unclassified|s__Eubacterium_rectale|t__SGB4933
#> 1               k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Phocaeicola|s__Phocaeicola_vulgatus|t__SGB1814
#> 2              k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_uniformis|t__SGB1836
#> 6                 k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_ovatus|t__SGB1871
#>        metadata  coef_core    pval_core     coef_LR      pval_LR   pval_joint
#> 16 dysbiosis_UC  1.9440197 4.197892e-10 -0.03729951 9.906647e-01 8.395783e-10
#> 15 dysbiosis_UC  2.6487670 5.932789e-07 -0.17424170 8.790036e-01 1.186563e-06
#> 3  dysbiosis_UC -2.0488947 1.422226e-02 -3.28080747 5.571692e-05 1.109993e-04
#> 1  dysbiosis_UC -1.8508152 7.016859e-03 -2.47057587 1.099657e-03 1.901367e-03
#> 2  dysbiosis_UC -1.3281813 5.366811e-02 -4.05438628 8.180096e-04 1.611677e-03
#> 6  dysbiosis_UC -0.3504128 5.555118e-01 -5.59457718 8.667992e-04 1.734418e-03
#>            qval
#> 16 1.679157e-08
#> 15 1.186563e-05
#> 3  7.399952e-04
#> 1  5.977412e-03
#> 2  5.977412e-03
#> 6  5.977412e-03
head(res$components$Maaslin2)
#>                                                                                                                                         feature
#> 1              k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Phocaeicola|s__Phocaeicola_vulgatus|t__SGB1814
#> 2             k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_uniformis|t__SGB1836
#> 3 k__Bacteria|p__Firmicutes|c__Clostridia|o__Eubacteriales|f__Lachnospiraceae|g__Lachnospiraceae_unclassified|s__Eubacterium_rectale|t__SGB4933
#> 4  k__Bacteria|p__Firmicutes|c__Clostridia|o__Eubacteriales|f__Oscillospiraceae|g__Faecalibacterium|s__Faecalibacterium_prausnitzii|t__SGB15318
#> 5                k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_caccae|t__SGB1877
#> 6                k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_ovatus|t__SGB1871
#>       metadata  coef_core   pval_core
#> 1 dysbiosis_UC -1.8508152 0.007016859
#> 2 dysbiosis_UC -1.3281813 0.053668107
#> 3 dysbiosis_UC -2.0488947 0.014222257
#> 4 dysbiosis_UC -0.7591049 0.258902795
#> 5 dysbiosis_UC -0.2689939 0.633510891
#> 6 dysbiosis_UC -0.3504128 0.555511835
head(res$components$LR)
#>                                                                                                                                                                                                                                                                 feature
#> UNCLASSIFIED                                                                                                                                                                                                                                               UNCLASSIFIED
#> k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Phocaeicola|s__Phocaeicola_vulgatus|t__SGB1814       k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Phocaeicola|s__Phocaeicola_vulgatus|t__SGB1814
#> k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_uniformis|t__SGB1836     k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_uniformis|t__SGB1836
#> k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Phocaeicola|s__Phocaeicola_dorei|t__SGB1815             k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Phocaeicola|s__Phocaeicola_dorei|t__SGB1815
#> k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_stercoris|t__SGB1830     k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_stercoris|t__SGB1830
#> k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Prevotellaceae|g__Prevotella|s__Prevotella_copri_clade_A|t__SGB1626 k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Prevotellaceae|g__Prevotella|s__Prevotella_copri_clade_A|t__SGB1626
#>                                                                                                                                         metadata
#> UNCLASSIFIED                                                                                                                        dysbiosis_UC
#> k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Phocaeicola|s__Phocaeicola_vulgatus|t__SGB1814    dysbiosis_UC
#> k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_uniformis|t__SGB1836   dysbiosis_UC
#> k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Phocaeicola|s__Phocaeicola_dorei|t__SGB1815       dysbiosis_UC
#> k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_stercoris|t__SGB1830   dysbiosis_UC
#> k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Prevotellaceae|g__Prevotella|s__Prevotella_copri_clade_A|t__SGB1626 dysbiosis_UC
#>                                                                                                                                         coef_LR
#> UNCLASSIFIED                                                                                                                        -0.03729951
#> k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Phocaeicola|s__Phocaeicola_vulgatus|t__SGB1814    -2.47057587
#> k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_uniformis|t__SGB1836   -4.05438628
#> k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Phocaeicola|s__Phocaeicola_dorei|t__SGB1815       -1.29992480
#> k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_stercoris|t__SGB1830   -3.80748387
#> k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Prevotellaceae|g__Prevotella|s__Prevotella_copri_clade_A|t__SGB1626 -1.14473113
#>                                                                                                                                          pval_LR
#> UNCLASSIFIED                                                                                                                        0.9906647216
#> k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Phocaeicola|s__Phocaeicola_vulgatus|t__SGB1814    0.0010996566
#> k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_uniformis|t__SGB1836   0.0008180096
#> k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Phocaeicola|s__Phocaeicola_dorei|t__SGB1815       0.0693521116
#> k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_stercoris|t__SGB1830   0.0297879923
#> k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Prevotellaceae|g__Prevotella|s__Prevotella_copri_clade_A|t__SGB1626 0.3859268822
```

The same example is also available as an executable script:

``` r

source(system.file(
  "extdata",
  "hmp2_adult_example.R",
  package = "DAssemble",
  mustWork = TRUE
))
```

## Microbiome Example

The `GlobalPatterns` dataset from `phyloseq` contains 16S rRNA profiles
from environmental and host-associated microbiome samples. The example
below compares fecal and soil samples using enhancer-only DAssemble with
CLR normalization.

``` r

data("GlobalPatterns", package = "phyloseq")
gp <- GlobalPatterns

otu <- as(phyloseq::otu_table(gp), "matrix")
metadata <- as.data.frame(phyloseq::sample_data(gp))

keep <- metadata$SampleType %in% c("Feces", "Soil")
otu <- otu[, keep, drop = FALSE]
keep_taxa <- order(rowSums(otu), decreasing = TRUE)[seq_len(250L)]
features <- as.data.frame(t(otu[keep_taxa, , drop = FALSE]))
metadata <- metadata[keep, , drop = FALSE]
metadata$group <- droplevels(factor(metadata$SampleType))
stopifnot(nlevels(metadata$group) == 2L)

se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = t(as.matrix(features)))
)

mae <- MultiAssayExperiment::MultiAssayExperiment(
  experiments = list(microbiome = se),
  colData = S4Vectors::DataFrame(metadata)
)

res <- DAssemble::DAssemble(
  features = mae,
  assay_name = "microbiome",
  core_method = NULL,
  enhancers = c("WLX", "LR", "KS"),
  expVar = "group",
  enhancer_norm = "clr",
  return_components = TRUE,
  return_subensembles = TRUE
)

head(res$res)
#>    feature metadata   pval_WLX  coef_LR     pval_LR    pval_KS pval_joint
#> 28  256977    group 0.05714286 5.615105 0.009103585 0.05714286 0.02074046
#> 37  309788    group 0.05714286 5.615105 0.009103585 0.05714286 0.02074046
#> 39  240687    group 0.05714286 5.615105 0.009103585 0.05714286 0.02074046
#> 40    3046    group 0.05714286 5.615105 0.009103585 0.05714286 0.02074046
#> 44  573135    group 0.05714286 5.615105 0.009103585 0.05714286 0.02074046
#> 50  220494    group 0.05714286 5.615105 0.009103585 0.05714286 0.02074046
#>         qval
#> 28 0.4518957
#> 37 0.4518957
#> 39 0.4518957
#> 40 0.4518957
#> 44 0.4518957
#> 50 0.4518957
names(res$components)
#> [1] "WLX" "LR"  "KS"
head(res$ensembles)
#> $WLX
#>     feature       Pval   adjPval
#> 1    256977 0.05714286 0.4230704
#> 2    309788 0.05714286 0.4230704
#> 3    240687 0.05714286 0.4230704
#> 4      3046 0.05714286 0.4230704
#> 5    573135 0.05714286 0.4230704
#> 6    220494 0.05714286 0.4230704
#> 7    306180 0.05714286 0.4230704
#> 8     27669 0.05714286 0.4230704
#> 9    565399 0.05714286 0.4230704
#> 10   111806 0.05714286 0.4230704
#> 11    65683 0.05714286 0.4230704
#> 12    71074 0.05714286 0.4230704
#> 13   509622 0.05714286 0.4230704
#> 14   360298 0.05714286 0.4230704
#> 15   223329 0.05714286 0.4230704
#> 16   313985 0.05714286 0.4230704
#> 17   574458 0.05714286 0.4230704
#> 18   247010 0.05714286 0.4230704
#> 19   154102 0.05714286 0.4230704
#> 20   129759 0.05714286 0.4230704
#> 21   113626 0.05714286 0.4230704
#> 22   273055 0.05714286 0.4230704
#> 23   551424 0.05714286 0.4230704
#> 24   262714 0.05714286 0.4230704
#> 25   340507 0.05714286 0.4230704
#> 26    84768 0.05714286 0.4230704
#> 27   158370 0.05714286 0.4230704
#> 28   566886 0.05714286 0.4230704
#> 29   216191 0.05714286 0.4230704
#> 30   352632 0.05714286 0.4230704
#> 31   245203 0.05714286 0.4230704
#> 32   265043 0.05714286 0.4230704
#> 33   242675 0.05714286 0.4230704
#> 34   559200 0.05714286 0.4230704
#> 35   367659 0.05714286 0.4230704
#> 36   161227 0.05714286 0.4230704
#> 37   277931 0.05714286 0.4230704
#> 38   527683 0.05714286 0.4230704
#> 39   113866 0.05714286 0.4230704
#> 40   356083 0.05714286 0.4230704
#> 41   200741 0.05714286 0.4230704
#> 42   308892 0.05714286 0.4230704
#> 43     1605 0.05714286 0.4230704
#> 44   250291 0.05714286 0.4230704
#> 45   218467 0.05714286 0.4230704
#> 46   286104 0.05714286 0.4230704
#> 47   511844 0.05714286 0.4230704
#> 48   543366 0.05714286 0.4230704
#> 49   535437 0.05714286 0.4230704
#> 50   159711 0.05714286 0.4230704
#> 51   144381 0.05714286 0.4230704
#> 52   584276 0.05714286 0.4230704
#> 53   570888 0.05714286 0.4230704
#> 54   278390 0.05714286 0.4230704
#> 55   170339 0.05714286 0.4230704
#> 56    93865 0.05714286 0.4230704
#> 57   251499 0.05714286 0.4230704
#> 58    70100 0.05714286 0.4230704
#> 59    97561 0.05714286 0.4230704
#> 60   278272 0.05714286 0.4230704
#> 61   130368 0.05714286 0.4230704
#> 62    76270 0.05714286 0.4230704
#> 63   190723 0.05714286 0.4230704
#> 64   208694 0.05714286 0.4230704
#> 65   222729 0.05714286 0.4230704
#> 66   203846 0.05714286 0.4230704
#> 67   156101 0.05714286 0.4230704
#> 68    13034 0.05714286 0.4230704
#> 69   109907 0.05714286 0.4230704
#> 70   326414 0.05714286 0.4230704
#> 71    36155 0.05714286 0.4230704
#> 72    47916 0.05714286 0.4230704
#> 73   279297 0.05714286 0.4230704
#> 74   150925 0.05714286 0.4230704
#> 75   578268 0.05714286 0.4230704
#> 76   208283 0.05714286 0.4230704
#> 77   546313 0.05714286 0.4230704
#> 78   154822 0.05714286 0.4230704
#> 79   565556 0.05714286 0.4230704
#> 80   100954 0.05714286 0.4230704
#> 81   166835 0.05714286 0.4230704
#> 82   150508 0.05714286 0.4230704
#> 83   113959 0.05714286 0.4230704
#> 84   544749 0.05714286 0.4230704
#> 85   306921 0.05714286 0.4230704
#> 86   114339 0.05714286 0.4230704
#> 87   203722 0.05714286 0.4230704
#> 88   184899 0.05714286 0.4230704
#> 89   235270 0.05714286 0.4230704
#> 90   151018 0.05714286 0.4230704
#> 91   252410 0.05714286 0.4230704
#> 92   331647 0.05714286 0.4230704
#> 93   244491 0.05714286 0.4230704
#> 94   211351 0.05714286 0.4230704
#> 95   154230 0.05714286 0.4230704
#> 96   549552 0.05714286 0.4230704
#> 97   265094 0.05714286 0.4230704
#> 98   215700 0.05714286 0.4230704
#> 99   582012 0.05714286 0.4230704
#> 100  204592 0.05714286 0.4230704
#> 101  252778 0.05714286 0.4230704
#> 102    1791 0.05714286 0.4230704
#> 103  206632 0.05714286 0.4230704
#> 104  575937 0.05714286 0.4230704
#> 105  349634 0.05714286 0.4230704
#> 106  553113 0.05714286 0.4230704
#> 107  210176 0.05714286 0.4230704
#> 108  308535 0.05714286 0.4230704
#> 109  113921 0.05714286 0.4230704
#> 110  534348 0.05714286 0.4230704
#> 111  189047 0.05714286 0.4230704
#> 112  192573 0.05714286 0.4230704
#> 113  180658 0.05714286 0.4230704
#> 114  325407 0.05714286 0.4230704
#> 115  361496 0.05714286 0.4230704
#> 116  316732 0.05714286 0.4230704
#> 117  291090 0.05714286 0.4230704
#> 118  196731 0.05714286 0.4230704
#> 119  290260 0.05714286 0.4230704
#> 120  361480 0.05714286 0.4230704
#> 121  193761 0.05714286 0.4230704
#> 122  357795 0.05714286 0.4230704
#> 123  223059 0.05714286 0.4230704
#> 124  287978 0.05714286 0.4230704
#> 125  357470 0.05714286 0.4230704
#> 126  367433 0.05714286 0.4230704
#> 127  187524 0.05714286 0.4230704
#> 128  518438 0.05714286 0.4230704
#> 129  327476 0.05714286 0.4230704
#> 130  199487 0.05714286 0.4230704
#> 131  191541 0.05714286 0.4230704
#> 132  365676 0.05714286 0.4230704
#> 133  269833 0.05714286 0.4230704
#> 134  194053 0.05714286 0.4230704
#> 135  190464 0.05714286 0.4230704
#> 136  288134 0.05714286 0.4230704
#> 137  194648 0.05714286 0.4230704
#> 138  138006 0.05714286 0.4230704
#> 139  212619 0.05714286 0.4230704
#> 140  293896 0.05714286 0.4230704
#> 141  195173 0.05714286 0.4230704
#> 142  192127 0.05714286 0.4230704
#> 143  193066 0.05714286 0.4230704
#> 144  177339 0.05714286 0.4230704
#> 145  193343 0.05714286 0.4230704
#> 146  203708 0.05714286 0.4230704
#> 147  261912 0.05714286 0.4230704
#> 148  175537 0.05714286 0.4230704
#> 149  557785 0.11428571 0.7995643
#> 150  179460 0.05714286 0.4230704
#> 151  181843 0.05714286 0.4230704
#> 152  310301 0.05714286 0.4230704
#> 153  162162 0.05714286 0.4230704
#> 154  190108 0.05714286 0.4230704
#> 155   82283 0.05714286 0.4230704
#> 156  262724 0.05714286 0.4230704
#> 157  111986 0.05714286 0.4230704
#> 158  250258 0.11428571 0.7995643
#> 159  188646 0.05714286 0.4230704
#> 160  469991 0.05714286 0.4230704
#> 161  346780 0.05714286 0.4230704
#> 162  174348 0.05714286 0.4230704
#> 163  185281 0.05714286 0.4230704
#> 164  167602 0.05714286 0.4230704
#> 165  145786 0.62857143 1.0000000
#> 166  182980 0.05714286 0.4230704
#> 167  204044 0.05714286 0.4230704
#> 168   51024 0.11428571 0.7995643
#> 169  255584 0.22857143 1.0000000
#> 170  216862 0.05714286 0.4230704
#> 171  154268 0.11428571 0.7995643
#> 172  547632 0.11428571 0.7995643
#> 173  294022 0.05714286 0.4230704
#> 174  259569 0.05714286 0.4230704
#> 175  348374 0.05714286 0.4230704
#> 176  470973 0.05714286 0.4230704
#> 177  589024 0.05714286 0.4230704
#> 178  268540 0.05714286 0.4230704
#> 179  302160 0.05714286 0.4230704
#> 180  352304 0.05714286 0.4230704
#> 181    2000 0.05714286 0.4230704
#> 182  109014 0.05714286 0.4230704
#> 183  366612 0.05714286 0.4230704
#> 184  170206 0.05714286 0.4230704
#> 185  104780 0.05714286 0.4230704
#> 186  511008 0.05714286 0.4230704
#> 187  278234 0.05714286 0.4230704
#> 188  213747 0.05714286 0.4230704
#> 189  190872 0.05714286 0.4230704
#> 190  178915 0.05714286 0.4230704
#> 191  191306 0.05714286 0.4230704
#> 192  510994 0.05714286 0.4230704
#> 193   21698 0.05714286 0.4230704
#> 194  193632 0.05714286 0.4230704
#> 195  469873 0.05714286 0.4230704
#> 196  328422 0.05714286 0.4230704
#> 197  405524 0.05714286 0.4230704
#> 198  367055 0.05714286 0.4230704
#> 199  470114 0.05714286 0.4230704
#> 200  544358 0.05714286 0.4230704
#> 201  204889 0.05714286 0.4230704
#> 202  470812 0.05714286 0.4230704
#> 203  363692 0.05714286 0.4230704
#> 204  171551 0.05714286 0.4230704
#> 205  322235 0.05714286 0.4230704
#> 206  331820 0.05714286 0.4230704
#> 207  158660 0.05714286 0.4230704
#> 208  244304 0.62857143 1.0000000
#> 209  263681 0.85714286 1.0000000
#> 210  326977 0.05714286 0.4230704
#> 211  278873 0.85714286 1.0000000
#> 212  357591 0.40000000 1.0000000
#> 213  127309 1.00000000 1.0000000
#> 214  517293 0.85714286 1.0000000
#> 215  248140 0.05714286 0.4230704
#> 216  211720 1.00000000 1.0000000
#> 217  470172 0.05714286 0.4230704
#> 218  549871 0.62857143 1.0000000
#> 219  347862 0.11428571 0.7995643
#> 220  249661 0.40000000 1.0000000
#> 221  183200 0.05714286 0.4230704
#> 222   42680 0.40000000 1.0000000
#> 223  326792 0.11428571 0.7995643
#> 224  171772 0.22857143 1.0000000
#> 225  541301 0.85714286 1.0000000
#> 226  469482 0.11428571 0.7995643
#> 227  181095 0.22857143 1.0000000
#> 228  299892 0.40000000 1.0000000
#> 229  175497 0.22857143 1.0000000
#> 230  366882 1.00000000 1.0000000
#> 231  109546 0.22857143 1.0000000
#> 232  320147 0.22857143 1.0000000
#> 233  243150 0.22857143 1.0000000
#> 234  191687 0.22857143 1.0000000
#> 235  302974 0.22857143 1.0000000
#> 236  181892 0.11428571 0.7995643
#> 237  181028 0.11428571 0.7995643
#> 238  550814 0.11428571 0.7995643
#> 239  360268 0.22857143 1.0000000
#> 240  263461 0.22857143 1.0000000
#> 241  275891 0.11428571 0.7995643
#> 242  368928 0.22857143 1.0000000
#> 243  230437 0.40000000 1.0000000
#> 244  212128 0.05714286 0.4230704
#> 245  469946 0.85714286 1.0000000
#> 246  103685 1.00000000 1.0000000
#> 247  565565 0.40000000 1.0000000
#> 248  197817 0.62857143 1.0000000
#> 249  173810 1.00000000 1.0000000
#> 250  215890 0.05714286 0.4230704
#> 
#> $LR
#>     feature        Pval   adjPval
#> 1    256977 0.009103585 0.1946615
#> 2    309788 0.009103585 0.1946615
#> 3    240687 0.009103585 0.1946615
#> 4      3046 0.009103585 0.1946615
#> 5    573135 0.009103585 0.1946615
#> 6    220494 0.009103585 0.1946615
#> 7    306180 0.009103585 0.1946615
#> 8     27669 0.009103585 0.1946615
#> 9    565399 0.009103585 0.1946615
#> 10   111806 0.009103585 0.1946615
#> 11    65683 0.009103585 0.1946615
#> 12    71074 0.009103585 0.1946615
#> 13   509622 0.009103585 0.1946615
#> 14   360298 0.009103585 0.1946615
#> 15   223329 0.009103585 0.1946615
#> 16   313985 0.009103585 0.1946615
#> 17   574458 0.009103585 0.1946615
#> 18   247010 0.009103585 0.1946615
#> 19   154102 0.009103585 0.1946615
#> 20   129759 0.009103585 0.1946615
#> 21   113626 0.009103585 0.1946615
#> 22   273055 0.009103585 0.1946615
#> 23   551424 0.009103585 0.1946615
#> 24   262714 0.009103585 0.1946615
#> 25   340507 0.009103585 0.1946615
#> 26    84768 0.009103585 0.1946615
#> 27   158370 0.009103585 0.1946615
#> 28   566886 0.009103585 0.1946615
#> 29   216191 0.009103585 0.1946615
#> 30   352632 0.009103585 0.1946615
#> 31   245203 0.009103585 0.1946615
#> 32   265043 0.009103585 0.1946615
#> 33   242675 0.009103585 0.1946615
#> 34   559200 0.009103585 0.1946615
#> 35   367659 0.009103585 0.1946615
#> 36   161227 0.009103585 0.1946615
#> 37   277931 0.009103585 0.1946615
#> 38   527683 0.009103585 0.1946615
#> 39   113866 0.009103585 0.1946615
#> 40   356083 0.009103585 0.1946615
#> 41   200741 0.009103585 0.1946615
#> 42   308892 0.009103585 0.1946615
#> 43     1605 0.009103585 0.1946615
#> 44   250291 0.009103585 0.1946615
#> 45   218467 0.009103585 0.1946615
#> 46   286104 0.009103585 0.1946615
#> 47   511844 0.009103585 0.1946615
#> 48   543366 0.009103585 0.1946615
#> 49   535437 0.009103585 0.1946615
#> 50   159711 0.009103585 0.1946615
#> 51   144381 0.009103585 0.1946615
#> 52   584276 0.009103585 0.1946615
#> 53   570888 0.009103585 0.1946615
#> 54   278390 0.009103585 0.1946615
#> 55   170339 0.009103585 0.1946615
#> 56    93865 0.009103585 0.1946615
#> 57   251499 0.009103585 0.1946615
#> 58    70100 0.009103585 0.1946615
#> 59    97561 0.009103585 0.1946615
#> 60   278272 0.009103585 0.1946615
#> 61   130368 0.009103585 0.1946615
#> 62    76270 0.009103585 0.1946615
#> 63   190723 0.009103585 0.1946615
#> 64   208694 0.009103585 0.1946615
#> 65   222729 0.009103585 0.1946615
#> 66   203846 0.009103585 0.1946615
#> 67   156101 0.009103585 0.1946615
#> 68    13034 0.009103585 0.1946615
#> 69   109907 0.009103585 0.1946615
#> 70   326414 0.009103585 0.1946615
#> 71    36155 0.022058541 0.2821997
#> 72    47916 0.022058541 0.2821997
#> 73   279297 0.022058541 0.2821997
#> 74   150925 0.022058541 0.2821997
#> 75   578268 0.022058541 0.2821997
#> 76   208283 0.022058541 0.2821997
#> 77   546313 0.022058541 0.2821997
#> 78   154822 0.022058541 0.2821997
#> 79   565556 0.022058541 0.2821997
#> 80   100954 0.022058541 0.2821997
#> 81   166835 0.022058541 0.2821997
#> 82   150508 0.022058541 0.2821997
#> 83   113959 0.022058541 0.2821997
#> 84   544749 0.022058541 0.2821997
#> 85   306921 0.022058541 0.2821997
#> 86   114339 0.022058541 0.2821997
#> 87   203722 0.022058541 0.2821997
#> 88   184899 0.022058541 0.2821997
#> 89   235270 0.022058541 0.2821997
#> 90   151018 0.022058541 0.2821997
#> 91   252410 0.022058541 0.2821997
#> 92   331647 0.022058541 0.2821997
#> 93   244491 0.022058541 0.2821997
#> 94   211351 0.022058541 0.2821997
#> 95   154230 0.022058541 0.2821997
#> 96   549552 0.022058541 0.2821997
#> 97   265094 0.022058541 0.2821997
#> 98   215700 0.022058541 0.2821997
#> 99   582012 0.022058541 0.2821997
#> 100  204592 0.022058541 0.2821997
#> 101  252778 0.022058541 0.2821997
#> 102    1791 0.022058541 0.2821997
#> 103  206632 0.022058541 0.2821997
#> 104  575937 0.022058541 0.2821997
#> 105  349634 0.022058541 0.2821997
#> 106  553113 0.022058541 0.2821997
#> 107  210176 0.022058541 0.2821997
#> 108  308535 0.022058541 0.2821997
#> 109  113921 0.022058541 0.2821997
#> 110  534348 0.022058541 0.2821997
#> 111  189047 0.187112927 1.0000000
#> 112  192573 0.187112927 1.0000000
#> 113  180658 0.187112927 1.0000000
#> 114  325407 0.187112927 1.0000000
#> 115  361496 0.187112927 1.0000000
#> 116  316732 0.187112927 1.0000000
#> 117  291090 0.187112927 1.0000000
#> 118  196731 0.187112927 1.0000000
#> 119  290260 0.187112927 1.0000000
#> 120  361480 0.187112927 1.0000000
#> 121  193761 0.187112927 1.0000000
#> 122  357795 0.187112927 1.0000000
#> 123  223059 0.187112927 1.0000000
#> 124  287978 0.187112927 1.0000000
#> 125  357470 0.187112927 1.0000000
#> 126  367433 0.187112927 1.0000000
#> 127  187524 0.187112927 1.0000000
#> 128  518438 0.187112927 1.0000000
#> 129  327476 0.187112927 1.0000000
#> 130  199487 0.187112927 1.0000000
#> 131  191541 0.187112927 1.0000000
#> 132  365676 0.187112927 1.0000000
#> 133  269833 0.187112927 1.0000000
#> 134  194053 0.187112927 1.0000000
#> 135  190464 0.187112927 1.0000000
#> 136  288134 0.187112927 1.0000000
#> 137  194648 0.187112927 1.0000000
#> 138  138006 0.187112927 1.0000000
#> 139  212619 0.187112927 1.0000000
#> 140  293896 0.187112927 1.0000000
#> 141  195173 0.187112927 1.0000000
#> 142  192127 0.187112927 1.0000000
#> 143  193066 0.187112927 1.0000000
#> 144  177339 0.187112927 1.0000000
#> 145  193343 0.187112927 1.0000000
#> 146  203708 0.187112927 1.0000000
#> 147  261912 0.187112927 1.0000000
#> 148  175537 0.187112927 1.0000000
#> 149  557785 0.022058541 0.2821997
#> 150  179460 0.187112927 1.0000000
#> 151  181843 0.187112927 1.0000000
#> 152  310301 0.187112927 1.0000000
#> 153  162162 0.187112927 1.0000000
#> 154  190108 0.187112927 1.0000000
#> 155   82283 0.187112927 1.0000000
#> 156  262724 0.187112927 1.0000000
#> 157  111986 0.093494567 1.0000000
#> 158  250258 0.022058541 0.2821997
#> 159  188646 0.187112927 1.0000000
#> 160  469991 0.187112927 1.0000000
#> 161  346780 0.187112927 1.0000000
#> 162  174348 0.187112927 1.0000000
#> 163  185281 0.187112927 1.0000000
#> 164  167602 0.187112927 1.0000000
#> 165  145786 0.022058541 0.2821997
#> 166  182980 0.187112927 1.0000000
#> 167  204044 0.187112927 1.0000000
#> 168   51024 0.022058541 0.2821997
#> 169  255584 0.022058541 0.2821997
#> 170  216862 0.187112927 1.0000000
#> 171  154268 0.022058541 0.2821997
#> 172  547632 0.022058541 0.2821997
#> 173  294022 0.187112927 1.0000000
#> 174  259569 0.471878302 1.0000000
#> 175  348374 0.408906218 1.0000000
#> 176  470973 0.471878302 1.0000000
#> 177  589024 0.408906218 1.0000000
#> 178  268540 0.408906218 1.0000000
#> 179  302160 0.408906218 1.0000000
#> 180  352304 0.471878302 1.0000000
#> 181    2000 0.408906218 1.0000000
#> 182  109014 0.408906218 1.0000000
#> 183  366612 0.408906218 1.0000000
#> 184  170206 0.408906218 1.0000000
#> 185  104780 0.408906218 1.0000000
#> 186  511008 0.471878302 1.0000000
#> 187  278234 0.408906218 1.0000000
#> 188  213747 0.408906218 1.0000000
#> 189  190872 0.408906218 1.0000000
#> 190  178915 0.408906218 1.0000000
#> 191  191306 0.408906218 1.0000000
#> 192  510994 0.471878302 1.0000000
#> 193   21698 0.408906218 1.0000000
#> 194  193632 0.408906218 1.0000000
#> 195  469873 0.696765470 1.0000000
#> 196  328422 0.696765470 1.0000000
#> 197  405524 0.696765470 1.0000000
#> 198  367055 0.696765470 1.0000000
#> 199  470114 0.696765470 1.0000000
#> 200  544358 0.696765470 1.0000000
#> 201  204889 0.696765470 1.0000000
#> 202  470812 0.696765470 1.0000000
#> 203  363692 0.696765470 1.0000000
#> 204  171551 0.915051949 1.0000000
#> 205  322235 0.915051949 1.0000000
#> 206  331820          NA        NA
#> 207  158660          NA        NA
#> 208  244304 0.165274665 1.0000000
#> 209  263681 0.241191857 1.0000000
#> 210  326977          NA        NA
#> 211  278873 0.950713866 1.0000000
#> 212  357591 0.408906218 1.0000000
#> 213  127309 0.241191857 1.0000000
#> 214  517293 0.950713866 1.0000000
#> 215  248140          NA        NA
#> 216  211720 0.950713866 1.0000000
#> 217  470172 0.950713866 1.0000000
#> 218  549871 0.696765470 1.0000000
#> 219  347862 0.408906218 1.0000000
#> 220  249661 0.950713866 1.0000000
#> 221  183200 0.950713866 1.0000000
#> 222   42680 0.696765470 1.0000000
#> 223  326792 0.408906218 1.0000000
#> 224  171772 0.696765470 1.0000000
#> 225  541301 0.696765470 1.0000000
#> 226  469482 0.696765470 1.0000000
#> 227  181095 0.696765470 1.0000000
#> 228  299892 0.696765470 1.0000000
#> 229  175497 0.696765470 1.0000000
#> 230  366882 0.950713866 1.0000000
#> 231  109546 0.696765470 1.0000000
#> 232  320147 0.950713866 1.0000000
#> 233  243150 0.950713866 1.0000000
#> 234  191687 0.408906218 1.0000000
#> 235  302974 0.696765470 1.0000000
#> 236  181892 0.408906218 1.0000000
#> 237  181028 0.950713866 1.0000000
#> 238  550814 0.950713866 1.0000000
#> 239  360268 0.950713866 1.0000000
#> 240  263461 0.408906218 1.0000000
#> 241  275891 0.950713866 1.0000000
#> 242  368928 0.950713866 1.0000000
#> 243  230437 0.950713866 1.0000000
#> 244  212128 0.950713866 1.0000000
#> 245  469946 0.950713866 1.0000000
#> 246  103685 0.950713866 1.0000000
#> 247  565565 0.950713866 1.0000000
#> 248  197817 0.950713866 1.0000000
#> 249  173810 0.950713866 1.0000000
#> 250  215890 0.950713866 1.0000000
#> 
#> $KS
#>     feature       Pval   adjPval
#> 1    256977 0.05714286 0.4230704
#> 2    309788 0.05714286 0.4230704
#> 3    240687 0.05714286 0.4230704
#> 4      3046 0.05714286 0.4230704
#> 5    573135 0.05714286 0.4230704
#> 6    220494 0.05714286 0.4230704
#> 7    306180 0.05714286 0.4230704
#> 8     27669 0.05714286 0.4230704
#> 9    565399 0.05714286 0.4230704
#> 10   111806 0.05714286 0.4230704
#> 11    65683 0.05714286 0.4230704
#> 12    71074 0.05714286 0.4230704
#> 13   509622 0.05714286 0.4230704
#> 14   360298 0.05714286 0.4230704
#> 15   223329 0.05714286 0.4230704
#> 16   313985 0.05714286 0.4230704
#> 17   574458 0.05714286 0.4230704
#> 18   247010 0.05714286 0.4230704
#> 19   154102 0.05714286 0.4230704
#> 20   129759 0.05714286 0.4230704
#> 21   113626 0.05714286 0.4230704
#> 22   273055 0.05714286 0.4230704
#> 23   551424 0.05714286 0.4230704
#> 24   262714 0.05714286 0.4230704
#> 25   340507 0.05714286 0.4230704
#> 26    84768 0.05714286 0.4230704
#> 27   158370 0.05714286 0.4230704
#> 28   566886 0.05714286 0.4230704
#> 29   216191 0.05714286 0.4230704
#> 30   352632 0.05714286 0.4230704
#> 31   245203 0.05714286 0.4230704
#> 32   265043 0.05714286 0.4230704
#> 33   242675 0.05714286 0.4230704
#> 34   559200 0.05714286 0.4230704
#> 35   367659 0.05714286 0.4230704
#> 36   161227 0.05714286 0.4230704
#> 37   277931 0.05714286 0.4230704
#> 38   527683 0.05714286 0.4230704
#> 39   113866 0.05714286 0.4230704
#> 40   356083 0.05714286 0.4230704
#> 41   200741 0.05714286 0.4230704
#> 42   308892 0.05714286 0.4230704
#> 43     1605 0.05714286 0.4230704
#> 44   250291 0.05714286 0.4230704
#> 45   218467 0.05714286 0.4230704
#> 46   286104 0.05714286 0.4230704
#> 47   511844 0.05714286 0.4230704
#> 48   543366 0.05714286 0.4230704
#> 49   535437 0.05714286 0.4230704
#> 50   159711 0.05714286 0.4230704
#> 51   144381 0.05714286 0.4230704
#> 52   584276 0.05714286 0.4230704
#> 53   570888 0.05714286 0.4230704
#> 54   278390 0.05714286 0.4230704
#> 55   170339 0.05714286 0.4230704
#> 56    93865 0.05714286 0.4230704
#> 57   251499 0.05714286 0.4230704
#> 58    70100 0.05714286 0.4230704
#> 59    97561 0.05714286 0.4230704
#> 60   278272 0.05714286 0.4230704
#> 61   130368 0.05714286 0.4230704
#> 62    76270 0.05714286 0.4230704
#> 63   190723 0.05714286 0.4230704
#> 64   208694 0.05714286 0.4230704
#> 65   222729 0.05714286 0.4230704
#> 66   203846 0.05714286 0.4230704
#> 67   156101 0.05714286 0.4230704
#> 68    13034 0.05714286 0.4230704
#> 69   109907 0.05714286 0.4230704
#> 70   326414 0.05714286 0.4230704
#> 71    36155 0.05714286 0.4230704
#> 72    47916 0.05714286 0.4230704
#> 73   279297 0.05714286 0.4230704
#> 74   150925 0.05714286 0.4230704
#> 75   578268 0.05714286 0.4230704
#> 76   208283 0.05714286 0.4230704
#> 77   546313 0.05714286 0.4230704
#> 78   154822 0.05714286 0.4230704
#> 79   565556 0.05714286 0.4230704
#> 80   100954 0.05714286 0.4230704
#> 81   166835 0.05714286 0.4230704
#> 82   150508 0.05714286 0.4230704
#> 83   113959 0.05714286 0.4230704
#> 84   544749 0.05714286 0.4230704
#> 85   306921 0.05714286 0.4230704
#> 86   114339 0.05714286 0.4230704
#> 87   203722 0.05714286 0.4230704
#> 88   184899 0.05714286 0.4230704
#> 89   235270 0.05714286 0.4230704
#> 90   151018 0.05714286 0.4230704
#> 91   252410 0.05714286 0.4230704
#> 92   331647 0.05714286 0.4230704
#> 93   244491 0.05714286 0.4230704
#> 94   211351 0.05714286 0.4230704
#> 95   154230 0.05714286 0.4230704
#> 96   549552 0.05714286 0.4230704
#> 97   265094 0.05714286 0.4230704
#> 98   215700 0.05714286 0.4230704
#> 99   582012 0.05714286 0.4230704
#> 100  204592 0.05714286 0.4230704
#> 101  252778 0.05714286 0.4230704
#> 102    1791 0.05714286 0.4230704
#> 103  206632 0.05714286 0.4230704
#> 104  575937 0.05714286 0.4230704
#> 105  349634 0.05714286 0.4230704
#> 106  553113 0.05714286 0.4230704
#> 107  210176 0.05714286 0.4230704
#> 108  308535 0.05714286 0.4230704
#> 109  113921 0.05714286 0.4230704
#> 110  534348 0.05714286 0.4230704
#> 111  189047 0.05714286 0.4230704
#> 112  192573 0.05714286 0.4230704
#> 113  180658 0.05714286 0.4230704
#> 114  325407 0.05714286 0.4230704
#> 115  361496 0.05714286 0.4230704
#> 116  316732 0.05714286 0.4230704
#> 117  291090 0.05714286 0.4230704
#> 118  196731 0.05714286 0.4230704
#> 119  290260 0.05714286 0.4230704
#> 120  361480 0.05714286 0.4230704
#> 121  193761 0.05714286 0.4230704
#> 122  357795 0.05714286 0.4230704
#> 123  223059 0.05714286 0.4230704
#> 124  287978 0.05714286 0.4230704
#> 125  357470 0.05714286 0.4230704
#> 126  367433 0.05714286 0.4230704
#> 127  187524 0.05714286 0.4230704
#> 128  518438 0.05714286 0.4230704
#> 129  327476 0.05714286 0.4230704
#> 130  199487 0.05714286 0.4230704
#> 131  191541 0.05714286 0.4230704
#> 132  365676 0.05714286 0.4230704
#> 133  269833 0.05714286 0.4230704
#> 134  194053 0.05714286 0.4230704
#> 135  190464 0.05714286 0.4230704
#> 136  288134 0.05714286 0.4230704
#> 137  194648 0.05714286 0.4230704
#> 138  138006 0.05714286 0.4230704
#> 139  212619 0.05714286 0.4230704
#> 140  293896 0.05714286 0.4230704
#> 141  195173 0.05714286 0.4230704
#> 142  192127 0.05714286 0.4230704
#> 143  193066 0.05714286 0.4230704
#> 144  177339 0.05714286 0.4230704
#> 145  193343 0.05714286 0.4230704
#> 146  203708 0.05714286 0.4230704
#> 147  261912 0.05714286 0.4230704
#> 148  175537 0.05714286 0.4230704
#> 149  557785 0.22857143 1.0000000
#> 150  179460 0.05714286 0.4230704
#> 151  181843 0.05714286 0.4230704
#> 152  310301 0.05714286 0.4230704
#> 153  162162 0.05714286 0.4230704
#> 154  190108 0.05714286 0.4230704
#> 155   82283 0.05714286 0.4230704
#> 156  262724 0.05714286 0.4230704
#> 157  111986 0.05714286 0.4230704
#> 158  250258 0.22857143 1.0000000
#> 159  188646 0.05714286 0.4230704
#> 160  469991 0.05714286 0.4230704
#> 161  346780 0.05714286 0.4230704
#> 162  174348 0.05714286 0.4230704
#> 163  185281 0.05714286 0.4230704
#> 164  167602 0.05714286 0.4230704
#> 165  145786 0.40000000 1.0000000
#> 166  182980 0.05714286 0.4230704
#> 167  204044 0.05714286 0.4230704
#> 168   51024 0.22857143 1.0000000
#> 169  255584 0.40000000 1.0000000
#> 170  216862 0.05714286 0.4230704
#> 171  154268 0.22857143 1.0000000
#> 172  547632 0.22857143 1.0000000
#> 173  294022 0.05714286 0.4230704
#> 174  259569 0.05714286 0.4230704
#> 175  348374 0.05714286 0.4230704
#> 176  470973 0.05714286 0.4230704
#> 177  589024 0.05714286 0.4230704
#> 178  268540 0.05714286 0.4230704
#> 179  302160 0.05714286 0.4230704
#> 180  352304 0.05714286 0.4230704
#> 181    2000 0.05714286 0.4230704
#> 182  109014 0.05714286 0.4230704
#> 183  366612 0.05714286 0.4230704
#> 184  170206 0.05714286 0.4230704
#> 185  104780 0.05714286 0.4230704
#> 186  511008 0.05714286 0.4230704
#> 187  278234 0.05714286 0.4230704
#> 188  213747 0.05714286 0.4230704
#> 189  190872 0.05714286 0.4230704
#> 190  178915 0.05714286 0.4230704
#> 191  191306 0.05714286 0.4230704
#> 192  510994 0.05714286 0.4230704
#> 193   21698 0.05714286 0.4230704
#> 194  193632 0.05714286 0.4230704
#> 195  469873 0.05714286 0.4230704
#> 196  328422 0.05714286 0.4230704
#> 197  405524 0.05714286 0.4230704
#> 198  367055 0.05714286 0.4230704
#> 199  470114 0.05714286 0.4230704
#> 200  544358 0.05714286 0.4230704
#> 201  204889 0.05714286 0.4230704
#> 202  470812 0.05714286 0.4230704
#> 203  363692 0.05714286 0.4230704
#> 204  171551 0.05714286 0.4230704
#> 205  322235 0.05714286 0.4230704
#> 206  331820 0.05714286 0.4230704
#> 207  158660 0.05714286 0.4230704
#> 208  244304 0.65714286 1.0000000
#> 209  263681 0.65714286 1.0000000
#> 210  326977 0.05714286 0.4230704
#> 211  278873 0.88571429 1.0000000
#> 212  357591 0.22857143 1.0000000
#> 213  127309 1.00000000 1.0000000
#> 214  517293 0.65714286 1.0000000
#> 215  248140 0.05714286 0.4230704
#> 216  211720 1.00000000 1.0000000
#> 217  470172 0.05714286 0.4230704
#> 218  549871 0.65714286 1.0000000
#> 219  347862 0.22857143 1.0000000
#> 220  249661 0.40000000 1.0000000
#> 221  183200 0.05714286 0.4230704
#> 222   42680 0.65714286 1.0000000
#> 223  326792 0.22857143 1.0000000
#> 224  171772 0.40000000 1.0000000
#> 225  541301 0.65714286 1.0000000
#> 226  469482 0.22857143 1.0000000
#> 227  181095 0.40000000 1.0000000
#> 228  299892 0.22857143 1.0000000
#> 229  175497 0.22857143 1.0000000
#> 230  366882 0.88571429 1.0000000
#> 231  109546 0.22857143 1.0000000
#> 232  320147 0.40000000 1.0000000
#> 233  243150 0.22857143 1.0000000
#> 234  191687 0.22857143 1.0000000
#> 235  302974 0.40000000 1.0000000
#> 236  181892 0.22857143 1.0000000
#> 237  181028 0.22857143 1.0000000
#> 238  550814 0.22857143 1.0000000
#> 239  360268 0.40000000 1.0000000
#> 240  263461 0.22857143 1.0000000
#> 241  275891 0.22857143 1.0000000
#> 242  368928 0.22857143 1.0000000
#> 243  230437 0.40000000 1.0000000
#> 244  212128 0.05714286 0.4230704
#> 245  469946 0.97142857 1.0000000
#> 246  103685 0.88571429 1.0000000
#> 247  565565 0.40000000 1.0000000
#> 248  197817 0.88571429 1.0000000
#> 249  173810 0.88571429 1.0000000
#> 250  215890 0.05714286 0.4230704
#> 
#> $`WLX+LR`
#>     feature       Pval   adjPval
#> 28   256977 0.01571932 0.3424946
#> 37   309788 0.01571932 0.3424946
#> 39   240687 0.01571932 0.3424946
#> 40     3046 0.01571932 0.3424946
#> 44   573135 0.01571932 0.3424946
#> 50   220494 0.01571932 0.3424946
#> 56   306180 0.01571932 0.3424946
#> 61    27669 0.01571932 0.3424946
#> 67   565399 0.01571932 0.3424946
#> 68   111806 0.01571932 0.3424946
#> 75    65683 0.01571932 0.3424946
#> 84    71074 0.01571932 0.3424946
#> 85   509622 0.01571932 0.3424946
#> 90   360298 0.01571932 0.3424946
#> 100  223329 0.01571932 0.3424946
#> 113  313985 0.01571932 0.3424946
#> 114  574458 0.01571932 0.3424946
#> 119  247010 0.01571932 0.3424946
#> 120  154102 0.01571932 0.3424946
#> 121  129759 0.01571932 0.3424946
#> 122  113626 0.01571932 0.3424946
#> 125  273055 0.01571932 0.3424946
#> 126  551424 0.01571932 0.3424946
#> 127  262714 0.01571932 0.3424946
#> 129  340507 0.01571932 0.3424946
#> 130   84768 0.01571932 0.3424946
#> 135  158370 0.01571932 0.3424946
#> 136  566886 0.01571932 0.3424946
#> 142  216191 0.01571932 0.3424946
#> 144  352632 0.01571932 0.3424946
#> 145  245203 0.01571932 0.3424946
#> 149  265043 0.01571932 0.3424946
#> 156  242675 0.01571932 0.3424946
#> 160  559200 0.01571932 0.3424946
#> 161  367659 0.01571932 0.3424946
#> 162  161227 0.01571932 0.3424946
#> 164  277931 0.01571932 0.3424946
#> 167  527683 0.01571932 0.3424946
#> 170  113866 0.01571932 0.3424946
#> 176  356083 0.01571932 0.3424946
#> 178  200741 0.01571932 0.3424946
#> 185  308892 0.01571932 0.3424946
#> 187    1605 0.01571932 0.3424946
#> 188  250291 0.01571932 0.3424946
#> 190  218467 0.01571932 0.3424946
#> 191  286104 0.01571932 0.3424946
#> 192  511844 0.01571932 0.3424946
#> 196  543366 0.01571932 0.3424946
#> 198  535437 0.01571932 0.3424946
#> 200  159711 0.01571932 0.3424946
#> 201  144381 0.01571932 0.3424946
#> 205  584276 0.01571932 0.3424946
#> 207  570888 0.01571932 0.3424946
#> 208  278390 0.01571932 0.3424946
#> 213  170339 0.01571932 0.3424946
#> 215   93865 0.01571932 0.3424946
#> 219  251499 0.01571932 0.3424946
#> 221   70100 0.01571932 0.3424946
#> 222   97561 0.01571932 0.3424946
#> 223  278272 0.01571932 0.3424946
#> 226  130368 0.01571932 0.3424946
#> 227   76270 0.01571932 0.3424946
#> 228  190723 0.01571932 0.3424946
#> 229  208694 0.01571932 0.3424946
#> 231  222729 0.01571932 0.3424946
#> 233  203846 0.01571932 0.3424946
#> 234  156101 0.01571932 0.3424946
#> 240   13034 0.01571932 0.3424946
#> 244  109907 0.01571932 0.3424946
#> 248  326414 0.01571932 0.3424946
#> 34    36155 0.03185593 0.4416879
#> 46    47916 0.03185593 0.4416879
#> 55   279297 0.03185593 0.4416879
#> 64   150925 0.03185593 0.4416879
#> 71   578268 0.03185593 0.4416879
#> 81   208283 0.03185593 0.4416879
#> 89   546313 0.03185593 0.4416879
#> 92   154822 0.03185593 0.4416879
#> 97   565556 0.03185593 0.4416879
#> 99   100954 0.03185593 0.4416879
#> 109  166835 0.03185593 0.4416879
#> 110  150508 0.03185593 0.4416879
#> 115  113959 0.03185593 0.4416879
#> 118  544749 0.03185593 0.4416879
#> 138  306921 0.03185593 0.4416879
#> 139  114339 0.03185593 0.4416879
#> 147  203722 0.03185593 0.4416879
#> 151  184899 0.03185593 0.4416879
#> 165  235270 0.03185593 0.4416879
#> 168  151018 0.03185593 0.4416879
#> 174  252410 0.03185593 0.4416879
#> 175  331647 0.03185593 0.4416879
#> 179  244491 0.03185593 0.4416879
#> 186  211351 0.03185593 0.4416879
#> 195  154230 0.03185593 0.4416879
#> 197  549552 0.03185593 0.4416879
#> 199  265094 0.03185593 0.4416879
#> 203  215700 0.03185593 0.4416879
#> 214  582012 0.03185593 0.4416879
#> 217  204592 0.03185593 0.4416879
#> 218  252778 0.03185593 0.4416879
#> 220    1791 0.03185593 0.4416879
#> 232  206632 0.03185593 0.4416879
#> 235  575937 0.03185593 0.4416879
#> 236  349634 0.03185593 0.4416879
#> 239  553113 0.03185593 0.4416879
#> 241  210176 0.03185593 0.4416879
#> 242  308535 0.03185593 0.4416879
#> 245  113921 0.03185593 0.4416879
#> 250  534348 0.03185593 0.4416879
#> 3    189047 0.08845351 0.7798066
#> 7    192573 0.08845351 0.7798066
#> 9    180658 0.08845351 0.7798066
#> 12   325407 0.08845351 0.7798066
#> 13   361496 0.08845351 0.7798066
#> 14   316732 0.08845351 0.7798066
#> 15   291090 0.08845351 0.7798066
#> 17   196731 0.08845351 0.7798066
#> 19   290260 0.08845351 0.7798066
#> 20   361480 0.08845351 0.7798066
#> 21   193761 0.08845351 0.7798066
#> 24   357795 0.08845351 0.7798066
#> 26   223059 0.08845351 0.7798066
#> 29   287978 0.08845351 0.7798066
#> 30   357470 0.08845351 0.7798066
#> 31   367433 0.08845351 0.7798066
#> 32   187524 0.08845351 0.7798066
#> 35   518438 0.08845351 0.7798066
#> 36   327476 0.08845351 0.7798066
#> 42   199487 0.08845351 0.7798066
#> 43   191541 0.08845351 0.7798066
#> 45   365676 0.08845351 0.7798066
#> 48   269833 0.08845351 0.7798066
#> 49   194053 0.08845351 0.7798066
#> 51   190464 0.08845351 0.7798066
#> 52   288134 0.08845351 0.7798066
#> 53   194648 0.08845351 0.7798066
#> 58   138006 0.08845351 0.7798066
#> 59   212619 0.08845351 0.7798066
#> 62   293896 0.08845351 0.7798066
#> 63   195173 0.08845351 0.7798066
#> 66   192127 0.08845351 0.7798066
#> 70   193066 0.08845351 0.7798066
#> 76   177339 0.08845351 0.7798066
#> 77   193343 0.08845351 0.7798066
#> 78   203708 0.08845351 0.7798066
#> 83   261912 0.08845351 0.7798066
#> 86   175537 0.08845351 0.7798066
#> 87   557785 0.03712191 0.4923233
#> 96   179460 0.08845351 0.7798066
#> 101  181843 0.08845351 0.7798066
#> 103  310301 0.08845351 0.7798066
#> 106  162162 0.08845351 0.7798066
#> 111  190108 0.08845351 0.7798066
#> 123   82283 0.08845351 0.7798066
#> 124  262724 0.08845351 0.7798066
#> 128  111986 0.07100559 0.7798066
#> 131  250258 0.03712191 0.4923233
#> 141  188646 0.08845351 0.7798066
#> 148  469991 0.08845351 0.7798066
#> 153  346780 0.08845351 0.7798066
#> 154  174348 0.08845351 0.7798066
#> 157  185281 0.08845351 0.7798066
#> 158  167602 0.08845351 0.7798066
#> 169  145786 0.04523193 0.5896267
#> 177  182980 0.08845351 0.7798066
#> 180  204044 0.08845351 0.7798066
#> 183   51024 0.03712191 0.4923233
#> 202  255584 0.04071219 0.5352842
#> 210  216862 0.08845351 0.7798066
#> 225  154268 0.03712191 0.4923233
#> 238  547632 0.03712191 0.4923233
#> 249  294022 0.08845351 0.7798066
#> 11   259569 0.10920535 0.8585391
#> 16   348374 0.10561819 0.8523046
#> 41   470973 0.10920535 0.8585391
#> 54   589024 0.10561819 0.8523046
#> 60   268540 0.10561819 0.8523046
#> 69   302160 0.10561819 0.8523046
#> 72   352304 0.10920535 0.8585391
#> 94     2000 0.10561819 0.8523046
#> 98   109014 0.10561819 0.8523046
#> 102  366612 0.10561819 0.8523046
#> 112  170206 0.10561819 0.8523046
#> 132  104780 0.10561819 0.8523046
#> 150  511008 0.10920535 0.8585391
#> 155  278234 0.10561819 0.8523046
#> 163  213747 0.10561819 0.8523046
#> 171  190872 0.10561819 0.8523046
#> 172  178915 0.10561819 0.8523046
#> 182  191306 0.10561819 0.8523046
#> 184  510994 0.10920535 0.8585391
#> 194   21698 0.10561819 0.8523046
#> 230  193632 0.10561819 0.8523046
#> 73   469873 0.12568184 0.9442661
#> 80   328422 0.12568184 0.9442661
#> 88   405524 0.12568184 0.9442661
#> 104  367055 0.12568184 0.9442661
#> 117  470114 0.12568184 0.9442661
#> 134  544358 0.12568184 0.9442661
#> 146  204889 0.12568184 0.9442661
#> 159  470812 0.12568184 0.9442661
#> 209  363692 0.12568184 0.9442661
#> 5    171551 0.26216189 1.0000000
#> 8    322235 0.26216189 1.0000000
#> 1    331820 1.00000000 1.0000000
#> 2    158660 1.00000000 1.0000000
#> 4    244304 0.31405659 1.0000000
#> 6    263681 0.65006789 1.0000000
#> 10   326977 1.00000000 1.0000000
#> 18   278873 0.92630130 1.0000000
#> 22   357591 0.40443382 1.0000000
#> 23   127309 1.00000000 1.0000000
#> 25   517293 0.92630130 1.0000000
#> 27   248140 1.00000000 1.0000000
#> 33   211720 1.00000000 1.0000000
#> 38   470172 0.63410287 1.0000000
#> 47   549871 0.66473219 1.0000000
#> 57   347862 0.18920661 1.0000000
#> 65   249661 0.89886959 1.0000000
#> 74   183200 0.63410287 1.0000000
#> 79    42680 0.56071994 1.0000000
#> 82   326792 0.18920661 1.0000000
#> 91   171772 0.43206312 1.0000000
#> 93   541301 0.80190256 1.0000000
#> 95   469482 0.25375509 1.0000000
#> 105  181095 0.43206312 1.0000000
#> 107  299892 0.56071994 1.0000000
#> 108  175497 0.43206312 1.0000000
#> 116  366882 1.00000000 1.0000000
#> 133  109546 0.43206312 1.0000000
#> 137  320147 0.88438586 1.0000000
#> 140  243150 0.88438586 1.0000000
#> 143  191687 0.30148798 1.0000000
#> 152  302974 0.43206312 1.0000000
#> 166  181892 0.18920661 1.0000000
#> 173  181028 0.84376712 1.0000000
#> 181  550814 0.84376712 1.0000000
#> 189  360268 0.88438586 1.0000000
#> 193  263461 0.30148798 1.0000000
#> 204  275891 0.84376712 1.0000000
#> 206  368928 0.88438586 1.0000000
#> 211  230437 0.89886959 1.0000000
#> 212  212128 0.63410287 1.0000000
#> 216  469946 0.92630130 1.0000000
#> 224  103685 1.00000000 1.0000000
#> 237  565565 0.89886959 1.0000000
#> 243  197817 0.90937751 1.0000000
#> 246  173810 1.00000000 1.0000000
#> 247  215890 0.63410287 1.0000000
#> 
#> $`WLX+KS`
#>     feature       Pval   adjPval
#> 28   256977 0.05714286 0.4230704
#> 37   309788 0.05714286 0.4230704
#> 39   240687 0.05714286 0.4230704
#> 40     3046 0.05714286 0.4230704
#> 44   573135 0.05714286 0.4230704
#> 50   220494 0.05714286 0.4230704
#> 56   306180 0.05714286 0.4230704
#> 61    27669 0.05714286 0.4230704
#> 67   565399 0.05714286 0.4230704
#> 68   111806 0.05714286 0.4230704
#> 75    65683 0.05714286 0.4230704
#> 84    71074 0.05714286 0.4230704
#> 85   509622 0.05714286 0.4230704
#> 90   360298 0.05714286 0.4230704
#> 100  223329 0.05714286 0.4230704
#> 113  313985 0.05714286 0.4230704
#> 114  574458 0.05714286 0.4230704
#> 119  247010 0.05714286 0.4230704
#> 120  154102 0.05714286 0.4230704
#> 121  129759 0.05714286 0.4230704
#> 122  113626 0.05714286 0.4230704
#> 125  273055 0.05714286 0.4230704
#> 126  551424 0.05714286 0.4230704
#> 127  262714 0.05714286 0.4230704
#> 129  340507 0.05714286 0.4230704
#> 130   84768 0.05714286 0.4230704
#> 135  158370 0.05714286 0.4230704
#> 136  566886 0.05714286 0.4230704
#> 142  216191 0.05714286 0.4230704
#> 144  352632 0.05714286 0.4230704
#> 145  245203 0.05714286 0.4230704
#> 149  265043 0.05714286 0.4230704
#> 156  242675 0.05714286 0.4230704
#> 160  559200 0.05714286 0.4230704
#> 161  367659 0.05714286 0.4230704
#> 162  161227 0.05714286 0.4230704
#> 164  277931 0.05714286 0.4230704
#> 167  527683 0.05714286 0.4230704
#> 170  113866 0.05714286 0.4230704
#> 176  356083 0.05714286 0.4230704
#> 178  200741 0.05714286 0.4230704
#> 185  308892 0.05714286 0.4230704
#> 187    1605 0.05714286 0.4230704
#> 188  250291 0.05714286 0.4230704
#> 190  218467 0.05714286 0.4230704
#> 191  286104 0.05714286 0.4230704
#> 192  511844 0.05714286 0.4230704
#> 196  543366 0.05714286 0.4230704
#> 198  535437 0.05714286 0.4230704
#> 200  159711 0.05714286 0.4230704
#> 201  144381 0.05714286 0.4230704
#> 205  584276 0.05714286 0.4230704
#> 207  570888 0.05714286 0.4230704
#> 208  278390 0.05714286 0.4230704
#> 213  170339 0.05714286 0.4230704
#> 215   93865 0.05714286 0.4230704
#> 219  251499 0.05714286 0.4230704
#> 221   70100 0.05714286 0.4230704
#> 222   97561 0.05714286 0.4230704
#> 223  278272 0.05714286 0.4230704
#> 226  130368 0.05714286 0.4230704
#> 227   76270 0.05714286 0.4230704
#> 228  190723 0.05714286 0.4230704
#> 229  208694 0.05714286 0.4230704
#> 231  222729 0.05714286 0.4230704
#> 233  203846 0.05714286 0.4230704
#> 234  156101 0.05714286 0.4230704
#> 240   13034 0.05714286 0.4230704
#> 244  109907 0.05714286 0.4230704
#> 248  326414 0.05714286 0.4230704
#> 34    36155 0.05714286 0.4230704
#> 46    47916 0.05714286 0.4230704
#> 55   279297 0.05714286 0.4230704
#> 64   150925 0.05714286 0.4230704
#> 71   578268 0.05714286 0.4230704
#> 81   208283 0.05714286 0.4230704
#> 89   546313 0.05714286 0.4230704
#> 92   154822 0.05714286 0.4230704
#> 97   565556 0.05714286 0.4230704
#> 99   100954 0.05714286 0.4230704
#> 109  166835 0.05714286 0.4230704
#> 110  150508 0.05714286 0.4230704
#> 115  113959 0.05714286 0.4230704
#> 118  544749 0.05714286 0.4230704
#> 138  306921 0.05714286 0.4230704
#> 139  114339 0.05714286 0.4230704
#> 147  203722 0.05714286 0.4230704
#> 151  184899 0.05714286 0.4230704
#> 165  235270 0.05714286 0.4230704
#> 168  151018 0.05714286 0.4230704
#> 174  252410 0.05714286 0.4230704
#> 175  331647 0.05714286 0.4230704
#> 179  244491 0.05714286 0.4230704
#> 186  211351 0.05714286 0.4230704
#> 195  154230 0.05714286 0.4230704
#> 197  549552 0.05714286 0.4230704
#> 199  265094 0.05714286 0.4230704
#> 203  215700 0.05714286 0.4230704
#> 214  582012 0.05714286 0.4230704
#> 217  204592 0.05714286 0.4230704
#> 218  252778 0.05714286 0.4230704
#> 220    1791 0.05714286 0.4230704
#> 232  206632 0.05714286 0.4230704
#> 235  575937 0.05714286 0.4230704
#> 236  349634 0.05714286 0.4230704
#> 239  553113 0.05714286 0.4230704
#> 241  210176 0.05714286 0.4230704
#> 242  308535 0.05714286 0.4230704
#> 245  113921 0.05714286 0.4230704
#> 250  534348 0.05714286 0.4230704
#> 3    189047 0.05714286 0.4230704
#> 7    192573 0.05714286 0.4230704
#> 9    180658 0.05714286 0.4230704
#> 12   325407 0.05714286 0.4230704
#> 13   361496 0.05714286 0.4230704
#> 14   316732 0.05714286 0.4230704
#> 15   291090 0.05714286 0.4230704
#> 17   196731 0.05714286 0.4230704
#> 19   290260 0.05714286 0.4230704
#> 20   361480 0.05714286 0.4230704
#> 21   193761 0.05714286 0.4230704
#> 24   357795 0.05714286 0.4230704
#> 26   223059 0.05714286 0.4230704
#> 29   287978 0.05714286 0.4230704
#> 30   357470 0.05714286 0.4230704
#> 31   367433 0.05714286 0.4230704
#> 32   187524 0.05714286 0.4230704
#> 35   518438 0.05714286 0.4230704
#> 36   327476 0.05714286 0.4230704
#> 42   199487 0.05714286 0.4230704
#> 43   191541 0.05714286 0.4230704
#> 45   365676 0.05714286 0.4230704
#> 48   269833 0.05714286 0.4230704
#> 49   194053 0.05714286 0.4230704
#> 51   190464 0.05714286 0.4230704
#> 52   288134 0.05714286 0.4230704
#> 53   194648 0.05714286 0.4230704
#> 58   138006 0.05714286 0.4230704
#> 59   212619 0.05714286 0.4230704
#> 62   293896 0.05714286 0.4230704
#> 63   195173 0.05714286 0.4230704
#> 66   192127 0.05714286 0.4230704
#> 70   193066 0.05714286 0.4230704
#> 76   177339 0.05714286 0.4230704
#> 77   193343 0.05714286 0.4230704
#> 78   203708 0.05714286 0.4230704
#> 83   261912 0.05714286 0.4230704
#> 86   175537 0.05714286 0.4230704
#> 87   557785 0.15390114 1.0000000
#> 96   179460 0.05714286 0.4230704
#> 101  181843 0.05714286 0.4230704
#> 103  310301 0.05714286 0.4230704
#> 106  162162 0.05714286 0.4230704
#> 111  190108 0.05714286 0.4230704
#> 123   82283 0.05714286 0.4230704
#> 124  262724 0.05714286 0.4230704
#> 128  111986 0.05714286 0.4230704
#> 131  250258 0.15390114 1.0000000
#> 141  188646 0.05714286 0.4230704
#> 148  469991 0.05714286 0.4230704
#> 153  346780 0.05714286 0.4230704
#> 154  174348 0.05714286 0.4230704
#> 157  185281 0.05714286 0.4230704
#> 158  167602 0.05714286 0.4230704
#> 169  145786 0.51629925 1.0000000
#> 177  182980 0.05714286 0.4230704
#> 180  204044 0.05714286 0.4230704
#> 183   51024 0.15390114 1.0000000
#> 202  255584 0.29829553 1.0000000
#> 210  216862 0.05714286 0.4230704
#> 225  154268 0.15390114 1.0000000
#> 238  547632 0.15390114 1.0000000
#> 249  294022 0.05714286 0.4230704
#> 11   259569 0.05714286 0.4230704
#> 16   348374 0.05714286 0.4230704
#> 41   470973 0.05714286 0.4230704
#> 54   589024 0.05714286 0.4230704
#> 60   268540 0.05714286 0.4230704
#> 69   302160 0.05714286 0.4230704
#> 72   352304 0.05714286 0.4230704
#> 94     2000 0.05714286 0.4230704
#> 98   109014 0.05714286 0.4230704
#> 102  366612 0.05714286 0.4230704
#> 112  170206 0.05714286 0.4230704
#> 132  104780 0.05714286 0.4230704
#> 150  511008 0.05714286 0.4230704
#> 155  278234 0.05714286 0.4230704
#> 163  213747 0.05714286 0.4230704
#> 171  190872 0.05714286 0.4230704
#> 172  178915 0.05714286 0.4230704
#> 182  191306 0.05714286 0.4230704
#> 184  510994 0.05714286 0.4230704
#> 194   21698 0.05714286 0.4230704
#> 230  193632 0.05714286 0.4230704
#> 73   469873 0.05714286 0.4230704
#> 80   328422 0.05714286 0.4230704
#> 88   405524 0.05714286 0.4230704
#> 104  367055 0.05714286 0.4230704
#> 117  470114 0.05714286 0.4230704
#> 134  544358 0.05714286 0.4230704
#> 146  204889 0.05714286 0.4230704
#> 159  470812 0.05714286 0.4230704
#> 209  363692 0.05714286 0.4230704
#> 5    171551 0.05714286 0.4230704
#> 8    322235 0.05714286 0.4230704
#> 1    331820 0.05714286 0.4230704
#> 2    158660 0.05714286 0.4230704
#> 4    244304 0.64316632 1.0000000
#> 6    263681 0.79214898 1.0000000
#> 10   326977 0.05714286 0.4230704
#> 18   278873 0.87293060 1.0000000
#> 22   357591 0.29829553 1.0000000
#> 23   127309 1.00000000 1.0000000
#> 25   517293 0.79214898 1.0000000
#> 27   248140 0.05714286 0.4230704
#> 33   211720 1.00000000 1.0000000
#> 38   470172 0.05714286 0.4230704
#> 47   549871 0.64316632 1.0000000
#> 57   347862 0.15390114 1.0000000
#> 65   249661 0.40000000 1.0000000
#> 74   183200 0.05714286 0.4230704
#> 79    42680 0.53380469 1.0000000
#> 82   326792 0.15390114 1.0000000
#> 91   171772 0.29829553 1.0000000
#> 93   541301 0.79214898 1.0000000
#> 95   469482 0.15390114 1.0000000
#> 105  181095 0.29829553 1.0000000
#> 107  299892 0.29829553 1.0000000
#> 108  175497 0.22857143 1.0000000
#> 116  366882 1.00000000 1.0000000
#> 133  109546 0.22857143 1.0000000
#> 137  320147 0.29829553 1.0000000
#> 140  243150 0.22857143 1.0000000
#> 143  191687 0.22857143 1.0000000
#> 152  302974 0.29829553 1.0000000
#> 166  181892 0.15390114 1.0000000
#> 173  181028 0.15390114 1.0000000
#> 181  550814 0.15390114 1.0000000
#> 189  360268 0.29829553 1.0000000
#> 193  263461 0.22857143 1.0000000
#> 204  275891 0.15390114 1.0000000
#> 206  368928 0.22857143 1.0000000
#> 211  230437 0.40000000 1.0000000
#> 212  212128 0.05714286 0.4230704
#> 216  469946 0.95209035 1.0000000
#> 224  103685 1.00000000 1.0000000
#> 237  565565 0.40000000 1.0000000
#> 243  197817 0.81724014 1.0000000
#> 246  173810 1.00000000 1.0000000
#> 247  215890 0.05714286 0.4230704
#> 
#> $`LR+KS`
#>     feature       Pval   adjPval
#> 28   256977 0.01571932 0.3424946
#> 37   309788 0.01571932 0.3424946
#> 39   240687 0.01571932 0.3424946
#> 40     3046 0.01571932 0.3424946
#> 44   573135 0.01571932 0.3424946
#> 50   220494 0.01571932 0.3424946
#> 56   306180 0.01571932 0.3424946
#> 61    27669 0.01571932 0.3424946
#> 67   565399 0.01571932 0.3424946
#> 68   111806 0.01571932 0.3424946
#> 75    65683 0.01571932 0.3424946
#> 84    71074 0.01571932 0.3424946
#> 85   509622 0.01571932 0.3424946
#> 90   360298 0.01571932 0.3424946
#> 100  223329 0.01571932 0.3424946
#> 113  313985 0.01571932 0.3424946
#> 114  574458 0.01571932 0.3424946
#> 119  247010 0.01571932 0.3424946
#> 120  154102 0.01571932 0.3424946
#> 121  129759 0.01571932 0.3424946
#> 122  113626 0.01571932 0.3424946
#> 125  273055 0.01571932 0.3424946
#> 126  551424 0.01571932 0.3424946
#> 127  262714 0.01571932 0.3424946
#> 129  340507 0.01571932 0.3424946
#> 130   84768 0.01571932 0.3424946
#> 135  158370 0.01571932 0.3424946
#> 136  566886 0.01571932 0.3424946
#> 142  216191 0.01571932 0.3424946
#> 144  352632 0.01571932 0.3424946
#> 145  245203 0.01571932 0.3424946
#> 149  265043 0.01571932 0.3424946
#> 156  242675 0.01571932 0.3424946
#> 160  559200 0.01571932 0.3424946
#> 161  367659 0.01571932 0.3424946
#> 162  161227 0.01571932 0.3424946
#> 164  277931 0.01571932 0.3424946
#> 167  527683 0.01571932 0.3424946
#> 170  113866 0.01571932 0.3424946
#> 176  356083 0.01571932 0.3424946
#> 178  200741 0.01571932 0.3424946
#> 185  308892 0.01571932 0.3424946
#> 187    1605 0.01571932 0.3424946
#> 188  250291 0.01571932 0.3424946
#> 190  218467 0.01571932 0.3424946
#> 191  286104 0.01571932 0.3424946
#> 192  511844 0.01571932 0.3424946
#> 196  543366 0.01571932 0.3424946
#> 198  535437 0.01571932 0.3424946
#> 200  159711 0.01571932 0.3424946
#> 201  144381 0.01571932 0.3424946
#> 205  584276 0.01571932 0.3424946
#> 207  570888 0.01571932 0.3424946
#> 208  278390 0.01571932 0.3424946
#> 213  170339 0.01571932 0.3424946
#> 215   93865 0.01571932 0.3424946
#> 219  251499 0.01571932 0.3424946
#> 221   70100 0.01571932 0.3424946
#> 222   97561 0.01571932 0.3424946
#> 223  278272 0.01571932 0.3424946
#> 226  130368 0.01571932 0.3424946
#> 227   76270 0.01571932 0.3424946
#> 228  190723 0.01571932 0.3424946
#> 229  208694 0.01571932 0.3424946
#> 231  222729 0.01571932 0.3424946
#> 233  203846 0.01571932 0.3424946
#> 234  156101 0.01571932 0.3424946
#> 240   13034 0.01571932 0.3424946
#> 244  109907 0.01571932 0.3424946
#> 248  326414 0.01571932 0.3424946
#> 34    36155 0.03185593 0.4416879
#> 46    47916 0.03185593 0.4416879
#> 55   279297 0.03185593 0.4416879
#> 64   150925 0.03185593 0.4416879
#> 71   578268 0.03185593 0.4416879
#> 81   208283 0.03185593 0.4416879
#> 89   546313 0.03185593 0.4416879
#> 92   154822 0.03185593 0.4416879
#> 97   565556 0.03185593 0.4416879
#> 99   100954 0.03185593 0.4416879
#> 109  166835 0.03185593 0.4416879
#> 110  150508 0.03185593 0.4416879
#> 115  113959 0.03185593 0.4416879
#> 118  544749 0.03185593 0.4416879
#> 138  306921 0.03185593 0.4416879
#> 139  114339 0.03185593 0.4416879
#> 147  203722 0.03185593 0.4416879
#> 151  184899 0.03185593 0.4416879
#> 165  235270 0.03185593 0.4416879
#> 168  151018 0.03185593 0.4416879
#> 174  252410 0.03185593 0.4416879
#> 175  331647 0.03185593 0.4416879
#> 179  244491 0.03185593 0.4416879
#> 186  211351 0.03185593 0.4416879
#> 195  154230 0.03185593 0.4416879
#> 197  549552 0.03185593 0.4416879
#> 199  265094 0.03185593 0.4416879
#> 203  215700 0.03185593 0.4416879
#> 214  582012 0.03185593 0.4416879
#> 217  204592 0.03185593 0.4416879
#> 218  252778 0.03185593 0.4416879
#> 220    1791 0.03185593 0.4416879
#> 232  206632 0.03185593 0.4416879
#> 235  575937 0.03185593 0.4416879
#> 236  349634 0.03185593 0.4416879
#> 239  553113 0.03185593 0.4416879
#> 241  210176 0.03185593 0.4416879
#> 242  308535 0.03185593 0.4416879
#> 245  113921 0.03185593 0.4416879
#> 250  534348 0.03185593 0.4416879
#> 3    189047 0.08845351 0.7798066
#> 7    192573 0.08845351 0.7798066
#> 9    180658 0.08845351 0.7798066
#> 12   325407 0.08845351 0.7798066
#> 13   361496 0.08845351 0.7798066
#> 14   316732 0.08845351 0.7798066
#> 15   291090 0.08845351 0.7798066
#> 17   196731 0.08845351 0.7798066
#> 19   290260 0.08845351 0.7798066
#> 20   361480 0.08845351 0.7798066
#> 21   193761 0.08845351 0.7798066
#> 24   357795 0.08845351 0.7798066
#> 26   223059 0.08845351 0.7798066
#> 29   287978 0.08845351 0.7798066
#> 30   357470 0.08845351 0.7798066
#> 31   367433 0.08845351 0.7798066
#> 32   187524 0.08845351 0.7798066
#> 35   518438 0.08845351 0.7798066
#> 36   327476 0.08845351 0.7798066
#> 42   199487 0.08845351 0.7798066
#> 43   191541 0.08845351 0.7798066
#> 45   365676 0.08845351 0.7798066
#> 48   269833 0.08845351 0.7798066
#> 49   194053 0.08845351 0.7798066
#> 51   190464 0.08845351 0.7798066
#> 52   288134 0.08845351 0.7798066
#> 53   194648 0.08845351 0.7798066
#> 58   138006 0.08845351 0.7798066
#> 59   212619 0.08845351 0.7798066
#> 62   293896 0.08845351 0.7798066
#> 63   195173 0.08845351 0.7798066
#> 66   192127 0.08845351 0.7798066
#> 70   193066 0.08845351 0.7798066
#> 76   177339 0.08845351 0.7798066
#> 77   193343 0.08845351 0.7798066
#> 78   203708 0.08845351 0.7798066
#> 83   261912 0.08845351 0.7798066
#> 86   175537 0.08845351 0.7798066
#> 87   557785 0.04071219 0.5399389
#> 96   179460 0.08845351 0.7798066
#> 101  181843 0.08845351 0.7798066
#> 103  310301 0.08845351 0.7798066
#> 106  162162 0.08845351 0.7798066
#> 111  190108 0.08845351 0.7798066
#> 123   82283 0.08845351 0.7798066
#> 124  262724 0.08845351 0.7798066
#> 128  111986 0.07100559 0.7798066
#> 131  250258 0.04071219 0.5399389
#> 141  188646 0.08845351 0.7798066
#> 148  469991 0.08845351 0.7798066
#> 153  346780 0.08845351 0.7798066
#> 154  174348 0.08845351 0.7798066
#> 157  185281 0.08845351 0.7798066
#> 158  167602 0.08845351 0.7798066
#> 169  145786 0.04295068 0.5598893
#> 177  182980 0.08845351 0.7798066
#> 180  204044 0.08845351 0.7798066
#> 183   51024 0.04071219 0.5399389
#> 202  255584 0.04295068 0.5598893
#> 210  216862 0.08845351 0.7798066
#> 225  154268 0.04071219 0.5399389
#> 238  547632 0.04071219 0.5399389
#> 249  294022 0.08845351 0.7798066
#> 11   259569 0.10920535 0.8585391
#> 16   348374 0.10561819 0.8523046
#> 41   470973 0.10920535 0.8585391
#> 54   589024 0.10561819 0.8523046
#> 60   268540 0.10561819 0.8523046
#> 69   302160 0.10561819 0.8523046
#> 72   352304 0.10920535 0.8585391
#> 94     2000 0.10561819 0.8523046
#> 98   109014 0.10561819 0.8523046
#> 102  366612 0.10561819 0.8523046
#> 112  170206 0.10561819 0.8523046
#> 132  104780 0.10561819 0.8523046
#> 150  511008 0.10920535 0.8585391
#> 155  278234 0.10561819 0.8523046
#> 163  213747 0.10561819 0.8523046
#> 171  190872 0.10561819 0.8523046
#> 172  178915 0.10561819 0.8523046
#> 182  191306 0.10561819 0.8523046
#> 184  510994 0.10920535 0.8585391
#> 194   21698 0.10561819 0.8523046
#> 230  193632 0.10561819 0.8523046
#> 73   469873 0.12568184 0.9442661
#> 80   328422 0.12568184 0.9442661
#> 88   405524 0.12568184 0.9442661
#> 104  367055 0.12568184 0.9442661
#> 117  470114 0.12568184 0.9442661
#> 134  544358 0.12568184 0.9442661
#> 146  204889 0.12568184 0.9442661
#> 159  470812 0.12568184 0.9442661
#> 209  363692 0.12568184 0.9442661
#> 5    171551 0.26216189 1.0000000
#> 8    322235 0.26216189 1.0000000
#> 1    331820 1.00000000 1.0000000
#> 2    158660 1.00000000 1.0000000
#> 4    244304 0.32663047 1.0000000
#> 6    263681 0.41920953 1.0000000
#> 10   326977 1.00000000 1.0000000
#> 18   278873 0.93092480 1.0000000
#> 22   357591 0.30148798 1.0000000
#> 23   127309 1.00000000 1.0000000
#> 25   517293 0.91074697 1.0000000
#> 27   248140 1.00000000 1.0000000
#> 33   211720 1.00000000 1.0000000
#> 38   470172 0.63410287 1.0000000
#> 47   549871 0.67772220 1.0000000
#> 57   347862 0.30148798 1.0000000
#> 65   249661 0.89886959 1.0000000
#> 74   183200 0.63410287 1.0000000
#> 79    42680 0.67772220 1.0000000
#> 82   326792 0.30148798 1.0000000
#> 91   171772 0.56071994 1.0000000
#> 93   541301 0.67772220 1.0000000
#> 95   469482 0.43206312 1.0000000
#> 105  181095 0.56071994 1.0000000
#> 107  299892 0.43206312 1.0000000
#> 108  175497 0.43206312 1.0000000
#> 116  366882 0.93092480 1.0000000
#> 133  109546 0.43206312 1.0000000
#> 137  320147 0.89886959 1.0000000
#> 140  243150 0.88438586 1.0000000
#> 143  191687 0.30148798 1.0000000
#> 152  302974 0.56071994 1.0000000
#> 166  181892 0.30148798 1.0000000
#> 173  181028 0.88438586 1.0000000
#> 181  550814 0.88438586 1.0000000
#> 189  360268 0.89886959 1.0000000
#> 193  263461 0.30148798 1.0000000
#> 204  275891 0.88438586 1.0000000
#> 206  368928 0.88438586 1.0000000
#> 211  230437 0.89886959 1.0000000
#> 212  212128 0.63410287 1.0000000
#> 216  469946 0.96381500 1.0000000
#> 224  103685 0.93092480 1.0000000
#> 237  565565 0.89886959 1.0000000
#> 243  197817 0.93092480 1.0000000
#> 246  173810 0.93092480 1.0000000
#> 247  215890 0.63410287 1.0000000
```

Because the core method is `NULL`, this example combines only enhancer
tests. The `$ensembles` element reports the CCT results for the enhancer
combinations when `return_subensembles = TRUE`.

## Single-Cell Example

The single-cell example uses a small `SingleCellExperiment` object with
two cell groups. The object is wrapped in a `MultiAssayExperiment`
before being passed to `DAssemble`, matching the same input architecture
used above.

``` r

set.seed(1)
n_genes <- 120L
n_cells <- 20L

cell_type <- factor(rep(c("control", "treated"), each = n_cells / 2L))
counts <- matrix(
  stats::rnbinom(n_genes * n_cells, mu = 20, size = 5),
  nrow = n_genes,
  dimnames = list(
    paste0("gene", seq_len(n_genes)),
    paste0("cell", seq_len(n_cells))
  )
)

counts[seq_len(12L), cell_type == "treated"] <-
  counts[seq_len(12L), cell_type == "treated"] + 15L

metadata <- data.frame(
  CellType = cell_type,
  row.names = colnames(counts)
)

sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = counts),
  colData = S4Vectors::DataFrame(metadata)
)

mae <- MultiAssayExperiment::MultiAssayExperiment(
  experiments = list(single_cell = sce),
  colData = S4Vectors::DataFrame(metadata)
)

res <- DAssemble::DAssemble(
  features = mae,
  assay_name = "single_cell",
  core_method = "edgeR",
  enhancers = c("WLX", "LR", "KS"),
  enhancer_norm = "TSS",
  expVar = "CellType",
  p_adj = "BH",
  return_components = TRUE,
  return_subensembles = TRUE
)
#> calcNormFactors has been renamed to normLibSizes

head(res$res)
#>   feature metadata    pval_core    pval_WLX coef_LR pval_LR      pval_KS
#> 1   gene1 CellType 0.1319973332 0.089209552      NA      NA 0.1678213427
#> 2   gene2 CellType 0.0127052773 0.006841456      NA      NA 0.0123406006
#> 3   gene3 CellType 0.0610776934 0.089209552      NA      NA 0.1678213427
#> 4   gene4 CellType 0.0001603293 0.001050034      NA      NA 0.0002165018
#> 5   gene5 CellType 0.0791618751 0.063012839      NA      NA 0.0524475524
#> 6   gene6 CellType 0.0032635317 0.008930698      NA      NA 0.0123406006
#>   pval_joint qval
#> 1          1    1
#> 2          1    1
#> 3          1    1
#> 4          1    1
#> 5          1    1
#> 6          1    1
names(res$components)
#> [1] "edgeR" "WLX"   "LR"    "KS"
```

## Session information

``` r

sessionInfo()
#> R version 4.6.2 (2026-10-29)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.5 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] BiocStyle_2.40.0
#> 
#> loaded via a namespace (and not attached):
#>   [1] Rdpack_2.6.6                DBI_1.3.0                  
#>   [3] pbapply_1.7-5               permute_0.10-0             
#>   [5] sandwich_3.1-3              rlang_1.3.0                
#>   [7] magrittr_2.0.5              ade4_1.7-24                
#>   [9] otel_0.3.0                  matrixStats_1.5.0          
#>  [11] compiler_4.6.2              mgcv_1.9-4                 
#>  [13] reshape2_1.4.5              systemfonts_1.3.2          
#>  [15] vctrs_0.7.3                 stringr_1.6.0              
#>  [17] pkgconfig_2.0.3             crayon_1.5.3               
#>  [19] fastmap_1.2.0               XVector_0.52.0             
#>  [21] rmarkdown_2.32              nloptr_2.2.1               
#>  [23] ragg_1.5.2                  xfun_0.62                  
#>  [25] MultiAssayExperiment_1.38.0 cachem_1.1.0               
#>  [27] jsonlite_2.0.0              biomformat_1.40.0          
#>  [29] DelayedArray_0.38.2         BiocParallel_1.46.0        
#>  [31] RTMB_2.0                    parallel_4.6.2             
#>  [33] Maaslin2_1.26.0             cluster_2.1.8.3            
#>  [35] R6_2.6.1                    biglm_0.9-3                
#>  [37] stringi_1.8.9               bslib_0.12.0               
#>  [39] RColorBrewer_1.1-3          limma_3.68.5               
#>  [41] boot_1.3-32                 GenomicRanges_1.64.0       
#>  [43] jquerylib_0.1.4             numDeriv_2016.8-1.1        
#>  [45] iterators_1.0.14            Rcpp_1.1.2                 
#>  [47] Seqinfo_1.2.0               bookdown_0.48              
#>  [49] SummarizedExperiment_1.42.0 knitr_1.52                 
#>  [51] zoo_1.9-1                   IRanges_2.46.0             
#>  [53] igraph_2.3.4                Matrix_1.7-6               
#>  [55] splines_4.6.2               glmmTMB_1.1.15.2           
#>  [57] tidyselect_1.2.1            abind_1.4-8                
#>  [59] yaml_2.3.12                 vegan_2.7-6                
#>  [61] TMB_1.9.25                  codetools_0.2-20           
#>  [63] plyr_1.8.9                  lattice_0.23-1             
#>  [65] tibble_3.3.1                lmerTest_3.2-1             
#>  [67] Biobase_2.72.0              withr_3.0.3                
#>  [69] S7_0.2.2                    evaluate_1.0.5             
#>  [71] survival_3.8-12             desc_1.4.3                 
#>  [73] Biostrings_2.80.2           phyloseq_1.56.0            
#>  [75] pillar_1.11.1               BiocManager_1.30.27        
#>  [77] MatrixGenerics_1.24.0       foreach_1.5.2              
#>  [79] stats4_4.6.2                reformulas_0.4.4           
#>  [81] pcaPP_2.0-7                 generics_0.1.4             
#>  [83] S4Vectors_0.50.3            ggplot2_4.0.3              
#>  [85] scales_1.4.0                minqa_1.2.8                
#>  [87] DAssemble_0.99.4            glue_1.8.1                 
#>  [89] tools_4.6.2                 data.table_1.18.6.1        
#>  [91] robustbase_0.99-7           lme4_2.1-0                 
#>  [93] locfit_1.5-9.12             fs_2.1.0                   
#>  [95] mvtnorm_1.4-2               grid_4.6.2                 
#>  [97] optparse_1.8.2              ape_5.8-1                  
#>  [99] rbibutils_2.4.1             SingleCellExperiment_1.34.0
#> [101] edgeR_4.10.7                nlme_3.1-171               
#> [103] cli_3.6.6                   textshaping_1.0.5          
#> [105] S4Arrays_1.12.1             dplyr_1.2.1                
#> [107] gtable_0.3.6                DEoptimR_1.2-2             
#> [109] logging_0.10-111            hash_2.2.6.4               
#> [111] DESeq2_1.52.0               sass_0.4.10                
#> [113] digest_0.6.39               BiocGenerics_0.58.1        
#> [115] SparseArray_1.12.3          farver_2.1.2               
#> [117] multtest_2.68.0             htmltools_0.5.9            
#> [119] pkgdown_2.2.1               lifecycle_1.0.5            
#> [121] statmod_1.5.2               MASS_7.3-66
```
