# DAssemble

DAssemble is an R package for ensemble differential-abundance / differential-expression (DA/DE) analysis across multiple omics domains. It wraps a collection of existing DA/DE methods as **core methods**, combines them via **Cauchy Combination Tests (CCT)**, and augments them with lightweight **enhancers**.

---

## Installation

### 1. Install the DAssemble package

You can install **DAssemble** directly from GitHub:

```r
# install.packages("remotes")  # if not already installed
remotes::install_github("himelmallick/DAssemble")
```

After installation, load the package with:

```r
library(DAssemble)
```

---

### 2. Install and load all method dependencies

DAssemble depends on several packages from **CRAN**, **Bioconductor**, and **GitHub**
for its core DA methods and enhancers. The script below will automatically
install any missing packages and then load them.

```r
## ===========================================
## Install + load all required DAssemble deps
## ===========================================

# All required packages
req_pkgs <- c(
  "MAST", "DESeq2", "edgeR", "limma", "metagenomeSeq",
  "dearseq", "SummarizedExperiment", "ALDEx2", "LinDA",
  "LOCOM", "Maaslin2", "maaslin3", "Tweedieverse",
  "Robseq", "ANCOMBC", "TreeSummarizedExperiment", "S4Vectors"
)

# GitHub-only packages
github_pkgs <- c(
  LinDA        = "zhouhj1994/LinDA",
  LOCOM        = "yijuanhu/LOCOM",
  Tweedieverse = "himelmallick/Tweedieverse",
  maaslin3     = "biobakery/maaslin3",
  Robseq       = "schatterjee30/Robseq"
)

# Install required managers
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")

# Install + load loop
for (pkg in req_pkgs) {

  # Install if missing
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing missing package: ", pkg)

    if (pkg %in% names(github_pkgs)) {
      repo <- github_pkgs[[pkg]]
      message("  -> Installing from GitHub: ", repo)
      remotes::install_github(repo, dependencies = TRUE)
    } else {
      tryCatch({
        install.packages(pkg)
      }, error = function(e) {
        message("  -> CRAN failed, installing via Bioconductor: ", pkg)
        BiocManager::install(pkg, ask = FALSE)
      })
    }
  }

  # Load package
  message("Loading package: ", pkg)
  ok <- suppressPackageStartupMessages(
    require(pkg, character.only = TRUE)
  )
  if (!ok)
    stop("Failed to load package: ", pkg, call. = FALSE)
}

cat("All required packages installed and loaded successfully.\n")
```


---

## Core Methods Supported

| Method       | Package(s)                      |
|-------------|----------------------------------|
| DESeq2      | DESeq2                          |
| edgeR       | edgeR                           |
| limma-voom  | limma + edgeR                   |
| metagenomeSeq | metagenomeSeq                 |
| MAST        | MAST                            |
| dearseq     | dearseq + SummarizedExperiment  |
| ALDEx2      | ALDEx2                          |
| LinDA       | LinDA (GitHub)                  |
| LOCOM      | LOCOM (GitHub)                  |
| Maaslin2    | Maaslin2                        |
| Maaslin3    | maaslin3 (GitHub)               |
| Tweedieverse | Tweedieverse (GitHub)          |
| Robseq      | Robseq (GitHub)                 |
| ANCOM-BC2   | ANCOMBC + TreeSummarizedExperiment + S4Vectors |



### Note on singleton core methods

Some methods already **combine two models internally** (e.g., hurdle or two‑part models) and/or directly output a **combined p‑value**.  These should normally be used **alone** in an ensemble, i.e. as the *only* core method.  Combining a hurdle model with other core methods can double‑count the same underlying signal and lead to inflated type I error.  The following methods fall into this category:

* `Maaslin3` – fits both a prevalence and abundance model under the hood;
* `MAST` – implements a hurdle model for single‑cell expression data;
* `metagenomeSeq` – fits a zero‑inflated log‑normal model.

You may still add enhancers to these models (e.g. to perform sensitivity analyses), but combining them with additional core methods is not recommended.

---

## Compatible combinations

DAssemble can be run in a number of configurations depending on whether you include a core method, enhancers, or both.  The current framework allows you
to combine one core method with up to **three** enhancers, or to run enhancer‑only analyses when `core_method = NULL`.  The table below lists
the valid combinations of enhancers by the number of enhancers used.  The currently available enhancers are `WLX` (Wilcoxon rank‑sum), `LR` logistic regression) and `KS` (Kolmogorov–Smirnov).

| # enhancers | With core | Example combinations |
|-------------|-----------|-----------------------|
| **0**       | Yes       | `core_method = "DESeq2"`, `enhancers = NULL` |
| **1**       | Yes / No  | `c("WLX")`, `c("LR")`, `c("KS")` |
| **2**       | Yes / No  | `c("WLX","LR")`, `c("WLX","KS")`, `c("LR","KS")` |
| **3**       | Yes / No  | `c("WLX","LR","KS")` |

For example, you might call:

```r
# core + two enhancers
DAssemble(X, metadata, core_method = "DESeq2", enhancers = c("WLX", "LR"), expVar = "group")

# enhancer‑only, all three
DAssemble(X, metadata, core_method = NULL, enhancers = c("WLX", "LR", "KS"), expVar = "group")
```

When `return_subensembles = TRUE` the returned `$ensembles` element will
include every sub‑combination of the requested methods.  For instance,
with three enhancers and no core, the sub‑ensembles will report CCT
results for each single test (`"WLX"`, `"LR"`, `"KS"`), each pair
(`"WLX+LR"`, etc.) and the full trio (`"WLX+LR+KS"`).

---

## Enhancers

DAssemble currently implements **three enhancer methods**:

- `DA_fit_enhancer_WLX()` – Wilcoxon rank-sum test (**WLX**)
- `DA_fit_enhancer_LR()` – presence–absence **logistic regression (LR)**
- `DA_fit_enhancer_KS()` – Kolmogorov–Smirnov test (**KS**)

All enhancers return:

```r
feature     # character
pval_<TAG>  # numeric, where <TAG> is WLX / LR / KS
```

Within `DAssemble()`, you usually specify enhancers using their short tags:

```r
enhancers = c("WLX", "LR", "KS")  # at most two at a time are currently supported
```

DAssemble then routes to the corresponding `DA_fit_enhancer_*()` functions.


---

## DAssemble main function

The main entry point is:

```r
DAssemble(
  X,
  expVar,
  coVars        = NULL,
  core_method   = "DESeq2",
  enhancers     = c("WLX", "KS"),
  domain        = c("bulkrnaseq", "singlecell", "microbiome", "none"),
  p_adj         = "BH",
  orientation   = c("features_by_samples", "samples_by_features", "auto"),
  return_components   = FALSE,
  return_subensembles = FALSE,
  ...
)
```

Key arguments:

- `X` – feature matrix of counts / abundances  
- `expVar` – name (or column) of the primary exposure variable (binary, 2 levels)
- `coVars` – optional covariates (used by some core methods)
- `core_method` – one of the supported core method names (see above), or `NULL` / `"none"` to run an **enhancer-only** analysis
- `enhancers` – `NULL` or a subset of `c("WLX", "LR", "KS")`
- `domain` – `"bulkrnaseq"`, `"singlecell"`, `"microbiome"`, or `"none"`; **used only by the enhancers** to choose modality-specific normalization
- `p_adj` – multiple testing correction method (passed to `p.adjust`)
- `orientation` – `"features_by_samples"`, `"samples_by_features"`, or `"auto"`
- `return_components` – if `TRUE`, return per-method results in `$components`
- `return_subensembles` – if `TRUE`, compute CCT **sub-ensembles** and return them in `$ensembles`

The main output includes:

- `$results` – data frame of features with ensemble p-values (raw and adjusted)
- `$components` – (optional) list of per-method result tables
- `$ensembles` – (optional) table of sub-ensemble CCT results

---

## Sub-ensembles and `core_method = NULL`

If `return_subensembles = TRUE` and the helper `DA_build_subensembles()` is available, DAssemble computes **sub-ensembles** in addition to the main joint CCT p-values.

- When **`core_method` is not `NULL`**, sub-ensembles typically correspond to:
  - the core alone,
  - each core + single-enhancer combination,
  - and (optionally) the full core + all-enhancers ensemble.

- When **`core_method = NULL`** (enhancer-only mode), the sub-ensembles are defined purely in terms of the enhancers.  
  For example, with `enhancers = c("WLX", "LR", "KS")`, the sub-ensembles can include:
  - WLX alone, LR alone, KS alone
  - WLX + LR, WLX + KS, LR + KS
  - WLX + LR + KS

In the resulting `$ensembles` object, the `sub_ensemble` (or similarly named) column labels each combination (e.g., `"WLX"`, `"WLX+LR"`, `"LR+KS"`, `"WLX+LR+KS"`).  

This is useful for sensitivity analysis: you can see how conclusions change as you move from individual tests to different CCT combinations, including the **enhancer-only** setting when no core method is used.

---

## Real data examples

To illustrate how to use DAssemble on real datasets, this section provides
two short examples: one for a bulk RNA‑Seq experiment and one for a
microbiome study.  These examples rely on publicly available data from
Bioconductor packages; citations are provided for a brief description of
each dataset.

### Bulk RNA‑Seq: Airway dexamethasone experiment

The **airway** dataset from the Bioconductor package `airway` is a
small RNA‑Seq experiment in which four airway smooth muscle cell lines
are treated with the asthma medication dexamethasone.  It is often used
as an introductory example for differential expression analysis; the
dataset comprises counts for ~63 k genes measured in eight samples
(`4 × 2` design).  A Bioconductor vignette notes that the airway
dataset provides a typical small‑scale RNA‑Seq experiment where four
ASM cell lines are treated with dexamethasone.

To run DAssemble on this dataset using the DESeq2 core and two enhancers (Wilcoxon and logistic regression):

```r
library(DESeq2)
library(airway)
library(DAssemble)

# load counts and metadata
data("airway")
counts   <- assay(airway, "counts")
metadata <- as.data.frame(colData(airway))

# the exposure variable is dex (treated vs untreated)
res <- DAssemble(
  features    = t(counts),
  metadata    = metadata,
  core_method = "DESeq2",
  enhancers   = c("WLX", "LR"),
  expVar      = "dex",
  p_adj       = "BH",
  enhancer_norm = "tmm",
  return_components = TRUE,
  return_subensembles = TRUE
)

# inspect the top differential genes
head(res$res)
# inspect per‑method p‑values
names(res$components)
```

This call normalizes the counts using the TMM scheme (`enhancer_norm = "tmm"`),
combines DESeq2 with the Wilcoxon and logistic regression enhancers, and
returns the full set of sub‑ensembles.  You can access individual
method outputs through the `$components` element of the result.

### Microbiome: Global Patterns study

The **Global Patterns** dataset, available via the `phyloseq` package,
contains 16S rRNA profiles for 25 environmental samples and three
synthetic “mock communities,” representing nine sample types in total at
an average sequencing depth of 3.1 million reads per sample【156744987103578†L128-L133】.
It has been used to explore diversity patterns across a variety of
ecosystems.

To run an enhancer‑only DAssemble analysis on a microbiome count table
extracted from `GlobalPatterns`:

```r
library(phyloseq)
library(DAssemble)

# load Global Patterns as a phyloseq object
data("GlobalPatterns")
gp <- GlobalPatterns

# extract count matrix and sample metadata
otu  <- as(otu_table(gp), "matrix")
meta <- as.data.frame(sample_data(gp))

# here we compare human stool samples against soil samples as an example
keep <- meta$SampleType %in% c("Stool", "Soil")
X    <- t(otu[, keep])
meta <- meta[keep, , drop = FALSE]
meta$group <- droplevels(factor(meta$SampleType))

res <- DAssemble(
  features    = X,
  metadata    = meta,
  core_method = NULL,
  enhancers   = c("WLX", "LR", "KS"),
  expVar      = "group",
  enhancer_norm = "clr",
  return_components = TRUE,
  return_subensembles = TRUE
)

head(res$res)
```

In this example the core is set to `NULL`, so the analysis combines
three nonparametric enhancers only.  Because microbiome data are
compositional, we use centered log‑ratio (CLR) normalization
(`enhancer_norm = "clr"`).  The `$components` element contains the
individual Wilcoxon, logistic regression and KS results, while
`$ensembles` summarises every sub‑combination of the three enhancers.


