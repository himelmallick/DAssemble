# DAssemble

DAssemble is an R framework for ensemble differential-abundance / differential-expression (DA/DE) analysis across multiple omics domains. It wraps a collection of existing DA/DE methods as **core methods**, combines them via **Cauchy Combination Tests (CCT)**, and augments them with lightweight **enhancers** that operate on normalized data.

This README describes:

- The **core methods** supported by DAssemble
- The three **enhancers**
- Package **dependencies** and required installations
- How **modality-specific normalization** is handled
- How **sub-ensembles** work, including the case `core_method = NULL`

---

## Implemented core methods

DAssemble currently supports the following **core methods**:

- `DESeq2`
- `edgeR`
- `limmaVOOM`
- `dearseq`
- `metagenomeSeq`
- `MAST`
- `Tweedieverse`
- `Maaslin2`
- `Maaslin3`
- `LOCOM`
- `LinDA`
- `ANCOMBC2`
- `Robseq`
- `ALDEx2`

Each core method has a corresponding wrapper function, for example:

- `DA_fit_core_DESeq2()` – DESeq2 for bulk / single-cell RNA-seq
- `DA_fit_core_edgeR()` – edgeR GLM pipeline
- `DA_fit_core_limmaVOOM()` – limma-voom for RNA-seq
- `DA_fit_core_dearseq()` – dearseq for RNA-seq
- `DA_fit_core_metagenomeSeq()` – metagenomeSeq for microbiome
- `DA_fit_core_MAST()` – MAST for single-cell data
- `DA_fit_core_Tweedieverse()` – Tweedieverse for count data
- `DA_fit_core_Maaslin2()` – MaAsLin2 for microbiome
- `DA_fit_core_Maaslin3()` – MaAsLin3 (abundance model)
- `DA_fit_core_LOCOM()` – LOCOM for microbiome
- `DA_fit_core_LinDA()` – LinDA for compositional microbiome data
- `DA_fit_core_ANCOMBC2()` – ANCOM-BC2
- `DA_fit_core_ALDEx2()` – ALDEx2
- `DA_fit_core_Robseq()` – Robseq

Each function returns a `data.frame` with at least:

```r
feature    # character
pval_core  # numeric
## plus any extra, method-specific columns
```

### Note on singleton core methods

Some methods already **combine two models internally** (e.g., hurdle or two-part models) and/or directly output a **combined p-value**. These should normally be used as **singletons**—that is, as the only core method in a CCT ensemble:

- `Maaslin3`
- `MAST`
- `metagenomeSeq`

You can still pair these with enhancers if desired, but you should avoid combining them with multiple additional core methods in ways that double-count the same underlying model components.

---

## Recommended core methods and enhancers by domain

The table below summarizes which core methods are most commonly used in DAssemble and the default enhancer choices by domain (matching the LaTeX table in the manuscript).

| Domain                | Core methods                                                                   | Enhancer 1                               | Enhancer 2 |
|-----------------------|--------------------------------------------------------------------------------|------------------------------------------|-----------|
| **Bulk RNA-Seq**      | DESeq2, edgeR, limmaVOOM, Robseq, dearseq                                     | Wilcoxon rank-sum test (**WLX**)         | KS        |
| **Single-cell RNA-Seq** | DESeq2, edgeR, limmaVOOM, Tweedieverse                                      | Presence–absence logistic regression (**LR**) | KS   |
| **Microbiome**        | ALDEx2, ANCOM-BC2, DESeq2, LinDA, LOCOM, MaAsLin2, MaAsLin3, metagenomeSeq    | Presence–absence logistic regression (**LR**) | KS   |

> **Dependencies:** DAssemble does **not** install these methods for you.  
> To use a given `core_method`, you must have the corresponding package(s) installed (e.g., **DESeq2**, **edgeR**, **limma**, **dearseq**, **metagenomeSeq**, **MAST**, **Tweedieverse**, **Maaslin2**, **Maaslin3**, **LOCOM**, **LinDA**, **ANCOMBC**, **ALDEx2**, **RobSeq**, etc.).  
> If a required package is missing, the corresponding wrapper will error with a clear message.

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

## Modality-specific normalization for enhancers

Enhancers operate on a **normalized version of the input features**.  
By default, DAssemble applies **modality-specific normalization** based on the `domain` argument:

- `domain = "bulkrnaseq"` – normalization appropriate for bulk RNA-seq  
- `domain = "singlecell"` – normalization appropriate for single-cell RNA-seq  
- `domain = "microbiome"` – normalization appropriate for microbiome data  
- `domain = "none"` – no domain-specific normalization

If you prefer to handle all preprocessing yourself, you can **turn off the internal normalization** in the enhancer wrappers:

> Set `normalization = NULL` and pass in already normalized features to the enhancer(s).  
> In this case, no additional modality-specific normalization is applied.

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


