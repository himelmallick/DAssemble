# DAssemble: Ensemble Models for Differential Analysis

DAssemble provides an ensemble framework for differential-abundance and
differential-expression analysis across bulk RNA-seq, single-cell
RNA-seq, and microbiome studies. The package combines a user-selected
core method with optional enhancer tests and aggregates evidence using
the Cauchy Combination Test.

## Details

The main user-facing function is
[`DAssemble()`](https://himelmallick.github.io/DAssemble/reference/DAssemble.md).
It accepts either a `MultiAssayExperiment` object or a feature table
paired with sample metadata. Core methods can be combined with enhancer
tests including Wilcoxon rank-sum, Kolmogorov-Smirnov, and
presence-absence logistic regression.

## See also

[`DAssemble()`](https://himelmallick.github.io/DAssemble/reference/DAssemble.md)

## Author

**Maintainer**: Nalin Arora <naa4050@med.cornell.edu>
([ORCID](https://orcid.org/0009-0009-1340-688X))

Authors:

- Nalin Arora <naa4050@med.cornell.edu>
  ([ORCID](https://orcid.org/0009-0009-1340-688X))

- Ziyu Liu <lclcencen@gmail.com>

- Chuxuan Gao <chg4016@alumni.weill.cornell.edu>

- Himel Mallick <him4004@med.cornell.edu>
  ([ORCID](https://orcid.org/0000-0003-4956-2429))
