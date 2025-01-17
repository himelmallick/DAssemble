# DAssemble
DAssemble is a lightweight implementation of the stacking method as applied to two or more differential analysis (DA) results tables. It takes as inputs a list of DA tables with p-values and returns a single table with omnibus p-values and q-values on a per-feature basis. Several p-value combination methods such as Stouffer's method, CCT, and Fisher are supported.

## Installation
To install the latest version of DAssemble from Github, run the following command:

``` r
install.packages('devtools')
library(devtools)
devtools::install_github("himelmallick/DAssemble")
library(DAssemble)
```

After installing `DAssemble`, please make sure the following package versions are also installed: "poolr", "metapod", "harmonicmeanp", which can be done as follows:
``` r
install.packages('poolr')
install.packages('metapod')
install.packages('harmonicmeanp')
```

## Input 
DAssemble requires a combined list of dataframes **dflist**, each containing at least two columns: **ID** for gene ID and **pvalue** for DA p-values corresponding to each gene. Only genes that appear in all dataframes will be combined.


## Output
A data frame containing gene ID, p-values, combined p-values and q-values (multiplicity-adjusted p-values) are returned.


## Example Usage
``` r
pvec1 <- data.frame(ID = c('gene1', 'gene2', 'gene3', 'gene4'),
                    pvalue = c(0.01, 0.09, 0.02, 0.06))

pvec2 <- data.frame(ID = c('gene1', 'gene2', 'gene3', 'gene4'),
                    pvalue = c(0.02, 0.03, 0.04, 0.07))

pvec3 <- data.frame(ID = c('gene1', 'gene2', 'gene3', 'gene4'),
                    pvalue = c(0.05, 0.01, 0.03, 0.08))

dflist <- list(pvec1, pvec2, pvec3)

DAssemble(dflist, combine.method = 'harmonic', correction = 'hommel')
```

