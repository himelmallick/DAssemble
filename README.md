# DAssemble

[![R build status](https://github.com/Ziyu-Liu-WCM/DAssemble/workflows/check-bioc/badge.svg)](https://github.com/Ziyu-Liu-WCM/DAssemble/actions)

DAssemble is a lightweight implementation of the stacking method as
applied to two or more differential analysis (DA) results tables. It
takes as inputs a list of DA tables with p-values and returns a single
table with omnibus p-values and q-values on a per-feature basis. Several
p-value combination methods such as Stouffer’s method, CCT, and Fisher
are supported.

## Installation

To install the latest version of DAssemble from Github, run the
following command:

    install.packages('devtools')
    library(devtools)
    devtools::install_github("himelmallick/DAssemble")
    library(DAssemble)

## Input

DAssemble requires a combined list of dataframes **dflist**, each
containing at least two columns: **ID** for gene ID and **pvalue** for
DA p-values corresponding to each gene. Only genes that appear in all
dataframes will be combined.

## Output

A data frame containing gene ID, p-values, combined p-values and
q-values (multiplicity-adjusted p-values) are returned.

## Example Usage

    # Install packages if not already installed (uncomment if needed)
    # install.packages('tidyverse')
    # BiocManager::install("airway")
    # BiocManager::install("edgeR")
    # BiocManager::install("DESeq2")
    # BiocManager::install("limma")

    library(tidyverse)
    library(airway)
    library(edgeR)
    library(DESeq2)
    library(limma)
    library(DAssemble)

    #############
    # Load Data #
    #############
    # Load the airway SummarizedExperiment object
    data("airway")

    # Filter the data for genes with at least 10 reads in at least 10 samples.
    filter <- filterByExpr(airway)
    filtered <- airway[filter,]

    # Extract the count matrix
    counts_full <- assay(filtered)
    # Extract the sample information (dex = treated or untreated)
    sample_info <- colData(filtered)

    # Set a seed for reproducibility, then randomly choose 1000 genes
    set.seed(123)
    genes_subset <- sample(rownames(counts_full), size = 1000)

    # Subset the counts
    counts_sub <- counts_full[genes_subset, ]
    treatment <- factor(sample_info$dex, levels = c("untrt", "trt"))


    #############
    # Run edgeR #
    #############
    # Create a DGEList
    dge <- DGEList(counts = counts_sub, group = treatment)

    # Calculate normalization factors
    dge <- calcNormFactors(dge)
    # Create a design matrix (~ group)
    design_edgeR <- model.matrix(~ treatment)


    # Estimate dispersion
    dge <- estimateDisp(dge, design_edgeR)

    # Fit a quasi-likelihood negative binomial generalized log-linear model
    fit_edgeR <- glmQLFit(dge, design_edgeR)

    # Test for differential expression
    # coef = 2 corresponds to the effect of "trt" (since untrt is the baseline)
    qlf_edgeR <- glmQLFTest(fit_edgeR, coef = 2)

    # Get all results as a data frame
    res_edgeR <- topTags(qlf_edgeR, n = Inf)
    res_edgeR_df <- as.data.frame(res_edgeR)



    ##############
    # Run DESeq2 #
    ##############
    # Create a DESeqDataSet from the subset count matrix
    dds <- DESeqDataSetFromMatrix(
      countData = counts_sub,
      colData   = as.data.frame(sample_info),
      design    = ~ dex
    )


    # Run DESeq (this estimates size factors, dispersion, and fits models)
    dds <- DESeq(dds)


    # Extract results for dex effect: "trt" vs "untrt"
    res_DESeq2 <- results(dds, contrast = c("dex", "trt", "untrt"))
    res_DESeq2_df <- as.data.frame(res_DESeq2)


    #################
    # Run limmavoom #
    #################
    design_limma <- model.matrix(~ treatment)

    # Convert counts to log2-cpm with precision weights using voom
    v <- voom(counts_sub, design_limma, plot = FALSE)

    # Fit linear model
    fit_limma <- lmFit(v, design_limma)

    # Empirical Bayes moderation
    fit_limma <- eBayes(fit_limma)


    # coef = 2 corresponds to "trt" vs "untrt" if the baseline is untrt
    res_limma <- topTable(fit_limma, coef = 2, n = Inf)

    # Convert to data frame if needed
    res_limma_df <- as.data.frame(res_limma)



    ###############################################
    # Convert the results to desired input format #
    ###############################################
    res_edgeR_df <- res_edgeR_df %>% 
      rownames_to_column("ID") %>% 
      dplyr::rename(pvalue = PValue)

    res_DESeq2_df <- res_DESeq2_df %>% 
      rownames_to_column("ID")

    res_limma_df <- res_limma_df %>% 
      rownames_to_column("ID") %>%
      dplyr::rename(pvalue = P.Value)

    ################
    # Run DAssemble#
    ################
    dflist <- list(res_edgeR_df, res_DESeq2_df, res_limma_df)
    assemble_res <- DAssemble(dflist, combine.method = "stouffer", correction = "BH")
    head(assemble_res)

    ##                  ID      pvalue1      pvalue2      pvalue3 pval.combined
    ## 544 ENSG00000162616 1.428166e-07 1.867256e-38 1.519030e-07  4.168862e-41
    ## 291 ENSG00000123562 9.479657e-07 5.606138e-32 5.402599e-06  9.411785e-34
    ## 669 ENSG00000172986 1.472552e-06 1.432363e-23 3.364826e-07  6.194061e-30
    ## 729 ENSG00000181061 2.048669e-06 7.628886e-22 2.669372e-06  1.898282e-27
    ## 817 ENSG00000196975 4.180457e-06 1.352799e-19 7.865789e-06  5.900290e-25
    ## 320 ENSG00000127824 7.971791e-06 7.738274e-17 9.705530e-06  1.210939e-22
    ##          qval_BH
    ## 544 4.168862e-38
    ## 291 4.705892e-31
    ## 669 2.064687e-27
    ## 729 4.745704e-25
    ## 817 1.180058e-22
    ## 320 2.018232e-20
