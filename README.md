# NullstrapDE

**NullstrapDE** is an R package providing an *add-on* calibration framework for popular differential expression (DE) methods such as **DESeq2** and **edgeR**. It implements a **synthetic-null-data-based FDR calibration** procedure to improve false discovery rate (FDR) control in settings where standard DESeq2 or edgeR may fail.

The core of Nullstrap-DE consists of three main components:

1.  **Synthetic null data generation**\
    Synthetic counts are generated using fitted model parameters under gene-specific null hypotheses.

2.  **Computation of test statistics**\
    Test statistics are computed on both the original data and synthetic-null data.

3.  **FDR calibration**\
    Test statistics from synthetic-null datasets are used to calibrate a significance threshold for the original data.

This package provides two wrapper functions:

-   `Nullstrap_DESeq2()` -- NullstrapDE calibration for DESeq2

-   `Nullstrap_edgeR()` -- NullstrapDE calibration for edgeR

Both are designed as add-on extensions to standard workflows.

------------------------------------------------------------------------

## Installation<a name="installation-"></a>

To install the development version from GitHub, please run:

``` r
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("chexjiang/NullstrapDE")
```

------------------------------------------------------------------------

## Quick Start<a name="quick-start"></a>

Below are minimal examples demonstrating how to use Nullstrap-DESeq2 and Nullstrap-edgeR for FDR calibration in DE analysis. For illustration, we use the synthetic dataset included with the package:

``` r
data("example_data", package = "NullstrapDE")
counts   <- example_data$counts     # gene count matrix
colData  <- example_data$colData    # condition labels
true_de  <- example_data$true_de    # true DE indicator (logical)
fdrcutoff <- 0.1
```

### Nullstrap-DESeq2

The following code is a quick example of how to implement the **Nullstrap-DESeq2** add-on.

``` r
  library(DESeq2)
  library(NullstrapDE)

  ### DESeq2
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = colData,
    design = ~ condition
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  res_deseq2 <- DESeq2::results(dds)
  message("FDR (DESeq2): ", fdr_deseq2)

  # Evaluate FDR for DESeq2
  sig_idx <- which(res_deseq2$padj < fdrcutoff & !is.na(res_deseq2$padj))
  tp <- sum(true_de[sig_idx])
  fp <- length(sig_idx) - tp
  fdr_deseq2 <- if (length(sig_idx) > 0) fp / length(sig_idx) else NA

  ### Nullstrap-DESeq2
  res_deseq2_nullstrap <- Nullstrap_DESeq2(dds = dds,
                                           fdrcutoff = fdrcutoff)
  # Evaluate FDR for Nullstrap-DESeq2
  true_de_select <- intersect(rownames(res_deseq2_nullstrap), rownames(dds)[true_de])
  tp <- length(true_de_select)
  fp <- nrow(res_deseq2_nullstrap) - tp
  fdr_deseq2_nullstrap <- if (nrow(res_deseq2_nullstrap) > 0) fp / nrow(res_deseq2_nullstrap) else NA
  message("FDR (Nullstrap-DESeq2): ", fdr_deseq2_nullstrap)
```

### Nullstrap-edgeR

``` r
  library(edgeR)
  library(NullstrapDE)
  
  ### edgeR
  group <- example_data$colData$condition
  design <- model.matrix(~group)
  dge <- DGEList(counts = example_data$counts, group = group)
  keep <- filterByExpr(dge, design = design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  res_fit <- glmQLFTest(fit)
  res_edger <- topTags(res_fit, n = Inf)$table

  # Evaluate FDR for edgeR
  sig_idx <- which(res_edger$FDR < fdrcutoff)
  sig_genes <- rownames(res_edger)[sig_idx]
  names(true_de) <- rownames(counts)
  true_de_filtered <- true_de[rownames(dge)]
  tp <- sum(true_de_filtered[sig_genes], na.rm = TRUE)
  fp <- length(sig_genes) - tp
  fdr_edger <- if (length(sig_genes) > 0) fp / length(sig_genes) else NA
  message("FDR (edgeR): ", fdr_edger)

  ### Nullstrap-edgeR
  res_edger_nullstrap <- Nullstrap_edgeR(dge = dge,
                                         design = design,
                                         fdrcutoff = fdrcutoff)
  true_de_select <- intersect(rownames(res_edger_nullstrap), rownames(counts)[true_de])
  tp <- length(true_de_select)
  fp <- nrow(res_edger_nullstrap) - tp
  fdr_edger_nullstrap <- if (nrow(res_edger_nullstrap) > 0) fp / nrow(res_edger_nullstrap) else NA
  message("FDR (Nullstrap-edgeR): ", fdr_edger_nullstrap)
```

## Related Manuscript

[Jiang, C., Wang, C., Li, J. J. (2025). Nullstrap-DE: A General Framework for Calibrating FDR and Preserving Power in DE Methods, with Applications to DESeq2 and edgeR. *arXiv*, 2025-07.](https://arxiv.org/abs/2507.20598)
