test_that("Nullstrap_DESeq2 runs without error", {
  data("example_data", package = "NullstrapDE")

  fdrcutoff <- 0.1

  ### DESeq2
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = example_data$counts,
    colData = example_data$colData,
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

})
