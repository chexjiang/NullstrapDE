test_that("Nullstrap_edgeR runs without error", {
  data("example_data", package = "NullstrapDE")

  fdrcutoff <- 0.1

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

})
