#' NullstrapDE calibration for DESeq2
#'
#' @param dds A DESeqDataSet (optional)
#' @param counts Count matrix (if dds not supplied)
#' @param colData Sample metadata (if dds not supplied)
#' @param fdrcutoff Target FDR level
#' @param correct Correction mode: "none", "ratio", "half"
#' @param stat Statistic for thresholding: "fc" or "pval"
#' @param sizeFactors_sample Whether to resample size factors
#' @param conservative Conservative mode for FDP estimation
#' @param seed Random seed
#'
#' @return A DESeq2 results table of significant genes
#' @export
Nullstrap_DESeq2 <- function(dds = NULL,
                             counts = NULL,
                             colData = NULL,
                             fdrcutoff = 0.05,
                             correct = "none",
                             stat = "fc",
                             sizeFactors_sample = TRUE,
                             conservative = FALSE,
                             seed = 5){

  #------------------------------
  # Step 0. Prepare DESeqDataSet
  #------------------------------
  if (is.null(dds)) {
    if (is.null(counts) || is.null(colData)) {
      stop("Provide either a DESeqDataSet (dds) or both counts and colData.")
    }
    dds <- suppressMessages(
      DESeqDataSetFromMatrix(countData = counts,
                             colData = colData,
                             design = ~ condition)
    )
  }
  if (!"mu" %in% assayNames(dds)) {
    message("running DESeq2")
    dds <- DESeq(dds, quiet = TRUE)
  }

  count_dim <- dim(dds@assays@data@listData$counts)
  if(is.null(colnames(dds))){
    colnames(dds) <- as.character(1:count_dim[2])
  }
  if(is.null(rownames(dds))){
    rownames(dds) <- as.character(1:count_dim[1])
  }

  res <- results(dds)

  ##### Run Nullstrap-DESeq2

  #-------------------------------------------------
  # Step 1. Extract Model Parameters from DESeq2 Fit
  #-------------------------------------------------

  # message("running Nullstrap-DESeq2")

  # Extract coefs: beta or pvals
  coefs <- coef(dds)
  coefs_pval <- -log(res$pvalue)

  # Extract nuisance parameters
  disp <- dispersions(dds)
  sf <- sizeFactors(dds)

  # Remove genes with NA/Inf in coefs or pvals
  keep <- which((complete.cases(coefs)) & (!is.na(coefs_pval)) & (!is.infinite(coefs_pval)))
  keep_disp <- which(!is.na(disp))
  keep <- intersect(keep, keep_disp)
  coefs <- coefs[keep, ]
  coefs_pval <- coefs_pval[keep]
  disp <- disp[keep]

  # Subset coefs columns
  coef_beta0 <- coefs[, 1]
  coef_beta1 <- coefs[, ncol(coefs)]

  if(ncol(coefs)>2){
    # Exist other covariates
    coef_others <- coefs[, -c(1, ncol(coefs))]
  }else{
    coef_others <- rep(0, nrow(coefs))
  }


  #---------------------------------------------
  # Step 2. Simulate Synthetic Null Data (β = 0)
  #---------------------------------------------
  n_genes_keep <- length(keep)
  n_samples <- dim(dds)[2]

  if(sizeFactors_sample){
    sf_range <- quantile(sf, c(0.1, 0.9))
    sf_in_range <- sf[sf >= sf_range[1] & sf <= sf_range[2]]
    set.seed(seed)
    sf_in_range_sampled <- sample(sf_in_range, length(sf), replace = TRUE)
  }else{
    sf_in_range_sampled <- sf
  }


  countData_syn <- matrix(nrow = n_genes_keep, ncol = n_samples)
  for (i in 1:n_genes_keep) {
    baseline_mu <- 2^(coef_beta0[i] + coef_others[i])
    for (j in 1:n_samples) {
      expected_mu <- baseline_mu * sf_in_range_sampled[j]
      size_param <- 1 / disp[i]      # Negative binomial size parameter
      countData_syn[i, j] <- rnbinom(1, mu = expected_mu, size = size_param)
    }
  }
  rownames(countData_syn) <- rownames(dds)[keep]
  colnames(countData_syn) <- colnames(dds)

  #----------------------------------------------
  # Step 3. Fit DESeq2 to the Synthetic Null Data
  #----------------------------------------------
  dds_syn <- suppressMessages(
    DESeqDataSetFromMatrix(countData = countData_syn,
                           colData = dds@colData,
                           design = dds@design)
  )
  # Override nuisance parameters with the previously estimated values
  sizeFactors(dds_syn) <- sf
  #sizeFactors(dds_syn) <- estimateSizeFactors(dds_syn)
  #sizeFactors(dds_syn) <- rep(1, length(sf))
  dispersions(dds_syn) <- disp

  # Run the Wald test
  dds_syn <- nbinomWaldTest(dds_syn)
  res_syn <- results(dds_syn)

  # Extract coefs from synthetic null
  coefs_syn <- coef(dds_syn)
  coefs_syn_pval <- -log(res_syn$pvalue)
  coefs_se_syn <- res_syn$lfcSE


  #------------------------------------
  # Step 4. Thresholding
  #------------------------------------
  if(stat == "pval"){
    # Use p values
    keep_syn <- which(!is.na(coefs_syn_pval))
    coef_real <- abs(coefs_pval[keep_syn])
    coef_syn  <- abs(coefs_syn_pval[keep_syn])
  }else if(stat == "fc"){
    # Use Wald stat
    keep_syn <- which(complete.cases(coefs_syn) & !is.na(coefs_se_syn))
    coef_real <- abs(coefs[keep_syn, 2]/coefs_se_syn[keep_syn])
    coef_syn  <- abs(coefs_syn[keep_syn, 2]/coefs_se_syn[keep_syn])
  }


  #correct_factor <- sqrt((n_samples-1)/n_samples)
  if(correct == "ratio"){
    correct_factor <- 1 / (1 + sqrt(log2(length(coef_real)) / n_samples))
  }else if(correct == "half"){
    correct_factor <- 1/2
  }else if(correct == "none"){
    correct_factor <- 1
  }

  threshold <- binary_search(coef_real, coef_syn,
                             q = correct_factor*fdrcutoff,
                             conservative = conservative)

  # Identify rejections
  sig_idx_deseq2_syn <- which(coef_real >= threshold)
  if (max(coef_real) < max(coef_syn) ) {
    sig_idx_deseq2_syn <- NULL  # safeguard
  }

  # Final results
  keep_index = rownames(dds)[keep][keep_syn]
  sig_index = keep_index[sig_idx_deseq2_syn]
  res_nullstrap <- res[sig_index,]
  res_nullstrap
}
