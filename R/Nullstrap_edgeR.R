
Nullstrap_edgeR <- function(counts = NULL,
                            colData = NULL,
                            dge = NULL,
                            design = NULL,
                            fdrcutoff = 0.05,
                            correct = "none",
                            stat = "fc",
                            test = "lrt",
                            sizeFactors_sample = TRUE,
                            conservative = FALSE,
                            seed = 5) {
  set.seed(seed)


  #######################################
  # Step 0. Prepare DGEList
  #######################################
  # if (!is.null(dge)) {
  #   # User provided pre-built DGEList
  #   if (is.null(dge$samples$group)) {
  #     stop("dge must include sample group info in dge$samples$group")
  #   }
  #   group <- dge$samples$group
  #   design <- model.matrix(~group)
  #   original_counts <- dge$counts
  #   if (is.null(dge$tagwise.dispersion)) {
  #     message("running edgeR")
  #     dge <- calcNormFactors(dge)
  #     dge <- estimateDisp(dge, design)
  #   }
  # } else {
  #   # Construct DGEList from counts and colData
  #   if (is.null(counts) || is.null(colData)) {
  #     stop("Provide either a DGEList (dge) or both counts and colData.")
  #   }
  #   message("running edgeR")
  #   group <- colData$condition
  #   design <- model.matrix(~group)
  #   dge <- DGEList(counts = counts, group = group)
  #   keep <- filterByExpr(dge, group = group)
  #   dge <- dge[keep, , keep.lib.sizes = FALSE]
  #   dge <- calcNormFactors(dge)
  #   dge <- estimateDisp(dge, design)
  #   original_counts <- counts
  # }

  if (!is.null(dge)) {
    if (is.null(dge$samples$group)) {
      stop("dge must include sample group info in dge$samples$group")
    }
    group <- dge$samples$group
    original_counts <- dge$counts

    # Build design from provided colData or dge$samples
    if (is.null(design)) {
      stop("design is not provided.")
    }

    if (is.null(dge$tagwise.dispersion)) {
      message("running edgeR")
      dge <- calcNormFactors(dge)
      dge <- estimateDisp(dge, design)
    }

  } else {
    if (is.null(counts) || is.null(colData)) {
      stop("Provide either a DGEList (dge) or both counts and colData.")
    }

    original_counts <- counts

    # If no design is provided, use all variables in colData
    if (is.null(design)) {
      design <- model.matrix(~ ., data = colData)
    }
    message("running edgeR")
    dge <- DGEList(counts = counts)
    keep <- filterByExpr(dge, design = design)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- calcNormFactors(dge)
    dge <- estimateDisp(dge, design)
    group <- dge$samples$group
  }


  #######################################
  # Step 1. Fit edgeR model to real data
  #######################################

  #message("running Nullstrap-edgeR")
  if(test == "lrt"){
    fit <- glmFit(dge, design)
    res_fit <- glmLRT(fit)
  }else if(test == "qlf"){
    fit <- glmQLFit(dge, design)
    res_fit <- glmQLFTest(fit)
  }

  res <- res_fit$table
  pval <- res$PValue
  disp <- dge$trended.dispersion
  logFC <- res$logFC
  coefs <- fit$coefficients

  # Subset coefs columns
  coef_beta0 <- coefs[, 1]
  coef_beta1 <- coefs[, ncol(coefs)]

  if(ncol(coefs)>2){
    # Exist other covariates
    coef_others <- coefs[, -c(1, ncol(coefs))]
  }else{
    coef_others <- rep(0, nrow(coefs))
  }
  log_Ni <- log(dge$samples$lib.size * dge$samples$norm.factors)


  ### Keep only complete cases
  valid <- which(!is.na(pval) & is.finite(pval) & (pval != 0))
  res <- res[valid, ]
  pval <- pval[valid]
  disp <- disp[valid]
  logFC <- logFC[valid]
  coef_beta0 <- coef_beta0[valid]
  coef_beta1 <- coef_beta1[valid]
  coef_others <- coef_others[valid]
  gene_names <- rownames(res)


  if(sizeFactors_sample){
    log_Ni_range <- quantile(log_Ni, c(0.1, 0.9))
    log_Ni_in_range <- log_Ni[log_Ni >= log_Ni_range[1] & log_Ni <= log_Ni_range[2]]
    set.seed(seed)
    log_Ni_in_range_sampled <- sample(log_Ni_in_range, length(log_Ni), replace = TRUE)
  }else{
    log_Ni_in_range_sampled <- log_Ni
  }


  mu <- matrix(nrow = length(coef_beta0), ncol = ncol(dge))
  for (j in 1:ncol(dge)) {
    mu[, j] <- exp(coef_beta0 + coef_others + log_Ni_in_range_sampled[j])
  }

  #######################################
  # Step 2. Simulate null counts
  #######################################
  n_genes <- length(valid)
  n_samples <- ncol(dge)
  syn_counts <- matrix(nrow = n_genes, ncol = n_samples)

  set.seed(seed)
  for (i in 1:n_genes) {
    for (j in 1:n_samples) {
      size_param <- 1 / disp[i]
      syn_counts[i, j] <- rnbinom(1, mu = mu[i, j], size = size_param)
    }
  }
  rownames(syn_counts) <- gene_names
  colnames(syn_counts) <- colnames(counts)

  #######################################
  # Step 3. Fit edgeR to synthetic null
  #######################################
  dge_syn <- DGEList(counts = syn_counts, group = group)
  dge_syn$samples$norm.factors <- dge$samples$norm.factors

  if(test == "lrt"){
    fit_syn <- glmFit(y=dge_syn, design=design, dispersion=disp)
    res_fit_syn <- glmLRT(fit_syn)
  }else if(test == "qlf"){
    fit_syn <- glmQLFit(y=dge_syn, design=design, dispersion=disp)
    res_fit_syn <- glmQLFTest(fit_syn)
  }


  res_syn <- res_fit_syn$table
  pval_syn <- res_syn$PValue
  logFC_syn <- res_syn$logFC

  if(stat == "pval"){
    stat_real <- -log(pval)
    stat_syn <- -log(pval_syn)
  }else if(stat == "fc") {
    stat_real <- abs(logFC)
    stat_syn <- abs(logFC_syn)
  }


  #######################################
  # Step 4. Thresholding
  #######################################
  #correct_factor <- 1 / (1 + sqrt(log2(nrow(original_counts)) / n_samples))
  if(correct == "ratio"){
    correct_factor <- 1 / (1 + sqrt(log2(nrow(original_counts)) / n_samples))
  }else if(correct == "half"){
    correct_factor <- 1/2
  }else if(correct == "none"){
    correct_factor <- 1
  }
  # correct_factor <- 1
  threshold <- binary_search(stat_real, stat_syn,
                             q = correct_factor * fdrcutoff,
                             conservative = conservative)

  sig_idx <- which(stat_real >= threshold)
  if (max(stat_real) < max(stat_syn)) {
    sig_idx <- NULL
  }

  sig_genes <- gene_names[sig_idx]
  res_nullstrap <- res[sig_genes, ]
  return(res_nullstrap)
}









