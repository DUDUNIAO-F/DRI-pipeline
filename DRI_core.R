#' Calculate the DRI (tRNA Regulatory Index) Score
#' @param expr_mat Gene expression matrix (rows = genes, columns = samples)
#' @param tds_genes Vector of genes used to calculate the differentiation score
#' @param anchor_gene The anchor gene for correlation (default = "NSUN2")
#' @param min_nsun2_rho Minimum positive correlation threshold with the anchor gene
#' @param max_tds_rho Maximum negative correlation threshold with TDS
#' @return A list containing DRI scores, core gene set, and TDS scores
calculate_DRI <- function(expr_mat, tds_genes, anchor_gene = "NSUN2",
                          min_nsun2_rho = 0.15, max_tds_rho = -0.10) {

  expr_mat <- as.matrix(expr_mat)

  # 1. Ensure required genes exist
  if (!anchor_gene %in% rownames(expr_mat)) {
    stop("Anchor gene not found in the expression matrix.")
  }
  tds_genes <- intersect(tds_genes, rownames(expr_mat))
  if (length(tds_genes) < 3) {
    stop("Too few TDS genes found in the expression matrix. At least 3 are required.")
  }

  # 2. TDS and per-gene correlations
  tds <- colMeans(expr_mat[tds_genes, , drop = FALSE], na.rm = TRUE)
  anchor_vec <- expr_mat[anchor_gene, ]

  nsun2_cor <- apply(expr_mat, 1, function(x)
    cor(x, anchor_vec, method = "spearman", use = "pairwise.complete.obs"))
  tds_cor <- apply(expr_mat, 1, function(x)
    cor(x, tds, method = "spearman", use = "pairwise.complete.obs"))

  # 3. Gene set selection
  core_genes <- names(which(nsun2_cor > min_nsun2_rho & tds_cor < max_tds_rho))
  supp_genes <- setdiff(names(which(tds_cor < -0.15)), core_genes)
  dri_genes <- union(core_genes, supp_genes)
  if (length(dri_genes) == 0) {
    stop("No genes passed the correlation filter. Please check your data or adjust thresholds.")
  }

  # 4. Per-sample DRI score
  sub <- expr_mat[dri_genes, , drop = FALSE]
  v <- apply(sub, 1, sd, na.rm = TRUE)
  sub <- sub[is.finite(v) & v > 0, , drop = FALSE]
  if (nrow(sub) == 0) stop("All selected DRI genes have zero variance.")

  z_expr <- scale(t(sub))                      # samples x genes
  dri_score <- rowMeans(z_expr, na.rm = TRUE)  # per-sample score
  names(dri_score) <- rownames(z_expr)

  return(list(
    DRI_score = dri_score,
    TDS_score = tds,
    DRI_genes = rownames(sub),
    N_genes = nrow(sub)
  ))
}
