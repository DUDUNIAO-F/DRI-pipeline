#' Calculate the DRI (tRNA Regulatory Index) Score
#' @param expr_mat Gene expression matrix (rows = genes, columns = samples)
#' @param tds_genes Vector of genes used to calculate the differentiation score (e.g., c("TG","TSHR","PAX8","SLC5A5","NKX2-1"))
#' @param anchor_gene The anchor gene for correlation (default = "NSUN2")
#' @param min_nsun2_rho Minimum positive correlation threshold with the anchor gene (default = 0.15)
#' @param max_tds_rho Maximum negative correlation threshold with TDS (default = -0.10)
#' @return A list containing DRI scores, core gene set, and TDS scores
calculate_DRI <- function(expr_mat, tds_genes, anchor_gene = "NSUN2", 
                          min_nsun2_rho = 0.15, max_tds_rho = -0.10) {
  
  # 1. Ensure genes exist in the matrix
  if (!anchor_gene %in% rownames(expr_mat)) {
    stop("Anchor gene not found in the expression matrix.")
  }
  
  tds_genes <- intersect(tds_genes, rownames(expr_mat))
  if (length(tds_genes) < 3) {
    stop("Too few TDS genes found in the expression matrix. At least 3 are required.")
  }
  
  # 2. Calculate TDS and correlations
  tds <- colMeans(expr_mat[tds_genes, , drop = FALSE], na.rm = TRUE)
  
  # Calculate Spearman correlation for each gene against anchor_gene and TDS
  nsun2_cor <- apply(expr_mat, 1, function(x) cor(x, expr_mat[anchor_gene, ], method = "spearman"))
  tds_cor <- apply(expr_mat, 1, function(x) cor(x, tds, method = "spearman"))
  
  # 3. Filter gene sets
  core_genes <- names(which(nsun2_cor > min_nsun2_rho & tds_cor < max_tds_rho))
  supp_genes <- setdiff(names(which(tds_cor < -0.15)), core_genes)
  dri_genes <- union(core_genes, supp_genes)
  
  if (length(dri_genes) == 0) {
    stop("No genes passed the correlation filter. Please check your data or adjust thresholds.")
  }
  
  # 4. Calculate DRI score
  # Z-score standardization across samples for the selected genes
  z_expr <- scale(t(expr_mat[dri_genes, , drop = FALSE]))
  dri_score <- colMeans(z_expr, na.rm = TRUE)
  
  return(list(
    DRI_score = dri_score,
    TDS_score = tds,
    DRI_genes = dri_genes,
    N_genes = length(dri_genes)
  ))
}